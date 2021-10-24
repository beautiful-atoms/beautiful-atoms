"""
"""
import bpy
import numpy as np
from time import time
from batoms.bondsetting import Setting

default_colors = [(1, 1, 0, 0.8), (0.0, 0.0, 1.0, 0.8)]

class IsosurfaceSetting(Setting):
    """
    IsosurfaceSetting object

    The IsosurfaceSetting object store the isosurface information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    def __init__(self, label, volume = None, isosurfacesetting = None) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bisosurface'
        # add a default level
        if len(self) == 0:
            self['1'] = {'level': 0.002, 'color': [1, 1, 0, 0.8]}
        if isosurfacesetting is not None:
            for key, data in isosurfacesetting.items():
                self[key] = data
        if volume is not None:
            self.volume = volume
    def draw_volume(self, volume):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        # remove old volume point
        tstart = time()
        if volume is None: return
        name = "volume_%s"%self.label
        if name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[name], do_unlink = True)
        shape = volume.shape
        volume = volume.reshape(-1, 1)
        npoint = len(volume)
        dn = 3 - npoint % 3
        verts = np.append(volume, np.zeros((dn, 1)), axis = 0)
        verts = verts.reshape(-1, 3)
        mesh = bpy.data.meshes.new("mesh_%s_volume"%self.label)
        mesh.from_pydata(verts, [], [])  
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        obj.data = mesh
        obj.bvolume.is_bvolume = True
        obj.bvolume.shape = shape
        obj.bvolume.npoint = npoint
        bpy.data.collections[self.label].children['%s_volume'%self.label].objects.link(obj)
        obj.hide_set(True)
        obj.hide_render = True
        print('Draw volume: {0:1.2f}'.format(time() - tstart))
    @property
    def npoint(self):
        return self.get_npoint()
    def get_npoint(self):
        if "volume_%s"%self.label not in bpy.data.objects:
            return 0
        return bpy.data.objects["volume_%s"%self.label].bvolume.npoint
    @npoint.setter
    def npoint(self, npoint):
        self.set_npoint(npoint)
    def set_npoint(self, npoint):
        bpy.data.objects["volume_%s"%self.label].bvolume.npoint = npoint
    @property
    def mesh(self):
        return self.get_mesh()
    def get_mesh(self):
        name = "volume_%s"%self.label
        if name not in bpy.data.objects:
            return None
        mesh = bpy.data.objects[name].data
        return mesh
    @property
    def shape(self):
        return self.get_shape()
    def get_shape(self):
        shape = bpy.data.objects["volume_%s"%self.label].bvolume.shape
        return shape
    @property
    def volume(self):
        return self.get_volume()
    def get_volume(self):
        tstart = time()
        if self.mesh is None:
            return None
        n = len(self.mesh.vertices)
        volume = np.empty(n*3, dtype=np.float64)
        self.mesh.vertices.foreach_get('co', volume)  
        volume = volume.reshape(-1, 1)
        volume = volume[:self.npoint]
        volume = volume.reshape(self.shape)
        # print('Read volume: {0:1.2f}'.format(time() - tstart))
        return volume
    @volume.setter
    def volume(self, volume):
        self.draw_volume(volume)
    def set_default(self):
        """
        """
        for sp, data in self.species.items():
            self[sp] = [np.append(data['color'][:3], 0.3), 0.005]
    def add(self, isosurfacepair):
        for key in isosurfacepair:
            self.set_default(key)
    def remove_isosurfaces(self, isosurfacepair):
        for key in isosurfacepair:
            name = '%s-%s'%(key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Center     level                color  \n'
        for iso in self.collection:
            s += '{0:10s}   {1:1.6f}  [{2:1.2f}  {3:1.2f}  {4:1.2f}  {5:1.2f}] \n'.format(\
                iso.name, iso.level, iso.color[0], iso.color[1], iso.color[2], iso.color[3])
        s += '-'*60 + '\n'
        return s
    def build_isosurface(self, cell):
        volume = self.volume
        isosurface = []
        for iso in self.collection:
            name = iso.name
            level = iso.level
            color = iso.color
            verts, faces = calc_isosurface(volume, cell, level)
            isosurface.append((name, verts, faces, color))
        return isosurface


def calc_isosurface(volume, cell, level,
                    gradient_direction = 'descent',
                    step_size = 1):
    """
    
    Computes an isosurface from a volume grid.
    
    Parameters:

    volume: np.array

    cell: np.array

    level: float
    
    """
    from batoms.tools import get_cell_vertices
    from skimage import measure

    cell_vertices = get_cell_vertices(cell)
    cell_vertices.shape = (2, 2, 2, 3)
    cell_origin = cell_vertices[0,0,0]
    #
    spacing = tuple(1.0/np.array(volume.shape))
    mlevel = np.mean(volume)
    if not level:
        level = mlevel*10
    # print('iso level: {0:1.9f}, iso mean: {1:1.9f}'.format(level, mlevel))
    scaled_verts, faces, normals, values = measure.marching_cubes(volume, level = level,
                    spacing=spacing,gradient_direction=gradient_direction , 
                    allow_degenerate = False, step_size=step_size)
    scaled_verts = scaled_verts.dot(cell)
    scaled_verts -= cell_origin
    faces = list(faces)
    return scaled_verts, faces
