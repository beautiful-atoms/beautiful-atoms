"""
"""
import bpy
import numpy as np
from time import time

class Isosurfacesetting():
    """
    Isosurfacesetting object

    The Isosurfacesetting object store the isosurface information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    def __init__(self, label, volume) -> None:
        self.label = label
        self.name = 'bisosurface'
        self.volume = volume
        self[1] = [0.002, [1, 1, 0, 0.5]]
    def draw_volume(self, volume):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        # remove old volume point
        if volume is None: return
        name = "volume_%s"%self.label
        if name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[name], du_unlink = True)
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
    @property
    def collection(self):
        return self.get_collection()
    def get_collection(self):
        collection = getattr(bpy.data.collections[self.label], self.name)
        return collection
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
        mesh = bpy.data.objects["volume_%s"%self.label].data
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
        n = len(self.mesh.vertices)
        volume = np.empty(n*3, dtype=np.float64)
        self.mesh.vertices.foreach_get('co', volume)  
        volume = volume.reshape(-1, 1)
        volume = volume[:self.npoint]
        volume = volume.reshape(self.shape)
        return volume
    @volume.setter
    def volume(self, volume):
        self.draw_volume(volume)
    def __setitem__(self, index, value):
        """
        Add isosurface one by one
        """
        p = self.find(index)
        if p is None:
            p = self.collection.add()
        p.name = str(index)
        if isinstance(value, (int, float)):
            value = [value]
        p.level = value[0]
        if len(value) == 2:
            p.color = value[1]
    def set_default(self):
        """
        """
        for sp, data in self.species.items():
            self[sp] = [np.append(data['color'][:3], 0.3), 0.005]
    def add_isosurfaces(self, isosurfacepair):
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
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    def __len__(self):
        return len(self.collection)
    def find(self, name):
        i = self.collection.find(str(name))
        if i == -1:
            return None
        else:
            return self.collection[i]
    def build_isosurface(self, cell):
        volume = self.volume
        isosurface = []
        for iso in self.collection:
            level = iso.level
            color = iso.color
            verts, faces = calc_isosurface(volume, cell, level)
            isosurface.append((verts, faces, color))
        return isosurface


def calc_isosurface(volume, cell, level,
                    gradient_direction = 'descent',
                    step_size = 1):
    """
    
    Computes an isosurface from a volume grid.
    
    Parameters:
    
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
