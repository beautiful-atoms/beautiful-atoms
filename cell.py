"""
"""
import bpy
import numpy as np
from ase.cell import Cell
from batoms.butils import object_mode

class Bcell():
    """
    Unit cell of three dimensions.

    """
    def __init__(self, label, 
                array = np.zeros([3, 3]), 
                location = np.array([0, 0, 0]),
                color = (0.0, 0.0, 0.0, 1.0),
                material_style = 'default',
                bsdf_inputs = None,
                ) -> None:
        """
        ver: 3x3 verlike object
          The three cell vectors: cell[0], cell[1], and cell[2].
        """
        self.label = label
        self.name = 'cell_%s_edge'%(self.label)
        self.edges = [[3, 0], [3, 1], [4, 0], [4, 1],
                    [2, 5], [2, 6], [7, 5], [7, 6], 
                    [3, 2], [0, 6], [1, 5], [4, 7]
            ]
        self.width = 0.05
        self.color = color
        self.draw_cell_edge(array, location)
    def draw_cell_edge(self, array, location):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        if array is None:
            array = np.zeros([3, 3])
        if not len(array) == 8:
            cell = Cell.new(array)
            verts = self.array2verts(cell.array)
        else:
            verts = array - array[3]
            location = array[3]
        if self.name not in bpy.data.objects:
            
            mesh = bpy.data.meshes.new("cell_%s_edge"%self.label)
            mesh.from_pydata(verts, self.edges, [])  
            mesh.update()
            for f in mesh.polygons:
                f.use_smooth = True
            obj_edge = bpy.data.objects.new("cell_%s_edge"%self.label, mesh)
            obj_edge.data = mesh
            obj_edge.location = location
            obj_edge.bcell.flag = True
            bpy.data.collections['Collection'].objects.link(obj_edge)
        elif bpy.data.objects[self.name].bcell.flag:
            # print('%s exist and is bcell, use it.'%self.name)
            pass
        else:
            raise Exception("Failed, the name %s already \
                in use and is not Bcell object!"%self.name)
        bpy.context.view_layer.update()
    def build_cell_cylinder(self):
         #
        cell_cylinder = {'lengths': [], 
                      'centers': [],
                      'normals': [],
                      'vertices': 16,
                      'width': self.width,
                      'color': self.color,
                      'battr_inputs': {},
                      }
        if np.max(abs(self.verts)) < 1e-6:
            return cell_cylinder
        for e in self.edges:
            center = (self.verts[e[0]] + self.verts[e[1]])/2.0
            vec = self.verts[e[0]] - self.verts[e[1]]
            length = np.linalg.norm(vec)
            nvec = vec/length
            cell_cylinder['lengths'].append(length)
            cell_cylinder['centers'].append(center)
            cell_cylinder['normals'].append(nvec)
        return cell_cylinder
    def __repr__(self) -> str:
        numbers = self.array.tolist()
        s = 'Cell({})'.format(numbers)
        return s
    def __getitem__(self, index):
        return self.array[index]
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        """Set unit cell vectors.

        Parameters:

        Examples:

        """
        bcell = self.bcell
        array = self.array
        array[index] = value
        verts = self.array2verts(array)
        for i in range(8):
            bcell.data.vertices[i].co = np.array(verts[i])
    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.array
    @property
    def bcell(self):
        return self.get_bcell()
    def get_bcell(self):
        return bpy.data.objects['cell_%s_edge'%(self.label)]
    @property
    def array(self):
        return self.get_array()
    def get_array(self):
        cell = np.array([self.verts[0] - self.verts[3],
                         self.verts[1] - self.verts[3],
                         self.verts[2] - self.verts[3]])
        return cell
    @property
    def local_verts(self):
        return self.get_local_verts()
    def get_local_verts(self):
        bcell = self.bcell
        return np.array([bcell.data.vertices[i].co for i in range(8)])
    @property
    def verts(self):
        return self.get_verts()
    def get_verts(self):
        return np.array([self.bcell.matrix_world @ \
                self.bcell.data.vertices[i].co for i in range(8)])
    @property
    def origin(self):
        return self.verts[3]
    def array2verts(self, array):
        """
        """
        verts = np.array([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, 0],
            [1, 1, 0],
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 1],
            ])
        verts = np.dot(verts, array)
        return verts
    @property
    def location(self):
        return self.get_location()
    def get_location(self):
        return np.array(self.bcell.location)
    def copy(self, label):
        object_mode()
        cell = Bcell(label, array = self.array, location = self.bcell.location)
        return cell
    def repeat(self, m):
        self[:] = np.array([m[c] * self.array[c] for c in range(3)])
    