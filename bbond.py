"""Definition of the Bbond class.

This module defines the Bbond object in the bbonds package.

"""

from time import time
import bpy
import bmesh
from batoms.butils import object_mode
from batoms.material import material_styles_dict
import numpy as np
from batoms.base import BaseObject

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE"]


class Bbond(BaseObject):
    """Bbond Class
    
    A Bbond object is linked to this main collection in Blender. 

    Parameters:

    label: str
        Name of the Bbonds.
    species: str
        species of the atoms.
    positions: array
        positions
    locations: array
        The object’s origin location in global coordinates.
    element: str or list
        element of the atoms, list for fractional Occupancy
    segments: list of 2 Int
        Number of segments used to draw the UV_Sphere
        Default: [32, 16]
    subdivisions: Int
        Number of subdivision used to draw the ICO_Sphere
        Default: 2
    color_style: str
        "JMOL", "ASE", "VESTA"
    radii_style: str
        "covelent", "vdw", "ionic"
    shape: Int
        0, 1, or 2. ["UV_SPHERE", "ICO_SPHERE", "CUBE"]

    Examples:

    >>> from batoms.bbond import Bbond
    >>> c = Bbond('C', [[0, 0, 0], [1.2, 0, 0]])

    """
    

    def __init__(self, 
                label = None,
                species = None,
                positions = None,
                width = 0.1,
                location = np.array([0, 0, 0]),
                segments = 16,
                color = (1, 0, 0, 1),
                material_style = 'default',
                material = None,
                node_inputs = None,
                battr_inputs = {},
                 ):
        #
        if species is not None:
            self.label = label
            self.species = species
            self.name = species
            obj_name = '%s_bond_%s'%(self.label, self.species)
        else:
            obj_name = label
        BaseObject.__init__(self, obj_name = obj_name)
        if positions is not None:
            positions = np.array(positions)
            if len(positions.shape) == 2:
                self._frames = np.array([positions])
            elif len(positions.shape) == 3:
                self._frames = positions
                positions = self._frames[0]
            else:
                raise Exception('Shape of positions is wrong!')
            self.set_material(color, node_inputs, material_style, material)
            self.set_object(positions, segments = segments, width = width, 
                            location = location, battr_inputs = battr_inputs)
            self.set_frames(self._frames, only_basis = True)
        else:
            self.from_bbond(label)
    def set_material(self, color, node_inputs = None, material_style = 'default', material = None):
        """
        """
        from batoms.material import create_material
        name = 'material_bond_{0}_{1}'.format(self.label, self.species)
        if material:
            material = material.copy()
            material.name = name
        elif name not in bpy.data.materials:
            material = create_material(name,
                        color,
                        node_inputs = node_inputs,
                        material_style = material_style,
                        backface_culling = True)
    def set_object(self, positions, width = 0.1, segments = 16,
                        location = (0, 0, 0), battr_inputs = {}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        from batoms.source_data import bond_source
        from batoms.bdraw import cylinder_mesh_from_vec
        tstart = time()
        mesh_data = bond_source[segments]
        vertices, faces = cylinder_mesh_from_vec(positions[:, 0:3], 
                            positions[:, 3:6], positions[:, 6:7], 
                            width, mesh_data)
        if self.obj_name in bpy.data.objects:
            obj = bpy.data.objects.get(self.obj_name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(self.obj_name)
        mesh.from_pydata(vertices, [], faces)
        mesh.update()
        mesh.polygons.foreach_set('use_smooth', [True]*len(mesh.polygons))
        obj = bpy.data.objects.new(self.obj_name, mesh)
        obj.location = location
        obj.data.materials.append(self.material)
        for name, inputs in battr_inputs.items():
            battr = getattr(obj, name)
            for key, value in inputs.items():
                setattr(battr, key, value)
        bpy.data.collections['Collection'].objects.link(obj)
        bpy.context.view_layer.update()
        # print('draw bond: {0:10.2f} s'.format(time() - tstart))
    def from_bbond(self, label):
        if label not in bpy.data.objects:
            raise Exception("%s is not a object!"%label)
        elif not bpy.data.objects[label].bbond.flag:
            raise Exception("%s is not Bbond object!"%label)
        obj = bpy.data.objects[label]
        self.species = obj.bbond.species
        self.label = obj.bbond.label
        self.element = obj.bbond.element
        self.species_data = {
            'radius':obj.bbond.radius,
            'scale':obj.scale,
        } 
    @property
    def mesh(self):
        return self.get_mesh()
    def get_mesh(self):
        return self.obj.children[0]
    @property
    def material(self):
        return self.get_material()
    def get_material(self):
        return bpy.data.materials['material_bond_%s_%s'%(self.label, self.species)]
    @property
    def scale(self):
        return self.get_scale()
    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)
    def get_scale(self):
        return np.array(self.mesh.scale)
    def set_scale(self, scale):
        if isinstance(scale, float) or isinstance(scale, int):
            scale = [scale]*3
        self.mesh.scale = scale
    @property
    def width(self):
        return self.get_width()
    def get_width(self):
        return np.array(self.obj.bbond.width)
    @property
    def segments(self):
        return self.get_segments()
    def get_segments(self):
        return self.obj.bbond.segments
    @segments.setter
    def segments(self, segments):
        self.set_segments(segments)
    def set_segments(self, segments):
        if not isinstance(segments, int):
            raise Exception('Segments should be int!')
        self.clean_bbonds_objects('mesh_bond_%s_%s'%(self.label, self.species))
        mesh = self.set_mesh(segments = segments)
        mesh.parent = self.obj
    @property
    def local_positions(self):
        return self.get_local_positions()
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        nvert = len(self.obj.data.vertices)
        local_positions = np.empty(nvert*3, dtype=np.float64)
        self.obj.data.vertices.foreach_get('co', local_positions)  
        local_positions = local_positions.reshape((nvert, 3))
        return local_positions
    @property
    def positions(self):
        return self.get_positions()
    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)
    def get_positions(self):
        """
        Get global positions.
        """
        from batoms.tools import local2global
        positions = local2global(self.local_positions, 
                np.array(self.obj.matrix_world))
        return positions
    def set_positions(self, positions):
        """
        Set global positions to local vertices
        """
        from batoms.tools import local2global
        natom = len(self)
        if len(positions) != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natom))
        positions = local2global(positions, 
                np.array(self.obj.matrix_world), reversed = True)
        # rashpe to (natoms*3, 1) and use forseach_set
        positions = positions.reshape((natom*3, 1))
        self.obj.data.vertices.foreach_set('co', positions)
        self.obj.data.update()
    def get_scaled_positions(self, cell):
        """
        Get array of scaled_positions.
        """
        from ase.cell import Cell
        cell = Cell.new(cell)
        scaled_positions = cell.scaled_positions(self.local_positions)
        return scaled_positions
    @property
    def nframe(self):
        return self.get_nframe()
    def get_nframe(self):
        if self.obj.data.shape_keys is None:
            return 0
        nframe = len(self.obj.data.shape_keys.key_blocks)
        return nframe
    @property
    def frames(self):
        return self.get_frames()
    @frames.setter
    def frames(self, frames):
        self.set_frames(frames)
    def get_frames(self):
        """
        read shape key
        """
        from batoms.tools import local2global
        obj = self.obj
        n = len(self)
        nframe = self.nframe
        frames = np.empty((nframe, n, 3), dtype=np.float64)
        for i in range(nframe):
            positions = np.empty(n*3, dtype=np.float64)
            sk = obj.data.shape_keys.key_blocks[i]
            sk.data.foreach_get('co', positions)
            local_positions = positions.reshape((n, 3))
            local_positions = local2global(local_positions, 
                            np.array(self.obj.matrix_world))
            frames[i] = local_positions
        return frames
    @property
    def color(self):
        return self.get_color()
    @color.setter
    def color(self, color):
        """
        >>> h.color = [0.8, 0.1, 0.3, 1.0]
        """
        self.set_color(color)
    def get_color(self):
        """

        """
        Viewpoint_color = self.material.diffuse_color
        for node in self.material.node_tree.nodes:
            if 'Base Color' in node.inputs:
                node_color = node.inputs['Base Color'].default_value[:]
            if 'Alpha' in node.inputs:
                Alpha = node.inputs['Alpha'].default_value
        color = [node_color[0], node_color[1], node_color[2], Alpha]
        return color
    def set_color(self, color):
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.material.diffuse_color = color
        for node in self.material.node_tree.nodes:
            if 'Base Color' in node.inputs:
                node.inputs['Base Color'].default_value = color
            if 'Alpha' in node.inputs:
                node.inputs['Alpha'].default_value = color[3]
    @property
    def node(self):
        return self.get_node()
    @node.setter
    def node(self, node):
        self.set_node(node)
    def get_node(self):
        return self.material.node_tree.nodes
    def set_node(self, node):
        for key, value in node.items():
            self.material.node_tree.nodes['Principled BSDF'].inputs[key].default_value = value
    @property
    def subdivisions(self):
        return self.get_subdivisions()
    @subdivisions.setter
    def subdivisions(self, subdivisions):
        self.set_subdivisions(subdivisions)
    def get_subdivisions(self):
        nverts = len(self.mesh.data.vertices)
        return nverts
    def set_subdivisions(self, subdivisions):
        if not isinstance(subdivisions, int):
            raise Exception('subdivisions should be int!')
        self.clean_bbonds_objects('mesh_bond_%s_%s'%(self.label, self.species))
        mesh = self.set_mesh(subdivisions = subdivisions, shape='ICO_SPHERE')
        mesh.parent = self.obj
    
    def clean_bbonds_objects(self, obj):
        obj = bpy.data.objects[obj]
        bpy.data.objects.remove(obj, do_unlink = True)
    
    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        """

        frames: list
            list of positions
        
        >>> from bbonds import Bbond
        >>> import numpy as np
        >>> positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
        >>> h = Bbond('h2o', 'H', positions)
        >>> frames = []
        >>> for i in range(10):
                frames.append(positions + [0, 0, i])
        >>> h.set_frames(frames)
        
        use shape_keys (faster)
        """
        from batoms.butils import add_keyframe_to_shape_key
        from batoms.source_data import bond_source
        from batoms.bdraw import cylinder_mesh_from_vec
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0 : return
        obj = self.obj
        base_name = 'Basis_%s'%self.species
        if obj.data.shape_keys is None:
            obj.shape_key_add(name = base_name)
        elif base_name not in obj.data.shape_keys.key_blocks:
            obj.shape_key_add(name = base_name)
        if only_basis:
            return
        nvert = len(obj.data.vertices)
        mesh_data = bond_source[self.segments]
        for i in range(1, nframe):
            sk = obj.shape_key_add(name = str(i))
            # Use the local position here
            positions = frames[i]
            vertices, faces = cylinder_mesh_from_vec(positions[:, 0:3], 
                                positions[:, 3:6], positions[:, 6:7], 
                                self.width, mesh_data)
            vertices = vertices.reshape((nvert*3, 1))
            sk.data.foreach_set('co', vertices)
            # Add Keyframes, the last one is different
            if i != nframe - 1:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1, 0], [frame_start + i - 1, 
                    frame_start + i, frame_start + i + 1])
            else:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1], [frame_start + i - 1, frame_start + i])

    def __len__(self):
        n = int(len(self.obj.data.vertices)/self.segments/2)
        return n
    
    def __getitem__(self, index):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system
        
        """
        return self.positions[index]
        
    def __setitem__(self, index, value):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[index] = value
        self.set_positions(positions)

    def repeat(self, m, cell):
        """
        In-place repeat of atoms.

        >>> from batoms.bbond import Bbond
        >>> c = Bbond('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
        >>> c.repeat([3, 3, 3], np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]))
        """
        if isinstance(m, int):
            m = (m, m, m)
        for x, vec in zip(m, cell):
            if x != 1 and not vec.any():
                raise ValueError('Cannot repeat along undefined lattice '
                                 'vector')
        M = np.product(m)
        n = len(self)
        positions = np.tile(self.positions, (M,) + (1,) * (len(self.positions.shape) - 1))
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
        self.add_vertices(positions[n:])
    def copy(self, label, species):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(label = 'h_new', species = 'H')

        """
        object_mode()
        bbond = Bbond(label, species, self.local_positions, 
                    location = self.obj.location, 
                    scale = self.scale, material=self.material)
        return bbond
    def extend(self, other):
        """
        Extend bbond object by appending bbond from *other*.
        
        >>> from batoms.bbonds import Bbond
        >>> h1 = Bbond('h2o', 'H_1', [[0, 0, 0], [2, 0, 0]])
        >>> h2 = Bbond('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
        >>> h = h1 + h2
        """
        # could also use self.add_vertices(other.positions)
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        other.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.join()
    def __iadd__(self, other):
        """
        >>> h1 += h2
        """
        self.extend(other)
        return self
    def __add__(self, other):
        """
        >>> h1 = h1 + h2
        """
        self += other
        return self
    def __iter__(self):
        bbond = self.obj
        for i in range(len(self)):
            yield bbond.matrix_world @ bbond.data.vertices[i].co
    def __repr__(self):
        s = "Bbond('%s', positions = %s" % (self.species, list(self.positions))
        return s
    def add_vertices(self, positions):
        """
        Todo: find a fast way.
        """
        object_mode()
        positions = positions - self.location
        bm = bmesh.new()
        bm.from_mesh(self.obj.data)
        bm.verts.ensure_lookup_table()
        verts = []
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.obj.data)