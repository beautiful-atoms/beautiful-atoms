"""Definition of the Batom class.

This module defines the Batom object in the batoms package.

"""

from time import time
import bpy
import bmesh
from batoms.butils import object_mode
from batoms.material import material_styles_dict
from batoms.tools import get_default_species_data
import numpy as np
from batoms.base import BaseObject

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE"]


class Batom(BaseObject):
    """Batom Class
    
    A Batom object is linked to this main collection in Blender. 

    Parameters:

    label: str
        Name of the Batoms.
    species: str or dict
        species of the atoms.
        str: 'O', 'O_1', 'Fe_up', 'Fe_3+',
        dict: {'Al': 0.9, 'Si': 0.1}
    positions: array
        positions
    locations: array
        The object’s origin location in global coordinates.
    element: str or dict
        element of the atoms, dict for fractional Occupancy
        str: 'O'
        dict: {'Al': 'Al', 'Si': 'Si'}
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

    >>> from batoms.batom import Batom
    >>> c = Batom('C', [[0, 0, 0], [1.2, 0, 0]])
    >>> al = Batom({'Al': 0.7, 'Si': 0.3}, [[0, 0, 0], [1.2, 0, 0]])

    """
    

    def __init__(self, 
                label = None,
                species = None,
                positions = None,
                location = np.array([0, 0, 0]),
                element = None,
                scale = 1.0, 
                segments = [32, 16],
                shape = 0,
                subdivisions = 2,
                props = {},
                color = None,
                radii_style = 'covalent',
                color_style = 'JMOL',
                material_style = 'default',
                material = None,
                node_inputs = None,
                instancer_data = None,
                 ):
        #
        if species:
            self.label = label
            self.species = species
            self.name = species
            obj_name = '%s_atom_%s'%(self.label, self.species)
        else:
            obj_name = label
        BaseObject.__init__(self, obj_name = obj_name)
        if species:
            positions = np.array(positions)
            if len(positions.shape) == 2:
                self._frames = np.array([positions])
            elif len(positions.shape) == 3:
                self._frames = positions
                positions = self._frames[0]
            else:
                raise Exception('Shape of positions is wrong!')
            if not element:
                self.element = species.split('_')[0]
            elif isinstance(element, str):
                self.element = element
            elif isinstance(element, dict):
                self.element = element
            self.props = props
            self.species_data = get_default_species_data(self.element, scale = scale, 
                                props = self.props, radii_style = radii_style, 
                                color_style = color_style)
            if color:
                self.species_data['color'] = color
            self.set_material(node_inputs, material_style, material)
            instancer = self.set_instancer(segments = segments, 
                            subdivisions = subdivisions, shape = shapes[shape],
                            instancer_data = instancer_data)
            self.set_object(positions, location, instancer)
            self.obj.batom.radius = self.species_data['radius']
            self.set_frames(self._frames, only_basis = True)
        else:
            self.from_batom(label)
    def set_material(self, node_inputs = None, material_style = 'default', material = None):
        """
        """
        from batoms.material import create_material
        name = '{0}_material_atom_{1}'.format(self.label, self.species)
        if material:
            material = material.copy()
            material.name = name
        elif name not in bpy.data.materials:
            material = create_material(name,
                        self.species_data['color'],
                        node_inputs = node_inputs,
                        material_style = material_style,
                        backface_culling = True)
    def set_instancer(self, segments = [32, 16], subdivisions = 2, 
                        shape = 'UV_SPHERE', shade_smooth = True, instancer_data = None):
        object_mode()
        name = '{0}_instancer_atom_{1}'.format(self.label, self.species)
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        if instancer_data is not None:
            obj = bpy.data.objects.new(name, instancer_data)
        else:
            if shape.upper() == 'UV_SPHERE':
                bpy.ops.mesh.primitive_uv_sphere_add(segments = segments[0], 
                                    ring_count = segments[1], 
                                    radius = self.species_data['radius'])
            if shape.upper() == 'ICO_SPHERE':
                shade_smooth = False
                bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = subdivisions, 
                            radius = self.species_data['radius'])
            if shape.upper() == 'CUBE':
                bpy.ops.mesh.primitive_cube_add(size = self.species_data['radius'])
                shade_smooth = False
            obj = bpy.context.view_layer.objects.active
            if isinstance(self.species_data['scale'], float):
                self.species_data['scale'] = [self.species_data['scale']]*3
            obj.scale = self.species_data['scale']
            obj.name = '{0}_instancer_atom_{1}'.format(self.label, self.species)
            obj.data.name = '{0}_instancer_atom_{1}'.format(self.label, self.species)
        obj.data.materials.append(self.material)
        if shade_smooth:
            bpy.ops.object.shade_smooth()
        obj.hide_set(True)
        return obj
    def set_object(self, positions, location, instancer):
        """
        build child object and add it to main objects.
        """
        if self.obj_name not in bpy.data.objects:
            mesh = bpy.data.meshes.new(self.obj_name)
            obj = bpy.data.objects.new(self.obj_name, mesh)
            obj.data.from_pydata(positions, [], [])
            obj.location = location
            obj.batom.flag = True
            obj.batom.species = self.species
            obj.batom.element = self.element
            obj.batom.label = self.label
            bpy.data.collections['Collection'].objects.link(obj)
        elif bpy.data.objects[self.obj_name].batom.flag:
            obj = bpy.data.objects[self.obj_name]
        else:
            raise Exception("Failed, the name %s already in use and is not Batom object!"%self.obj_name)
        instancer.parent = obj
        obj.instance_type = 'VERTS'
        bpy.context.view_layer.update()
    def from_batom(self, label):
        if label not in bpy.data.objects:
            raise Exception("%s is not a object!"%label)
        elif not bpy.data.objects[label].batom.flag:
            raise Exception("%s is not Batom object!"%label)
        obj = bpy.data.objects[label]
        self.species = obj.batom.species
        self.label = obj.batom.label
        self.element = obj.batom.element
        self.species_data = {
            'radius':obj.batom.radius,
            'scale':obj.scale,
        } 
    @property
    def instancer(self):
        return self.get_instancer()
    def get_instancer(self):
        return self.obj.children[0]
    @property
    def material(self):
        return self.get_material()
    def get_material(self):
        return bpy.data.materials['%s_material_atom_%s'%(self.label, self.species)]
    @property
    def scale(self):
        return self.get_scale()
    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)
    def get_scale(self):
        return np.array(self.instancer.scale)
    def set_scale(self, scale):
        if isinstance(scale, float) or isinstance(scale, int):
            scale = [scale]*3
        self.instancer.scale = scale
    @property
    def radius(self):
        return self.get_radius()
    def get_radius(self):
        return np.array(self.obj.batom.radius)
    @property
    def size(self):
        return self.get_size()
    @size.setter
    def size(self, size):
        self.set_size(size)
    def get_size(self):
        return np.array(self.instancer.scale*self.radius)
    def set_size(self, size):
        scale = size/self.radius
        self.scale = [scale]*3
    @property
    def local_positions(self):
        return self.get_local_positions()
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = np.empty(n*3, dtype=np.float64)
        self.obj.data.vertices.foreach_get('co', local_positions)  
        local_positions = local_positions.reshape((n, 3))
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
    def segments(self):
        return self.get_segments()
    @segments.setter
    def segments(self, segments):
        self.set_segments(segments)
    def get_segments(self):
        nverts = len(self.instancer.data.vertices)
        return nverts
    def set_segments(self, segments):
        if not isinstance(segments, int):
            raise Exception('Segments should be int!')
        self.clean_batoms_objects('instancer_atom_%s_%s'%(self.label, self.species))
        instancer = self.set_instancer(segments = segments)
        instancer.parent = self.obj
    @property
    def subdivisions(self):
        return self.get_subdivisions()
    @subdivisions.setter
    def subdivisions(self, subdivisions):
        self.set_subdivisions(subdivisions)
    def get_subdivisions(self):
        nverts = len(self.instancer.data.vertices)
        return nverts
    def set_subdivisions(self, subdivisions):
        if not isinstance(subdivisions, int):
            raise Exception('subdivisions should be int!')
        self.clean_batoms_objects('instancer_atom_%s_%s'%(self.label, self.species))
        instancer = self.set_instancer(subdivisions = subdivisions, shape='ICO_SPHERE')
        instancer.parent = self.obj
    @property
    def shape(self):
        return self.get_shape()
    @shape.setter
    def shape(self, shape):
        self.set_shape(shape)
    def get_shape(self):
        """
        todo"""
        # nverts = len(self.instancer.data.vertices)
        return 'to do'
    def set_shape(self, shape):
        """
        "UV_SPHERE", "ICO_SPHERE", "CUBE"
        """
        scale = self.scale
        if shape not in [0, 1, 2]:
            raise Exception('Shape %s is not supported!'%shape)
        self.clean_batoms_objects('instancer_atom_%s_%s'%(self.label, self.species))
        instancer = self.set_instancer(shape = shapes[shape])
        instancer.parent = self.obj
        self.scale = scale
    def clean_batoms_objects(self, obj):
        obj = bpy.data.objects[obj]
        bpy.data.objects.remove(obj, do_unlink = True)
    def delete_verts(self, index = []):
        """
        delete verts
        """
        object_mode()
        obj = self.obj
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        verts_select = [bm.verts[i] for i in index] 
        bmesh.ops.delete(bm, geom=verts_select, context='VERTS')
        if len(bm.verts) == 0:
            bpy.data.objects.remove(self.instancer)
            bpy.data.objects.remove(obj)
        else:
            bm.to_mesh(obj.data)
    def delete(self, index = []):
        """
        delete atom.

        index: list
            index of atoms to be delete
        
        For example, delete the second atom in h species. 
        Please note that index start from 0.

        >>> h.delete([1])

        """
        if isinstance(index[0], (bool, np.bool_)):
            index = np.where(index)[0]
        if isinstance(index, int):
            index = [index]
        self.delete_verts(index)
    def __delitem__(self, index):
        """
        """
        self.delete(index)
    def draw_constraints(self):
        """
        """
        #
        constr = self.atoms.constraints
        self.constrainatom = []
        for c in constr:
            if isinstance(c, FixAtoms):
                for n, i in enumerate(c.index):
                    self.constrainatom += [i]
    
    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        """

        frames: list
            list of positions
        
        >>> from batoms import Batom
        >>> import numpy as np
        >>> positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
        >>> h = Batom('h2o', 'H', positions)
        >>> frames = []
        >>> for i in range(10):
                frames.append(positions + [0, 0, i])
        >>> h.set_frames(frames)
        
        use shape_keys (faster)
        """
        from batoms.butils import add_keyframe_to_shape_key
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
        for i in range(1, nframe):
            sk = obj.data.shape_keys.key_blocks.get(str(i))
            if sk is None:
                sk = obj.shape_key_add(name = str(i))
            # Use the local position here
            positions = frames[i].reshape((nvert*3, 1))
            sk.data.foreach_set('co', positions)
            # Add Keyframes, the last one is different
            if i != nframe - 1:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1, 0], [frame_start + i - 1, 
                    frame_start + i, frame_start + i + 1])
            else:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1], [frame_start + i - 1, frame_start + i])

    def __len__(self):
        return len(self.obj.data.vertices)
    
    def __getitem__(self, index):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system
        
        """
        return self.positions[index]
        # batom = self.obj
        # if isinstance(index, int):
        #     natom = len(self)
        #     if index < -natom or index >= natom:
        #         raise IndexError('Index out of range.')
        #     return batom.matrix_world @ batom.data.vertices[index].co
        # if isinstance(index, list):
        #     positions = np.array([self[i] for i in index])
        #     return positions
        # if isinstance(index, slice):
        #     start, stop, step = index.indices(len(self))
        #     index = list(range(start, stop, step))
        #     return self[index]
        # if isinstance(index, tuple):
        #     i, j = index
        #     positions = self[i]
        #     return positions[:, j]

    def __setitem__(self, index, value):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[index] = value
        self.set_positions(positions)

        # batom  =self.obj
        # if isinstance(index, int):
        #     natom = len(self)
        #     if index < -natom or index >= natom:
        #         raise IndexError('Index out of range.')
        #     batom.data.vertices[index].co = np.array(value) - np.array(batom.location)
        # if isinstance(index, list):
        #     for i in index:
        #         self[i] = value[i]
        # if isinstance(index, tuple):
        #     i, j = index
        #     batom.data.vertices[i].co[j] = np.array(value) - np.array(batom.location[j])

    def repeat(self, m, cell):
        """
        In-place repeat of atoms.

        >>> from batoms.batom import Batom
        >>> c = Batom('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
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
        frames = self.frames
        positions = self.positions
        positions = np.tile(positions, (M,) + (1,) * (len(positions.shape) - 1))
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
        self.add_vertices(positions[n:])
        # repeat frames
        frames_new = []
        if self.nframe > 1:
            for i in range(0, self.nframe):
                positions = np.tile(frames[i], (M,) + (1,) * (len(frames[i].shape) - 1))
                i0 = 0
                for m0 in range(m[0]):
                    for m1 in range(m[1]):
                        for m2 in range(m[2]):
                            i1 = i0 + n
                            positions[i0:i1] += np.dot((m0, m1, m2), cell)
                            i0 = i1
                frames_new.append(positions)
        self.set_frames(frames_new)

    def copy(self, label, species):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(label = 'h_new', species = 'H')

        """
        object_mode()
        batom = Batom(label, species, self.local_positions, 
                    location = self.obj.location, 
                    scale = self.scale, material=self.material)
        return batom
    def extend(self, other):
        """
        Extend batom object by appending batom from *other*.
        
        >>> from batoms.batoms import Batom
        >>> h1 = Batom('h2o', 'H_1', [[0, 0, 0], [2, 0, 0]])
        >>> h2 = Batom('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
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
        batom = self.obj
        for i in range(len(self)):
            yield batom.matrix_world @ batom.data.vertices[i].co
    def __repr__(self):
        s = "Batom('%s', positions = %s" % (self.species, list(self.positions))
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
    def get_cell(self):
        if not self.label in bpy.data.collections:
            return None
        bcell = bpy.data.collections['%s_cell'%self.label]
        cell = np.array([bcell.matrix_world @ bcell.data.vertices[i].co for i in range(3)])
        return cell