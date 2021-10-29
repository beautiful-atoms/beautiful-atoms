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

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE"]


class Batom():
    """Batom Class
    
    A Batom object is linked to this main collection in Blender. 

    Parameters:

    label: str
        Name of the Batoms.
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

    >>> from batoms.batom import Batom
    >>> c = Batom('C', [[0, 0, 0], [1.2, 0, 0]])
    >>> c.draw_atom()

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
                 ):
        #
        if species:
            self.label = label
            self.species = species
            if not element:
                self.element = species.split('_')[0]
            elif isinstance(element, str):
                self.element = element
            elif isinstance(element, dict):
                self.element = element
            self.name = 'atom_%s_%s'%(self.label, self.species)
            self.props = props
            self.species_data = get_default_species_data(self.element, scale = scale, 
                                props = self.props, radii_style = radii_style, 
                                color_style = color_style)
            if color:
                self.species_data['color'] = color
            self.set_material(node_inputs, material_style, material)
            self.set_instancer(segments = segments, 
                            subdivisions = subdivisions, shape = shapes[shape])
            self.set_object(positions, location)
            self.batom.batom.radius = self.species_data['radius']
        else:
            self.from_batom(label)
    def set_material(self, node_inputs = None, material_style = 'default', material = None):
        """
        """
        from batoms.material import create_material
        name = 'material_atom_{0}_{1}'.format(self.label, self.species)
        if material:
            material = material.copy()
            material.name = name
        elif name not in bpy.data.materials:
            material = create_material(name,
                        self.species_data['color'],
                        node_inputs = node_inputs,
                        material_style = material_style,
                        backface_culling = True)
    def object_mode(self):
        for object in bpy.data.objects:
            if object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode = 'OBJECT')
    def set_instancer(self, segments = [32, 16], subdivisions = 2, 
                        shape = 'UV_SPHERE', shade_smooth = True):
        object_mode()
        name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
        if name not in bpy.data.objects:
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
            obj.name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
            obj.data.materials.append(self.material)
            obj.data.name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
            if shade_smooth:
                bpy.ops.object.shade_smooth()
            obj.hide_set(True)
    def set_object(self, positions, location):
        """
        build child object and add it to main objects.
        """
        if self.name not in bpy.data.objects:
            mesh = bpy.data.meshes.new(self.name)
            obj_atom = bpy.data.objects.new(self.name, mesh)
            obj_atom.data.from_pydata(positions, [], [])
            obj_atom.location = location
            obj_atom.batom.flag = True
            obj_atom.batom.species = self.species
            obj_atom.batom.element = self.element
            obj_atom.batom.label = self.label
            bpy.data.collections['Collection'].objects.link(obj_atom)
        elif bpy.data.objects[self.name].batom.flag:
            obj_atom = bpy.data.objects[self.name]
        else:
            raise Exception("Failed, the name %s already in use and is not Batom object!"%self.name)
        self.instancer.parent = obj_atom
        obj_atom.instance_type = 'VERTS'
        bpy.context.view_layer.update()
    def from_batom(self, label):
        if label not in bpy.data.objects:
            raise Exception("%s is not a object!"%label)
        elif not bpy.data.objects[label].batom.flag:
            raise Exception("%s is not Batom object!"%label)
        ba = bpy.data.objects[label]
        self.species = ba.batom.species
        self.label = ba.batom.label
        self.element = ba.batom.element
        self.species_data = {
            'radius':ba.batom.radius,
            'scale':ba.scale,
        }
    @property
    def scene(self):
        return self.get_scene()
    def get_scene(self):
        return bpy.data.scenes['Scene']
    @property
    def batom(self):
        return self.get_batom()
    def get_batom(self):
        return bpy.data.objects['atom_%s_%s'%(self.label, self.species)]
    @property
    def instancer(self):
        return self.get_instancer()
    def get_instancer(self):
        return bpy.data.objects['instancer_atom_%s_%s'%(self.label, self.species)]
    @property
    def material(self):
        return self.get_material()
    def get_material(self):
        return bpy.data.materials['material_atom_%s_%s'%(self.label, self.species)]
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
        return np.array(self.batom.batom.radius)
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
    def location(self):
        return self.get_location()
    def get_location(self):
        return np.array(self.batom.location)
    @property
    def local_positions(self):
        return self.get_local_positions()
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = np.empty(n*3, dtype=np.float64)
        self.batom.data.vertices.foreach_get('co', local_positions)  
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
        local_positions = self.local_positions
        n = len(local_positions)
        # positions (natom, 3) to (natom, 4)
        local_positions = np.append(local_positions, np.ones((n, 1)), axis = 1)
        mat= np.array(self.batom.matrix_world)
        positions = mat.dot(local_positions.T).T
        # (natom, 4) back to (natom, 3)
        positions = positions[:, :3]
        return positions
    def set_positions(self, positions):
        """
        Set global positions to local vertices
        """
        natom = len(self)
        if len(positions) != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natom))
        mat= np.array(self.batom.matrix_world)
        # positions (natom, 3) to (natom, 4)
        n = len(positions)
        positions = np.append(positions, np.ones((n, 1)), axis = 1)
        # inverse of transformation matrix
        mat = np.linalg.inv(mat)
        local_positions = mat.dot(positions.T).T
        # (natom, 4) back to (natom, 3)
        local_positions = local_positions[:, :3]
        # rashpe to (natoms*3, 1) and use forseach_set
        local_positions = local_positions.reshape((natom*3, 1))
        self.batom.data.vertices.foreach_set('co', local_positions)
        self.batom.data.update()
    def get_scaled_positions(self, cell):
        """
        Get array of scaled_positions.
        """
        from ase.cell import Cell
        cell = Cell.new(cell)
        scaled_positions = cell.scaled_positions(self.local_positions)
        return scaled_positions
    @property
    def frames(self):
        return self.get_frames()
    @frames.setter
    def frames(self, frames):
        self.set_frames(frames)
    def get_frames(self):
        frames = []
        for f in range(bpy.context.scene.frame_start, bpy.context.scene.frame_end):
            bpy.context.scene.frame_set(f)
            frames.append(self.positions)
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
        self.set_instancer(segments = segments)
        self.instancer.parent = self.batom
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
        self.set_instancer(subdivisions = subdivisions, shape='ICO_SPHERE')
        self.instancer.parent = self.batom
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
        if shape not in [0, 1, 2]:
            raise Exception('Shape %s is not supported!'%shape)
        self.clean_batoms_objects('instancer_atom_%s_%s'%(self.label, self.species))
        self.set_instancer(shape = shapes[shape])
        self.instancer.parent = self.batom
    def clean_batoms_objects(self, obj):
        obj = bpy.data.objects[obj]
        bpy.data.objects.remove(obj, do_unlink = True)
    def delete_verts(self, index = []):
        """
        delete verts
        """
        object_mode()
        obj = self.batom
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        verts_select = [bm.verts[i] for i in index] 
        bmesh.ops.delete(bm, geom=verts_select, context='VERTS')
        if len(bm.verts) == 0:
            bpy.data.objects.remove(obj)
            bpy.data.objects.remove(self.instancer)
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
    
    def set_frames(self, frames = []):
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
        >>> h.load_frames(frames)
        
        Becareful!
        Only one frame, we will not register keyframe.
        If we registered keyframes here. Any change or the batoms has to be registered.
        Register new change before rendering, or delete existing keyframes.
        """
        nframe = len(frames)
        if nframe <= 1: return
        batom = self.batom
        nverts = len(batom.data.vertices)
        for i in range(0, nframe):
            positions = frames[i]
            for j in range(nverts):
                p = np.array(positions[j])
                batom.data.vertices[j].co = np.array(positions[j]) - np.array(batom.location)
                batom.data.vertices[j].keyframe_insert('co', frame=i + 1)
    
    def __len__(self):
        return len(self.batom.data.vertices)
    
    def __getitem__(self, index):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system
        
        """
        return self.positions[index]
        # batom = self.batom
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

        # batom  =self.batom
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
        batom = Batom(label, species, self.local_positions, 
                    location = self.batom.location, 
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
        self.batom.select_set(True)
        other.batom.select_set(True)
        bpy.context.view_layer.objects.active = self.batom
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
        batom = self.batom
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
        bm.from_mesh(self.batom.data)
        bm.verts.ensure_lookup_table()
        verts = []
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.batom.data)
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h.translate([0, 0, 5])
        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.batom.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.batom.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), 
                        orient_type = orient_type)
    def get_cell(self):
        if not self.label in bpy.data.collections:
            return None
        bcell = bpy.data.collections['%s_cell'%self.label]
        cell = np.array([bcell.matrix_world @ bcell.data.vertices[i].co for i in range(3)])
        return cell
    @property
    def hide(self):
        return self.get_hide()
    @hide.setter
    def hide(self, state):
        self.set_hide(state)
    def get_hide(self):
        return 'Unknown'
    def set_hide(self, state):
        self.batom.hide_render = state
        self.batom.hide_set(state)
    @property
    def select(self):
        return self.get_select()
    @select.setter
    def select(self, state):
        self.set_select(state)
    def get_select(self):
        return 'Unknown'
    def set_select(self, state):
        self.batom.select_set(state)