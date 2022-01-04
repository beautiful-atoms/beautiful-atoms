"""Definition of the Batom class.

This module defines the Batom object in the batoms package.

"""

from time import time

import bpy
import bmesh
from batoms.butils import object_mode
from batoms.tools import get_default_species_data
import numpy as np
from batoms.base import BaseObject

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]


class Batom(BaseObject):
    """Batom Class
    
    A Batom object is linked to this Object in Blender. 

    Parameters:

    label: str
        Name of the Batoms.
    species: str
        species of the atoms.
        str: 'O', 'O_1', 'Fe_up', 'Fe_3+',
    positions: array
        positions
    attributes: array
        eg. index, chianid
    locations: array
        The object's origin location in global coordinates.
    element: str or dict
        element of the atoms, dict for fractional Occupancy
        str: 'O'
        dict: {'Al': 0.887, 'Si': 0.333}
    segments: list of 2 Int
        Number of segments used to draw the UV_Sphere
        Default: [32, 16]
    subdivisions: Int
        Number of subdivision used to draw the ICO_Sphere
        Default: 2
    color_style: str
        "JMOL", "ASE", "VESTA"
    radius_style: str
        "covelent", "vdw", "ionic"
    shape: Int
        0, 1, or 2. ["UV_SPHERE", "ICO_SPHERE", "CUBE"]

    Examples:

    >>> from batoms.batom import Batom
    >>> c = Batom('C', [[0, 0, 0], [1.2, 0, 0]])
    >>> si = Batom('zsm', 'si', element = {'Si': 0.833, 'Al': 0.167}, 
                positions = [[0, 0, 0], [1.52, 0, 0]])

    """
    def __init__(self, 
                name,
                species = None,
                positions = None,
                attributes = {},
                elements = None,
                location = np.array([0, 0, 0]),
                label = 'batoms',
                scale = 1.0, 
                segments = [32, 16],
                shape = 0,
                subdivisions = 2,
                props = {},
                color = None,
                radius_style = 'covalent',
                color_style = 'JMOL',
                material_style = 'default',
                materials = None,
                node_inputs = None,
                instancer_data = None,
                metaball = False,
                 ):
        #
        self.name = name
        obj_name = name
        self.instancer_name = 'instancer_%s'%name
        if species is None:
            species = name
        bobj_name = 'batom'
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        if positions is not None:
            positions = np.array(positions)
            if len(positions.shape) == 2:
                self._frames = np.array([positions])
            elif len(positions.shape) == 3:
                self._frames = positions
                positions = self._frames[0]
            else:
                raise Exception('Shape of positions is wrong!')
            if not elements:
                elements = species.split('_')[0]
            elements = self.check_elements(elements)
            self.props = props
            self.species_data = get_default_species_data(elements,
                                radius_style = radius_style, 
                                color_style = color_style,
                                props = self.props)
            if color:
                self.species_data['color'] = color
            self.build_object(species, elements = elements, 
                        positions = positions, 
                        attributes = attributes, 
                        location = location,
                        label = label,
                        metaball = metaball)
            self.build_material(species, label = label, 
                            node_inputs = node_inputs, 
                            material_style = material_style, 
                            materials = materials)
            self.build_instancer(radius = self.species_data['radius'], 
                            scale = scale,
                            segments = segments, 
                            subdivisions = subdivisions, shape = shapes[shape],
                            instancer_data = instancer_data)
            self.set_frames(self._frames, only_basis = True)
        else:
            self.from_batom(name)
            self.instancer_name = 'instancer_{0}'.format(self.name)
    
    def build_material(self, species, label = 'batoms', node_inputs = None, 
                material_style = 'default', materials = None):
        """
        """
        from batoms.material import create_material
        if materials is not None:
            for name, mat in materials.items():
                mat1 = mat.copy()
                label = name.split('_')[0]
                mat1.name = name.replace(label, label)
        for ele, color in self.species_data['color'].items():
            name = '{0}_material_atom_{1}_{2}'.format(label, species, ele)
            if name not in bpy.data.materials:
                create_material(name,
                            color = color,
                            node_inputs = node_inputs,
                            material_style = material_style,
                            backface_culling = True)
    
    def build_object(self, species, elements, positions, attributes, location, label = 'batoms', metaball = False):
        """
        build child object and add it to main objects.
        """
        if self.obj_name not in bpy.data.objects:
            mesh = bpy.data.meshes.new(self.obj_name)
            obj = bpy.data.objects.new(self.obj_name, mesh)
            obj.data.from_pydata(positions, [], [])
            obj.location = location
            obj.batoms.batom.flag = True
            obj.batoms.batom.species = species
            obj.batoms.batom.label = label
            for ele, occupancy in elements.items():
                eledata = obj.batoms.batom.elements.add()
                eledata.name = ele
                eledata.occupancy = occupancy
            bpy.data.collections['Collection'].objects.link(obj)
            # Assigning attributes to vertexes
            self.set_attributes(attributes)
        elif bpy.data.objects[self.obj_name].batoms.batom.flag:
            obj = bpy.data.objects[self.obj_name]
        else:
            raise Exception("Failed, the name %s already in use and is not Batom object!"%self.obj_name)
    
    def build_instancer(self, radius = 1.0, scale = [1, 1, 1], segments = [32, 16], subdivisions = 2, 
                        shape = 'UV_SPHERE', shade_smooth = True, instancer_data = None):
        object_mode()
        name = self.instancer_name
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        if instancer_data is not None:
            obj = bpy.data.objects.new(name, instancer_data)
        else:
            if shape.upper() == 'UV_SPHERE':
                bpy.ops.mesh.primitive_uv_sphere_add(segments = segments[0], 
                                    ring_count = segments[1], 
                                    radius = radius)
            if shape.upper() == 'ICO_SPHERE':
                shade_smooth = False
                bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = subdivisions, 
                            radius = radius)
            if shape.upper() == 'CUBE':
                bpy.ops.mesh.primitive_cube_add(size = radius)
                shade_smooth = False
            if shape.upper() == 'METABALL':
                bpy.ops.object.metaball_add(type = 'BALL', location = [0, 0, 0])
            obj = bpy.context.view_layer.objects.active
            if isinstance(scale, float):
                scale = [scale]*3
            obj.scale = scale
            obj.batoms.batom.radius = radius
            obj.name = self.instancer_name
            obj.data.name = self.instancer_name
        #
        # obj.data.materials.append(self.material)
        if shade_smooth:
            bpy.ops.object.shade_smooth()
        if shape.upper() != 'METABALL':
            obj.hide_set(True)
        obj.parent = self.obj
        self.obj.instance_type = 'VERTS'
        # self.obj.show_instancer_for_viewport = False
        # self.obj.show_instancer_for_render = False
        #
        self.assign_materials()
        bpy.context.view_layer.update()
        return obj
    
    def draw_SAS_mb(self, render_resolution = 0.2,
                resolution = 0.4,
                threshold = 1e-4,
                stiffness = 1,
                indices = None,
                update_method = 'FAST'
                ):
        """
        Algorithm: Metaball
        Computes density from given metaball at given position.
        Metaball equation is: 
        dens = `(1 - r^2 / R^2)^3 * s`
        r = distance from center
        R = metaball radius
        s - metaball stiffness
        field = threshold - dens;
        """
        from batoms.material import create_material
        stiffness = min(stiffness, 10)
        frames = self.batoms.frames
        n = len(frames[0])
        if indices is None:
            indices = range(n)
        radii = self.batoms.radii_vdw
        tstart = time()
        #
        self.build_metaball(render_resolution = render_resolution,
                resolution = resolution,
                threshold = threshold,
                update_method = update_method)
        color = default_colors[0]
        mat = create_material(self.sas_name,
                        # material_style = 'plastic',
                        color = color)
        obj = self.sas_obj
        mb = obj.data
        obj.data.materials.append(mat)
        atoms = frames[0]
        positions = atoms.positions
        scale = (1 - (threshold/stiffness)**(1.0/3))**(0.5)
        # print('scale: %s', scale)
        mb.elements.clear()
        for i in indices:
            mbele = mb.elements.new(type = 'BALL')
            mbele.co = positions[i]
            mbele.radius = (radii[i] + self.probe)/scale
            mbele.stiffness = stiffness
        # add frames
        nframe = len(frames)
        self.nframe = nframe
        for i in range(1, nframe):
            positions = frames[i].positions
            for j in indices:
                mb.elements[j].co = positions[j]
                mb.elements[j].keyframe_insert(data_path='co', frame = i)
        print('Build_SAS: %s'%(time() - tstart))
        return obj

    def assign_materials(self):
        # sort element by occu
        mesh = self.instancer.data
        mesh.materials.clear()
        sorted_ele = sorted(self.elements.items(), key=lambda x: -x[1])
        materials = self.materials
        for data in sorted_ele:
            mesh.materials.append(materials[data[0]])
        # find the face index for ele
        nele = len(sorted_ele)
        # for occupancy
        if nele > 1:
            # calc the angles
            npoly = len(mesh.polygons)
            normals = np.zeros(npoly*3)
            material_indexs = np.zeros(npoly, dtype='int')
            mesh.polygons.foreach_get('normal', normals)
            mesh.polygons.foreach_get('material_index', material_indexs)
            normals = normals.reshape(-1, 3)
            xy = normals - np.dot(normals, [0, 0, 1])[:, None]*[0, 0, 1]
            angles = (np.arctan2(xy[:, 1], xy[:, 0]) + np.pi)/np.pi/2
            # 
            tos = 0
            for i in range(1, nele):
                toe = tos + sorted_ele[i][1]
                index = np.where((angles > tos) & (angles < toe))[0]
                material_indexs[index] = i
                tos = toe
            mesh.polygons.foreach_set('material_index', material_indexs)

    def from_batom(self, name):
        if name not in bpy.data.objects:
            raise Exception("%s is not a object!"%name)
        elif not bpy.data.objects[name].batoms.batom.flag:
            raise Exception("%s is not Batom object!"%name)
        obj = bpy.data.objects[name]
        # self.label = obj.batoms.batom.label
        self.species_data = {
            'radius':obj.batoms.batom.radius,
            'scale':obj.scale,
        }
    
    @property
    def instancer(self):
        return self.get_instancer()
    
    def get_instancer(self):
        instancer = bpy.data.objects.get(self.instancer_name)
        return instancer
    
    @property
    def materials(self):
        return self.get_materials()
    
    def get_materials(self):
        materials = {}
        for ele in self.elements:
            name = '%s_material_atom_%s_%s'%(self.label, self.species, ele)
            mat = bpy.data.materials.get(name)
            materials[ele] = mat
        return materials
    
    @property
    def label(self):
        return self.get_label()
    
    @label.setter
    def label(self, label):
        self.set_label(label)
    
    def get_label(self):
        return self.obj.batoms.batom.label
    
    def set_label(self, label):
        self.obj.batoms.batom.label = label

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
    def species(self):
        return self.get_species()
    
    def get_species(self):
        return self.obj.batoms.batom.species
    
    @species.setter
    def species(self, species):
        self.obj.batoms.batom.species = species
    
    @staticmethod
    def check_elements(elements):
        if isinstance(elements, str):
            elements = {elements: 1.0}
        elif isinstance(elements, dict):
            elements = elements
            occu = sum(elements.values())
            # not fully occupied. 
            if occu < 1 - 1e-6:
                elements['X'] = 1 - occu
            elif occu > 1 + 1e-6:
                raise ValueError("Total occumpancy should be smaller than 1!")
        return elements

    @property
    def elements(self):
        return self.get_elements()
    
    def get_elements(self):
        elements = {}
        collection = self.obj.batoms.batom.elements
        for eledata in collection:
            elements[eledata.name] = round(eledata.occupancy, 3)
        return elements
    
    @elements.setter
    def elements(self, elements):
        self.set_elements(elements)
    
    def set_elements(self, elements):
        elements = self.check_elements(elements)
        collection = self.obj.batoms.batom.elements
        collection.clear()
        for ele, occupancy in elements.items():
            eledata = collection.add()
            eledata.name = ele
            eledata.occupancy = occupancy
        # build materials
        self.species_data = get_default_species_data(elements)
        self.build_material(self.species)
        self.assign_materials()

    @property
    def main_element(self):
        sorted_ele = sorted(self.elements.items(), key=lambda x: -x[1])
        if sorted_ele[0][0] == 'X':
            ele = sorted_ele[1][0]
        else:
            ele = sorted_ele[0][0]
        return ele
    
    @property
    def radius(self):
        return self.get_radius()
    
    def get_radius(self):
        return np.array(self.instancer.batoms.batom.radius)
    
    @property
    def radius_style(self):
        return self.get_radius_style()
    
    @radius_style.setter
    def radius_style(self, radius_style):
        self.set_radius_style(radius_style)
    
    def get_radius_style(self):
        return self.obj.batoms.batom.radius_style
    
    def set_radius_style(self, radius_style):
        self.obj.batoms.batom.radius_style = str(radius_style)
        scale = self.scale
        self.clean_batom_objects(self.instancer_name)
        species_data = get_default_species_data(self.elements,
                                radius_style = radius_style)
        # print(species_data)
        instancer = self.build_instancer(radius = species_data['radius'], 
                                scale = scale)
        instancer.parent = self.obj

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
        object_mode()
        from batoms.tools import local2global
        natom = len(self)
        if len(positions) != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natom))
        positions = local2global(positions, 
                np.array(self.obj.matrix_world), reversed = True)
        # rashpe to (natoms*3, 1) and use forseach_set
        positions = positions.reshape((natom*3, 1))
        # I don't know why 'Basis' shape keys is not updated when editing mesh,
        # so we edit the 'Basis' shape keys directly.
        # self.obj.data.vertices.foreach_set('co', positions)
        self.obj.data.shape_keys.key_blocks[0].data.foreach_set('co', positions)
        self.obj.data.update()
        # bpy.context.view_layer.update()
        # I don't why this is need to update the mesh positions
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
    
    def get_scaled_positions(self, cell):
        """
        Get array of scaled_positions.
        """
        from ase.cell import Cell
        cell = Cell.new(cell)
        scaled_positions = cell.scaled_positions(self.local_positions)
        return scaled_positions
    
    @property
    def attributes(self):
        return self.get_attributes()
    
    @attributes.setter
    def attributes(self, attributes):
        self.set_attributes(attributes)

    def get_attributes(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        # attributes
        me = self.obj.data
        nvert = len(me.vertices)
        attributes = {}
        for key in me.vertex_layers_int.keys():
            attributes[key] = np.zeros(nvert, dtype = int)
            layer = me.vertex_layers_int.get(key)
            layer.data.foreach_get("value", attributes[key])
        for key in me.vertex_layers_float.keys():
            attributes[key] = np.zeros(nvert)
            layer = me.vertex_layers_float.get(key)
            layer.data.foreach_get("value", attributes[key])
        for key in me.vertex_layers_string.keys():
            attributes[key] = []
            layer = me.vertex_layers_string.get(key)
            for d in layer.data:
                attributes[key].append(d.value.decode())
            attributes[key] = np.array(attributes[key])
        return attributes
    
    def set_attributes(self, attributes):
        me = self.obj.data
        for key, data in attributes.items():
            # print(key)
            if key == 'positions':
                self.positions = data
                continue
            dtype = type(attributes[key][0])
            if dtype == np.int_:
                layer = me.vertex_layers_int.get(key)
                if layer is None:
                    layer = me.vertex_layers_int.new(name = key)
                layer.data.foreach_set("value", data)
            elif dtype == np.float_:
                layer = me.vertex_layers_float.get(key)
                if layer is None:
                    layer = me.vertex_layers_float.new(name = key)
                layer.data.foreach_set("value", data)
            elif dtype == np.str_:
                # foreach_set not work with str
                # can not set len of bytes when using me.vertex_layers_string,
                # so use bmesh
                data = data.astype(bytes)
                bm = bmesh.new()
                bm.from_mesh(me)
                layer = bm.verts.layers.string.get(key)
                if layer is None:
                    layer = bm.verts.layers.string.new(key)
                bm.verts.ensure_lookup_table()
                nvert = len(bm.verts)
                for i in range(nvert):
                    v = bm.verts[i]
                    v[layer] = data[i]
                bm.to_mesh(me)
                bm.clear()             
            else:
                raise Exception('Data type %s not supported!'%dtype)
    
    @property
    def index(self):
        return self.get_index()
    
    def get_index(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        # index
        me = self.obj.data
        nvert = len(me.vertices)
        index = np.zeros(nvert, dtype = int)
        layer = me.vertex_layers_int.get('index')
        layer.data.foreach_get("value", index)
        return index

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
        Viewpoint_color = self.materials[self.main_element].diffuse_color
        for node in self.materials[self.main_element].node_tree.nodes:
            if 'Base Color' in node.inputs:
                node_color = node.inputs['Base Color'].default_value[:]
            if 'Alpha' in node.inputs:
                Alpha = node.inputs['Alpha'].default_value
        color = [node_color[0], node_color[1], node_color[2], Alpha]
        return color
    
    def set_color(self, color):
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.materials[self.main_element].diffuse_color = color
        for node in self.materials[self.main_element].node_tree.nodes:
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
        return self.materials[self.main_element].node_tree.nodes
    
    def set_node(self, node):
        for key, value in node.items():
            self.materials[self.main_element].node_tree.nodes['Principled BSDF'].inputs[key].default_value = value
    
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
        scale = self.scale
        radius = self.radius
        self.clean_batom_objects(self.instancer_name)
        instancer = self.build_instancer(radius = radius, 
                                scale = scale, 
                                segments = segments)
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
        scale = self.scale
        radius = self.radius
        self.clean_batom_objects(self.instancer_name)
        instancer = self.build_instancer(radius = radius, 
                                scale = scale, 
                                subdivisions = subdivisions, shape='ICO_SPHERE')
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
        scale = self.scale
        radius = self.radius
        self.clean_batom_objects(self.instancer_name)
        instancer = self.build_instancer(radius = radius, 
                                scale = scale, 
                                shape = shapes[shape])
        instancer.parent = self.obj
    
    def clean_batom_objects(self, obj):
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
            bm.clear()
    
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

    
    def repeat(self, m, cell, natom):
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
        index = self.index
        positions = np.tile(positions, (M,) + (1,) * (len(positions.shape) - 1))
        index = np.tile(index, (M,) + (1,) * (len(index.shape) - 1))
        i0 = 0
        n1 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    index[i0:i1] += n1*natom
                    i0 = i1
                    n1 += 1
        self.add_vertices(positions[n:])
        self.set_attributes({'index': index})
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

    
    def copy(self, name):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(name = 'h_new', species = 'H')

        """
        object_mode()
        #
        obj = self.obj.copy()
        obj.data = self.obj.data.copy()
        obj.name = obj.data.name = name
        #
        instancer = self.instancer.copy()
        bpy.data.collections['Collection'].objects.link(instancer)
        instancer.hide_set(True)
        instancer.name = 'instancer_%s'%(name)
        instancer.parent = obj
        batom = self.__class__(name)
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
        s = "Batom('%s', species = '%s', elements = %s, \
                    positions = %s" % (self.name, self.species,  \
                    str(self.elements), list(self.positions))
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
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.obj.data)
        bm.clear()
    
    def get_cell(self):
        if not self.label in bpy.data.collections:
            return None
        bcell = bpy.data.collections['%s_cell'%self.label]
        cell = np.array([bcell.matrix_world @ bcell.data.vertices[i].co for i in range(3)])
        return cell
    
    def make_real(self):
        """
        """
        self.select = True
        bpy.ops.object.duplicates_make_real()
    