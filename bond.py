"""Definition of the Bbond class.

This module defines the Bbond object in the bbonds package.

"""

from time import time
import bpy
import bmesh
from batoms.butils import object_mode
from batoms.tools import get_default_species_data, string2Number, read_from_ase, read_from_pymatgen
import numpy as np
from batoms.base import BaseObject

default_attributes = [
        ['style', 'INT'], 
        ['species_index1', 'INT'], 
        ['species_index2', 'INT'], 
        ['show', 'BOOLEAN'], 
        ['select', 'INT'],
        ['order', 'INT'],
        ['model_style', 'INT'],
        ]

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

    >>> from batoms.bond import Bbond
    >>> c = Bbond('C', [[0, 0, 0], [1.2, 0, 0]])

    """
    

    
    def __init__(self, 
                label = None,
                species = None,
                datas = None,
                location = np.array([0, 0, 0]),
                battr_inputs = {},
                attributes = {},
                batoms = None,
                bondsetting = None,
                 ):
        #
        self.batoms = batoms
        self.bondsetting = bondsetting
        if species is not None:
            self.label = label
            self.species = species
            self.name = species
            obj_name = '%s_bond_center'%(self.label)
        else:
            obj_name = label
        bobj_name = 'bbond'
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        if datas is not None:
            centers = datas['centers']
            if len(centers.shape) == 2:
                self._frames = np.array([centers])
            elif len(centers.shape) == 3:
                self._frames = centers
                centers = self._frames[0]
            else:
                raise Exception('Shape of centers is wrong!')
            self.build_object(centers, datas['normals'], datas['lengths'],
                         width = datas['widths'])
            nbond = len(centers)
            show = np.ones(nbond, dtype = int)
            # elements = self.check_elements(elements)
            species_index1 = [string2Number(sp) for sp in datas['species1']]
            species_index2 = [string2Number(sp) for sp in datas['species2']]
            attributes.update({
                                # 'species': species, 
                               'species_index1': species_index1, 
                               'species_index2': species_index2, 
                               'show': show,
                               'model_style': datas['model_styles'],
                               'style': datas['styles'],
                               'order': datas['orders'],
                               })
            self.set_attributes(attributes)
            self.build_geometry_node()
            self.set_frames(self._frames, only_basis = True)
        else:
            self.from_bbond(label)
    
    def build_object(self, positions, normals, scales, width, battr_inputs = {}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        from scipy.spatial.transform import Rotation as R
        tstart = time()
        name = '%s_bond_center'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(positions, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = 'POINT')
        self.bondsetting.coll.objects.link(obj)
        #
        name = '%s_bond_normal'%self.label
        vec = np.cross([0.0, 0.0, 1], normals) + np.array([1e-6, 0, 0])
        vec = vec/np.linalg.norm(vec, axis = 1)[:, None]
        # print(np.arccos(normals[0, 2]*0.999999))
        ang = np.arccos(normals[:, 2]*0.999999)
        # print(-1*ang[0]*vec[0])
        vec = -1*(vec.T*ang).T
        # print(R.from_rotvec(vec[0]).as_matrix())
        r = R.from_rotvec(vec)
        rotations = r.as_euler('xyz') #, degrees=True)
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(normals, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.bondsetting.coll.objects.link(obj)
        obj.hide_set(True)
        #
        scales = np.concatenate((width.reshape(-1, 1), 
                        width.reshape(-1, 1),
                        scales.reshape(-1, 1)), axis = 1)
        name = '%s_bond_scale'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(scales, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        # obj.data.materials.append(self.material)
        # for name, inputs in battr_inputs.items():
        #     battr = getattr(obj.batoms, name)
        #     for key, value in inputs.items():
        #         setattr(battr, key, value)
        self.bondsetting.coll.objects.link(obj)
        obj.hide_set(True)
        bpy.context.view_layer.update()
        # print('draw bond: {0:10.2f} s'.format(time() - tstart))
    
    def from_bbond(self, label):
        if label not in bpy.data.objects:
            raise Exception("%s is not a object!"%label)
        elif not bpy.data.objects[label].batoms.bbond.flag:
            raise Exception("%s is not Bbond object!"%label)
        obj = bpy.data.objects[label]
        self.species = obj.batoms.bbond.species
        self.label = obj.batoms.bbond.label
        self.element = obj.batoms.bbond.element
        self.species_data = {
            'radius':obj.batoms.bbond.radius,
            'scale':obj.scale,
        } 
    
    @property
    def gnodes(self):
        return self.get_gnodes()
    
    @gnodes.setter
    def gnodes(self, gnodes):
        self.set_gnodes(gnodes)
    
    def get_gnodes(self):
        name = 'GeometryNodes_%s_bond'%self.label
        modifier = self.obj.modifiers.get(name)
        if modifier is None:
            self.build_geometry_node()
        return modifier
    
    def set_gnodes(self, gnodes):
        pass

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        from batoms.butils import get_nodes_by_name
        from batoms.tools import string2Number
        tstart = time()
        name = 'GeometryNodes_%s_bond'%self.label
        modifier = self.obj.modifiers.new(name = name, type = 'NODES')
        modifier.node_group.name = name
        GroupInput = modifier.node_group.nodes.get('Group Input')
        GroupInput.outputs.new(type = 'INT', name = 'species1')
        GroupInput.outputs.new(type = 'INT', name = 'species2')
        GroupInput.outputs.new(type = 'INT', name = 'order')
        GroupInput.outputs.new(type = 'INT', name = 'style')
        GroupInput.outputs.new(type = 'INT', name = 'show')
        GroupInput.outputs.new(type = 'BOOLEAN', name = 'model_style')
        # the above codes not works. maybe bug in blender, 
        # we add this, maybe deleted in the future
        for i in range(1, 7):
            test = get_nodes_by_name(modifier.node_group.nodes, 
                            'BooleanMath_%s'%i,
                            'FunctionNodeCompareFloats')
            modifier.node_group.links.new(GroupInput.outputs[i], test.inputs[0])
        # the Input_1 is Geometry
        modifier['Input_2_use_attribute'] = 1
        modifier['Input_2_attribute_name'] = 'species_index1'
        modifier['Input_3_use_attribute'] = 1
        modifier['Input_3_attribute_name'] = 'species_index2'
        modifier['Input_4_use_attribute'] = 1
        modifier['Input_4_attribute_name'] = 'order'
        modifier['Input_5_use_attribute'] = 1
        modifier['Input_5_attribute_name'] = 'style'
        modifier['Input_6_use_attribute'] = 1
        modifier['Input_6_attribute_name'] = 'show'
        modifier['Input_7_use_attribute'] = 1
        modifier['Input_7_attribute_name'] = 'model_style'
        gn = modifier
        # print(gn.name)
        # print(GroupInput.outputs[:])
        GroupOutput = gn.node_group.nodes.get('Group Output')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        'JoinGeometry_%s'%self.label, 
                        'GeometryNodeJoinGeometry')
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])
        #
        ObjectRotation = get_nodes_by_name(gn.node_group.nodes, 
                        'ObjectInfo_%s_rotation'%(self.label),
                        'GeometryNodeObjectInfo')
        ObjectRotation.inputs['Object'].default_value = self.obj_r
        PositionRotation = get_nodes_by_name(gn.node_group.nodes, 
                        'Position%s_Rotation'%(self.label),
                        'GeometryNodeInputPosition')
        ObjectScale = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s_scale'%(self.label),
                    'GeometryNodeObjectInfo')
        ObjectScale.inputs['Object'].default_value = self.obj_s
        PositionScale = get_nodes_by_name(gn.node_group.nodes, 
                        'Position%s_Scale'%(self.label),
                        'GeometryNodeInputPosition')

        for sp in self.batoms.species:
            CompareSpecies1 = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareSpecies1_%s_%s'%(self.label, sp.name),
                        'FunctionNodeCompareFloats')
            CompareSpecies1.operation = 'EQUAL'
            CompareSpecies1.inputs[1].default_value = string2Number(sp.name)
            CompareSpecies2 = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareSpecies2_%s_%s'%(self.label, sp.name),
                        'FunctionNodeCompareFloats')
            CompareSpecies2.operation = 'EQUAL'
            CompareSpecies2.inputs[1].default_value = string2Number(sp.name)
            gn.node_group.links.new(GroupInput.outputs[1], CompareSpecies1.inputs[0])
            gn.node_group.links.new(GroupInput.outputs[2], CompareSpecies2.inputs[0])
        # order 
        for order in [1, 2, 3]:
            CompareOrder = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareFloats_%s_%s_order'%(self.label, order),
                        'FunctionNodeCompareFloats')
            CompareOrder.operation = 'EQUAL'
            CompareOrder.inputs[1].default_value = order
            gn.node_group.links.new(GroupInput.outputs[3], CompareOrder.inputs[0])
        # style 
        for style in [0, 1, 2]:
            CompareStyle = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareFloats_%s_%s_style'%(self.label, style),
                        'FunctionNodeCompareFloats')
            CompareStyle.operation = 'EQUAL'
            CompareStyle.inputs[1].default_value = style
            gn.node_group.links.new(GroupInput.outputs[4], CompareStyle.inputs[0])
        #
        for sp in self.bondsetting:
            self.add_geometry_node(sp)
        
        print('Build geometry nodes for bonds: %s'%(time() - tstart))

    def add_geometry_node(self, sp):
        from batoms.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        'JoinGeometry_%s'%self.label, 
                        'GeometryNodeJoinGeometry')
        #
        ObjectRotation = get_nodes_by_name(gn.node_group.nodes, 
                        'ObjectInfo_%s_rotation'%(self.label),
                        'GeometryNodeObjectInfo')
        ObjectRotation.inputs['Object'].default_value = self.obj_r
        PositionRotation = get_nodes_by_name(gn.node_group.nodes, 
                        'Position%s_Rotation'%(self.label),
                        'GeometryNodeInputPosition')
        ObjectScale = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s_scale'%(self.label),
                    'GeometryNodeObjectInfo')
        ObjectScale.inputs['Object'].default_value = self.obj_s
        PositionScale = get_nodes_by_name(gn.node_group.nodes, 
                        'Position%s_Scale'%(self.label),
                        'GeometryNodeInputPosition')
        #
        order = sp.order
        style = int(sp.style)
        name = '%s_%s_%s_%s'%(self.label, sp.name, order, style)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                    'InstanceOnPoint_%s'%name,
                    'GeometryNodeInstanceOnPoints')
        ObjectInstancer = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s'%name,
                    'GeometryNodeObjectInfo')
        ObjectInstancer.inputs['Object'].default_value = self.bondsetting.instancers[sp.name]['%s_%s'%(order, style)]
        #
        BoolSpecies = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_species'%name,
                        'FunctionNodeBooleanMath')
        BoolOrder = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_order'%name,
                        'FunctionNodeBooleanMath')
        BoolStyle = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_style'%name,
                        'FunctionNodeBooleanMath')
        BoolModelStyle = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_modelstyle'%name,
                        'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes, 
                    'BooleanMath_%s_show'%name,
                    'FunctionNodeBooleanMath')
        #
        AttributeTransferRotation = get_nodes_by_name(gn.node_group.nodes, 
                    'AttributeTransfer_%s_Rotation'%name,
                    'GeometryNodeAttributeTransfer')
        AttributeTransferRotation.mapping = 'INDEX'
        AttributeTransferRotation.data_type = 'FLOAT_VECTOR'
        AlignEulerToVector = get_nodes_by_name(gn.node_group.nodes, 
                    'AlignEulerToVector_%s'%name,
                    'FunctionNodeAlignEulerToVector')
        AlignEulerToVector.axis = 'Z'
        #
        AttributeTransferScale = get_nodes_by_name(gn.node_group.nodes, 
                    'AttributeTransfer_%s_scale'%name,
                    'GeometryNodeAttributeTransfer')
        AttributeTransferScale.mapping = 'INDEX'
        AttributeTransferScale.data_type = 'FLOAT_VECTOR'
        # BooleanMath.inputs[1].default_value = True
        CompareSpecies1 = get_nodes_by_name(gn.node_group.nodes, 
                    'CompareSpecies1_%s_%s'%(self.label, sp.species1))
        CompareSpecies2 = get_nodes_by_name(gn.node_group.nodes, 
                    'CompareSpecies2_%s_%s'%(self.label, sp.species2))
        CompareOrder = get_nodes_by_name(gn.node_group.nodes, 
                'CompareFloats_%s_%s_order'%(self.label, order))
        CompareStyle = get_nodes_by_name(gn.node_group.nodes, 
                'CompareFloats_%s_%s_style'%(self.label, style))
        #
        gn.node_group.links.new(GroupInput.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[5], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[6], BoolModelStyle.inputs[0])
        gn.node_group.links.new(CompareSpecies1.outputs[0], BoolSpecies.inputs[0])
        gn.node_group.links.new(CompareSpecies2.outputs[0], BoolSpecies.inputs[1])
        gn.node_group.links.new(BoolSpecies.outputs[0], BoolOrder.inputs[0])
        gn.node_group.links.new(CompareOrder.outputs[0], BoolOrder.inputs[1])
        gn.node_group.links.new(BoolOrder.outputs[0], BoolStyle.inputs[0])
        gn.node_group.links.new(CompareStyle.outputs[0], BoolStyle.inputs[1])
        gn.node_group.links.new(BoolStyle.outputs[0], BoolModelStyle.inputs[1])
        gn.node_group.links.new(BoolModelStyle.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectScale.outputs['Geometry'], AttributeTransferScale.inputs['Target'])
        gn.node_group.links.new(PositionScale.outputs[0], AttributeTransferScale.inputs[1])
        gn.node_group.links.new(AttributeTransferScale.outputs[0], InstanceOnPoint.inputs['Scale'])
        #
        gn.node_group.links.new(ObjectRotation.outputs['Geometry'], AttributeTransferRotation.inputs['Target'])
        gn.node_group.links.new(PositionRotation.outputs['Position'], AttributeTransferRotation.inputs['Attribute'])
        gn.node_group.links.new(AttributeTransferRotation.outputs[0], AlignEulerToVector.inputs['Vector'])
        gn.node_group.links.new(AlignEulerToVector.outputs[0], InstanceOnPoint.inputs['Rotation'])
        #
        gn.node_group.links.new(ObjectInstancer.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])

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
        for key in me.attributes.keys():
            att = me.attributes.get(key)
            dtype = att.data_type
            if dtype == 'STRING':
                attributes[key] = np.zeros(nvert, dtype = 'U20')
                for i in range(nvert):
                    attributes[key][i] = att.data[i].value
            elif dtype == 'INT':
                attributes[key] = np.zeros(nvert, dtype = int)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'FLOAT':
                attributes[key] = np.zeros(nvert, dtype = float)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'BOOLEAN':
                attributes[key] = np.zeros(nvert, dtype = bool)
                att.data.foreach_get("value", attributes[key])
            else:
                raise KeyError('%s is not support.'%dtype)
            attributes[key] = np.array(attributes[key])
        return attributes
    def set_attributes(self, attributes):
        tstart = time()
        me = self.obj.data
        for key, data in attributes.items():
            # print(key)
            att = me.attributes.get(key)
            if att is None:
                dtype = type(attributes[key][0])
                if np.issubdtype(dtype, int):
                    dtype = 'INT'
                elif np.issubdtype(dtype, float):
                    dtype = 'FLOAT'
                elif np.issubdtype(dtype, str):
                    dtype = 'STRING'
                else:
                    raise KeyError('%s is not supported.'%dtype)
                att = me.attributes.new(name = key, type = dtype, domain = 'POINT')
            if att.data_type == 'STRING':
                nvert = len(me.vertices)
                for i in range(nvert):
                    att.data[i].value = data[i]
            else:
                att.data.foreach_set("value", data)
        me.update()
        print('set_attributes: %s'%(time() - tstart))

    
    def set_attribute_with_indices(self, name, indices, data):
        data0 = self.attributes[name]
        data0[indices] = data
        self.set_attributes({name: data0})
    
    @property
    def obj_r(self):
        return self.get_obj_r()
    
    def get_obj_r(self):
        name = '%s_bond_normal'%self.label
        obj_r = bpy.data.objects.get(name)
        if obj_r is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_r
    
    @property
    def obj_s(self):
        return self.get_obj_s()
    
    def get_obj_s(self):
        name = '%s_bond_scale'%self.label
        obj_s = bpy.data.objects.get(name)
        if obj_s is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_s


    @property
    def arrays(self):
        return self.get_arrays()
    
    @arrays.setter
    def arrays(self, arrays):
        self.set_arrays(arrays)

    def get_arrays(self, batoms = None, local = False, X = False, sort = True):
        """
        """
        object_mode()
        tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        # radius
        radius = self.radius
        arrays.update({'radius': np.zeros(len(self))})
        for sel, data in radius.items():
            for sp, value in data.items():
                mask = np.where((arrays['species'] == sp) & (arrays['select'] == string2Number(sel)))
                arrays['radius'][mask] = value
        # main elements
        main_elements = self.species.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        arrays.update({'elements': np.array(elements, dtype='U20')})
        print('get_arrays: %s'%(time() - tstart))
        return arrays


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
        return np.array(self.obj.batoms.bbond.width)
    
    @property    
    def segments(self):
        return self.get_segments()
    
    def get_segments(self):
        return self.obj.batoms.bbond.segments
    
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
        name = '%s_bond_center'%self.label
        obj = bpy.data.objects.get(name)
        base_name = 'Basis_%s_center'%self.label
        if obj.data.shape_keys is None:
            obj.shape_key_add(name = base_name)
        elif base_name not in obj.data.shape_keys.key_blocks:
            obj.shape_key_add(name = base_name)
        if only_basis:
            return
        return
        nvert = len(obj.data.vertices)
        mesh_data = bond_source[self.segments]
        for i in range(1, nframe):
            sk = obj.shape_key_add(name = str(i))
            # Use the local position here
            positions = frames[i][:, 0:3]
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

        >>> from batoms.bond import Bbond
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
        
        >>> from batoms.bonds import Bbond
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