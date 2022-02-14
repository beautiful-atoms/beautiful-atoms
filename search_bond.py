from pandas import array
import bpy
from ase import Atoms
import numpy as np
from time import time
from batoms.base import ObjectGN
# from batoms.tools import build_search_bond
from batoms.butils import object_mode
from batoms.tools import number2String, string2Number
from ase.geometry import wrap_positions, complete_cell


shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]

default_attributes = [
        ['atoms_index', 'INT'], 
        ['species_index', 'INT'], 
        ['show', 'BOOLEAN'], 
        ['select', 'INT'],
        ['model_style', 'INT'],
        ['scale', 'FLOAT'], 
        ['radius_style', 'INT'],
        ]

default_search_bond_datas = {
        'atoms_index': np.ones(0, dtype = int),
        'species_index': np.ones(0, dtype = int),
        # 'species': np.ones(0, dtype = 'U4'),
        'positions':np.zeros((0, 3)),
        'scales':np.zeros(0),
        'offsets':np.zeros((0, 3)),
        'model_styles':np.ones(0, dtype = int),
        'radius_styles':np.ones(0, dtype = int),
        'shows':np.ones(0, dtype = int),
        'selects':np.ones(0, dtype = int),
        }

class SearchBond(ObjectGN):
    """SearchBond Class

    """
    def __init__(self, 
                label = None,
                search_bond_datas = None,
                batoms = None,
                 ):
        #
        self.batoms = batoms
        self.label = label
        name = 'search_bond'
        ObjectGN.__init__(self, label, name)
        if search_bond_datas is not None:
            self.build_object(search_bond_datas)
        else:
            self.load(label)

    def build_object(self, search_bond_datas, attributes = {}):
        """
        build child object and add it to main objects.
        """
        tstart = time()
        if len(search_bond_datas['positions'].shape) == 2:
            self._frames = {'positions': np.array([search_bond_datas['positions']]),
                            'offsets': np.array([search_bond_datas['offsets']]),
                        }
            positions = search_bond_datas['positions']
            offsets = search_bond_datas['offsets']
        elif len(search_bond_datas['positions'].shape) == 3:
            self._frames = {'positions': search_bond_datas['positions'],
                            'offsets': search_bond_datas['offsets'],
                            }
            positions = search_bond_datas['positions'][0]
            offsets = search_bond_datas['offsets'][0]
        else:
            raise Exception('Shape of positions is wrong!')
        #
        attributes.update({
                            'atoms_index': search_bond_datas['atoms_index'], 
                            'species_index': search_bond_datas['species_index'], 
                            'show': search_bond_datas['shows'],
                            'model_style': search_bond_datas['model_styles'],
                            'select': search_bond_datas['selects'],
                            'scale': search_bond_datas['scales'],
                            'radius_style': search_bond_datas['radius_styles'],
                            })
        name = '%s_search_bond'%self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(positions, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = 'POINT')
        self.batoms.coll.objects.link(obj)
        #
        name = '%s_search_bond_offset'%self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.batoms.coll.objects.link(obj)
        obj.hide_set(True)
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis = True)
        # print('boundary: build_object: {0:10.2f} s'.format(time() - tstart))
    
    def build_geometry_node(self):
        """
        """
        from batoms.butils import get_nodes_by_name
        name = 'GeometryNodes_%s_search_bond'%self.label
        modifier = self.obj.modifiers.new(name = name, type = 'NODES')
        modifier.node_group.name = name
        #------------------------------------------------------------------
        # select attributes
        GroupInput = modifier.node_group.nodes.get('Group Input')
        # add new output sockets
        for att in default_attributes:
            GroupInput.outputs.new(type = att[1], name = att[0])
        # the above codes not works. maybe bug in blender, 
        # we add this, maybe deleted in the future
        for i in range(1, 7):
            test = get_nodes_by_name(modifier.node_group.nodes, 
                            'BooleanMath_%s'%i,
                            'FunctionNodeCompareFloats')
            modifier.node_group.links.new(GroupInput.outputs[i], test.inputs[0])
        #
        i = 2
        for att in default_attributes:
            modifier['Input_%s_use_attribute'%i] = 1
            modifier['Input_%s_attribute_name'%i] = att[0]
            i += 1
        gn = modifier
        #------------------------------------------------------------------
        GroupOutput = gn.node_group.nodes.get('Group Output')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])
        #------------------------------------------------------------------
        # calculate bond vector, length, rotation based on the index
        # Get four positions from batoms, bond and the second bond for high order bond plane
        ObjectBatoms = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_ObjectBatoms'%self.label,
                    'GeometryNodeObjectInfo')
        ObjectBatoms.inputs['Object'].default_value = self.batoms.obj
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_PositionBatoms'%(self.label),
                        'GeometryNodeInputPosition')
        TransferBatoms = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_TransferBatoms'%(self.label),
                    'GeometryNodeAttributeTransfer')
        TransferBatoms.mapping = 'INDEX'
        TransferBatoms.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectBatoms.outputs['Geometry'], TransferBatoms.inputs['Target'])
        gn.node_group.links.new(PositionBatoms.outputs['Position'], TransferBatoms.inputs['Attribute'])
        gn.node_group.links.new(GroupInput.outputs[1], TransferBatoms.inputs['Index'])
        #------------------------------------------------------------------
        # add positions with offsets
        # transfer offsets from object self.obj_o
        ObjectOffsets = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_ObjectOffsets'%(self.label),
                        'GeometryNodeObjectInfo')
        ObjectOffsets.inputs['Object'].default_value = self.obj_o
        PositionOffsets = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_PositionOffsets'%(self.label),
                        'GeometryNodeInputPosition')
        TransferOffsets = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_TransferOffsets'%self.label,
                    'GeometryNodeAttributeTransfer')
        TransferOffsets.mapping = 'INDEX'
        TransferOffsets.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectOffsets.outputs['Geometry'], TransferOffsets.inputs['Target'])
        gn.node_group.links.new(PositionOffsets.outputs['Position'], TransferOffsets.inputs['Attribute'])
        OffsetNode = self.vectorDotMatrix(gn, TransferOffsets, self.batoms.cell, '')
        # we need one add operation to get the positions with offset
        VectorAdd = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorAdd'%(self.label),
                    'ShaderNodeVectorMath')
        VectorAdd.operation = 'ADD'
        gn.node_group.links.new(TransferBatoms.outputs[0], VectorAdd.inputs[0])
        gn.node_group.links.new(OffsetNode.outputs[0], VectorAdd.inputs[1])
        # set positions
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label, 
                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(VectorAdd.outputs[0], SetPosition.inputs['Position'])
        
    def add_geometry_node(self, spname, selname):
        """
        """
        from batoms.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        CompareSelect = get_nodes_by_name(gn.node_group.nodes, 
                    'select_%s_%s'%(self.label, selname),
                    'FunctionNodeCompareFloats')
        CompareSelect.operation = 'EQUAL'
        # CompareSelect.data_type = 'INT'
        CompareSelect.inputs[1].default_value = string2Number(selname)
        gn.node_group.links.new(GroupInput.outputs[4], CompareSelect.inputs[0])
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes, 
                    'CompareFloats_%s_%s'%(self.label, spname),
                    'FunctionNodeCompareFloats')
        CompareSpecies.operation = 'EQUAL'
        # CompareSpecies.data_type = 'INT'
        CompareSpecies.inputs[1].default_value = string2Number(spname)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                    'InstanceOnPoint_%s_%s_%s'%(self.label, selname, spname), 
                    'GeometryNodeInstanceOnPoints')
        ObjectInfo = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s_%s_%s'%(self.label, selname, spname),
                    'GeometryNodeObjectInfo')
        ObjectInfo.inputs['Object'].default_value = self.batoms.species.instancers[selname][spname]
        #
        BoolSelectSpecies = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_%s_%s_0'%(self.label, selname, spname),
                        'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes, 
                    'BooleanMath_%s_%s_%s_1'%(self.label, selname, spname),
                    'FunctionNodeBooleanMath')
        #
        gn.node_group.links.new(SetPosition.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[2], CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[6], InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSelect.outputs[0], BoolSelectSpecies.inputs[0])
        gn.node_group.links.new(CompareSpecies.outputs[0], BoolSelectSpecies.inputs[1])
        gn.node_group.links.new(BoolSelectSpecies.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])        
    
    @property
    def obj_o(self):
        return self.get_obj_o()
    
    def get_obj_o(self):
        name = '%s_search_bond_offset'%self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_o
    
    def get_search_bond(self):
        boundary = np.array(self.batoms.coll.batoms.boundary)
        return boundary.reshape(3, -1)

    def set_arrays(self, arrays):
        """
        """
        if len(arrays['positions']) == 0:
            return
        attributes = self.attributes
        # same length
        if len(arrays['positions']) == len(attributes['show']):
            self.positions = arrays['positions']
            self.offsets = arrays['offsets']
            species_index = [string2Number(sp) for sp in arrays['species']]
            self.set_attributes({'species_index': species_index,
                                'scale': arrays['scales'],
                                'show': arrays['shows'],
                                })
        else:
            # add or remove vertices
            self.build_object(arrays)
            species = np.unique(arrays['species'])
            for sp in species:
                self.add_geometry_node(sp, 'sel0')

    def get_arrays(self):
        """
        """
        object_mode()
        tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        arrays.update({'offsets': self.offsets})
        # radius
        radius = self.batoms.radius
        arrays.update({'radius': np.zeros(len(self))})
        species = np.array([number2String(i) for i in arrays['species_index']], dtype = 'U20')
        arrays['species'] = species
        for sel, data in radius.items():
            for sp, value in data.items():
                mask = np.where((arrays['species'] == sp) & (arrays['select'] == string2Number(sel)))
                arrays['radius'][mask] = value
        # size
        arrays['size'] = arrays['radius']*arrays['scale']
        # main elements
        main_elements = self.batoms.species.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        arrays.update({'elements': np.array(elements, dtype='U20')})
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    @property
    def search_bond_data(self):
        return self.get_search_bond_data()

    def get_search_bond_data(self, include_batoms = False):
        """
        using foreach_get and foreach_set to improve performance.
        """
        arrays = self.arrays
        search_bond_data = {'positions': arrays['positions'],
                        'species': arrays['species'],
                        'indices': arrays['atoms_index'],
                        'offsets': arrays['offsets']}
        return search_bond_data

    @property
    def offsets(self):
        return self.get_offsets()
    
    def get_offsets(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        offsets = np.empty(n*3, dtype = int)
        self.obj_o.data.vertices.foreach_get('co', offsets)  
        return offsets.reshape((n, 3))
    
    @offsets.setter
    def offsets(self, offsets):
        self.set_offsets(offsets)
    
    def set_offsets(self, offsets):
        """
        Set global offsets to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        n = len(self.obj_o.data.vertices)
        if len(offsets) != n:
            raise ValueError('offsets has wrong shape %s != %s.' %
                                (len(offsets), n))
        offsets = offsets.reshape((n*3, 1))
        self.obj_o.data.shape_keys.key_blocks[0].data.foreach_set('co', offsets)
        self.obj_o.data.update()
        bpy.context.view_layer.objects.active = self.obj_o
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
    
    @property
    def bondlists(self):
        return self.get_bondlists()
    
    def get_bondlists(self):
        try:
            bondlists = self.batoms.bonds.arrays
        except:
            bondlists = None
        return bondlists

    def get_frames(self):
        """
        """
        frames = {}
        frames['positions'] = self.get_obj_frames(self.obj)
        frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames

    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0 : return
        name = '%s_search_bond'%(self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['positions'])
        #
        name = '%s_search_bond_offset'%(self.label)
        obj = self.obj_o
        self.set_obj_frames(name, obj, frames['offsets'])