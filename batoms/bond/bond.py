"""Definition of the Bond class.

This module defines the Bond object in the Batoms package.

"""

import bpy
import bmesh
from time import time
from batoms.attribute import Attributes
from batoms.utils.butils import object_mode, compareNodeType, get_nodes_by_name
from batoms.utils import string2Number, number2String
import numpy as np
from batoms.base.object import ObjectGN
from batoms.base.collection import BaseCollection
from .setting import BondSettings
from batoms.bond.search_bond import SearchBond, default_search_bond_datas
# from pprint import pprint
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

default_attributes = [
    {"name": 'atoms_index1', "data_type": 'INT', "dimension": 0},
    {"name": 'atoms_index2', "data_type": 'INT', "dimension": 0},
    {"name": 'atoms_index3', "data_type": 'INT', "dimension": 0},
    {"name": 'atoms_index4', "data_type": 'INT', "dimension": 0},
    {"name": 'species_index1', "data_type": 'INT', "dimension": 0},
    {"name": 'species_index2', "data_type": 'INT', "dimension": 0},
    {"name": 'order', "data_type": 'INT', "dimension": 0},
    {"name": 'style', "data_type": 'INT', "dimension": 0},
    {"name": 'show', "data_type": 'INT', "dimension": 0},
    {"name": 'model_style', "data_type": 'INT', "dimension": 0},
    {"name": 'polyhedra', "data_type": 'INT', "dimension": 0},
    {"name": 'second_bond', "data_type": 'INT', "dimension": 0},
]

default_GroupInput = [
    ['atoms_index1', 'NodeSocketInt'],
    ['atoms_index2', 'NodeSocketInt'],
    ['atoms_index3', 'NodeSocketInt'],
    ['atoms_index4', 'NodeSocketInt'],
    ['species_index1', 'NodeSocketInt'],
    ['species_index2', 'NodeSocketInt'],
    ['order', 'NodeSocketInt'],
    ['style', 'NodeSocketInt'],
    ['show', 'NodeSocketBool'],
    ['model_style', 'NodeSocketInt'],
    ['second_bond', 'NodeSocketInt'],
]

default_bond_datas = {
    'atoms_index1': np.ones(0, dtype=int),
    'atoms_index2': np.ones(0, dtype=int),
    'atoms_index3': np.ones(0, dtype=int),
    'atoms_index4': np.ones(0, dtype=int),
    'species_index1': np.ones(0, dtype=int),
    'species_index2': np.ones(0, dtype=int),
    'centers': np.zeros((1, 0, 3)),
    # 'vectors':np.zeros((0, 3)),
    'offsets1': np.zeros((0, 3)),
    'offsets2': np.zeros((0, 3)),
    'offsets3': np.zeros((0, 3)),
    'offsets4': np.zeros((0, 3)),
    # 'eulers':np.eye(3),
    # 'lengths':np.zeros((0, 3)),
    # 'widths': np.ones(0, dtype=float),
    'show': np.zeros(0, dtype=int),
    'order': np.zeros(0, dtype=int),
    'style': np.zeros(0, dtype=int),
    'model_style': np.ones(0, dtype=int),
    'polyhedra': np.ones(0, dtype=float),
    'second_bond': np.ones(0, dtype=int),
}


class Bond(BaseCollection, ObjectGN):
    """Bbond Class

    A Bbond object is linked to this main collection in Blender.

    Parameters:


    """

    def __init__(self,
                 label=None,
                 batoms=None,
                 bond_datas=None,
                 location=np.array([0, 0, 0]),
                 ):
        #
        self.batoms = batoms
        self.label = label
        name = 'bond'
        ObjectGN.__init__(self, label, name)
        BaseCollection.__init__(self, coll_name=label)
        if bond_datas is not None:
            # overwrite old object
            self.build_object(bond_datas)
            self.settings = BondSettings(self.label, batoms=batoms, bonds=self)
            self.update_geometry_nodes()
        elif self.loadable():
            # load old object
            self.settings = BondSettings(self.label, batoms=batoms, bonds=self)
            self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        else:
            # create empty new object
            self.build_object(default_bond_datas)
            self.settings = BondSettings(self.label, batoms=batoms, bonds=self)
            self.update_geometry_nodes()
        self._search_bond = None

    def loadable(self):
        """Check loadable or not
        """
        # object exist
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            return False
        # batoms exist, and flag is True
        coll = bpy.data.collections.get(self.label)
        if coll is None:
            return False
        return coll.Bbond.flag

    def build_object(self, bond_datas, attributes={}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        tstart = time()
        if len(bond_datas['centers'].shape) == 2:
            self._frames = {'centers': np.array([bond_datas['centers']]),
                            'offsets1': np.array([bond_datas['offsets1']]),
                            'offsets2': np.array([bond_datas['offsets2']]),
                            'offsets3': np.array([bond_datas['offsets3']]),
                            'offsets4': np.array([bond_datas['offsets4']]),
                            }
            centers = bond_datas['centers']
            offsets1 = bond_datas['offsets1']
            offsets2 = bond_datas['offsets2']
            offsets3 = bond_datas['offsets3']
            offsets4 = bond_datas['offsets4']
        elif len(bond_datas['centers'].shape) == 3:
            self._frames = {'centers': bond_datas['centers'],
                            'offsets1': bond_datas['offsets1'],
                            'offsets2': bond_datas['offsets2'],
                            'offsets3': bond_datas['offsets3'],
                            'offsets4': bond_datas['offsets4'],
                            }
            centers = bond_datas['centers'][0]
        else:
            raise Exception('Shape of centers is wrong!')
        offsets1 = bond_datas['offsets1']
        offsets2 = bond_datas['offsets2']
        offsets3 = bond_datas['offsets3']
        offsets4 = bond_datas['offsets4']
        attributes.update({
            'atoms_index1': bond_datas['atoms_index1'],
            'atoms_index2': bond_datas['atoms_index2'],
            'atoms_index3': bond_datas['atoms_index3'],
            'atoms_index4': bond_datas['atoms_index4'],
            'species_index1': bond_datas['species_index1'],
            'species_index2': bond_datas['species_index2'],
            'show': bond_datas['show'],
            'model_style': bond_datas['model_style'],
            'style': bond_datas['style'],
            'order': bond_datas['order'],
            'second_bond': bond_datas['second_bond'],
        })
        name = self.obj_name
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(centers, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        # Add attributes
        self._attributes.add(default_attributes)
        self.coll.objects.link(obj)
        obj.parent = self.batoms.obj
        obj.batoms.type = 'BOND'
        obj.batoms.label = self.batoms.label
        obj.Bbond.label = self.batoms.label
        self.batoms.coll.Bbond.flag = True
        # The reason we use object to store offset instead of using vector attribute
        # of the mesh, is that we need shape key for offsets too.
        offsets = [offsets1, offsets2, offsets3, offsets4]
        for i in range(4):
            name = '%s_bond_offset%s' % (self.label, i)
            self.delete_obj(name)
            mesh = bpy.data.meshes.new(name)
            mesh.from_pydata(offsets[i], [], [])
            mesh.update()
            obj = bpy.data.objects.new(name, mesh)
            self.coll.objects.link(obj)
            obj.hide_set(True)
            obj.parent = self.obj
        #
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis=False)
        #
        logger.debug('bonds: build_object: {0:10.2f} s'.format(time() - tstart))

    def build_geometry_node(self):
        """
        todo: add width to nodes
        """
        from batoms.utils.butils import get_nodes_by_name, compareNodeType, \
                            build_modifier
        from batoms.utils import string2Number
        tstart = time()
        name = 'GeometryNodes_%s_bond' % self.label
        modifier = build_modifier(self.obj, name)
        # ------------------------------------------------------------------
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = modifier.node_group.nodes[1]
        # Blender 3.1.2 crashed on Win10
        # separate this from the aboved for loop avoid the crash
        # I don't know why, could be a bug of Blender
        for att in default_GroupInput:
            inputs.new(att[1], att[0])
            id = inputs[att[0]].identifier
            modifier['%s_use_attribute' % id] = True
            modifier['%s_attribute_name' % id] = att[0]
        gn = modifier
        # ------------------------------------------------------------------
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])
        # ------------------------------------------------------------------
        # calculate bond vector, length, rotation based on the index
        # Get four positions from batoms, bond and the second bond
        # for high order bond plane
        ObjectBatoms = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_ObjectBatoms' % self.label,
                                         'GeometryNodeObjectInfo')
        ObjectBatoms.inputs['Object'].default_value = self.batoms.obj
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_PositionBatoms' % (self.label),
                                           'GeometryNodeInputPosition')
        TransferBatoms = []
        for i in range(4):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                '%s_TransferBatoms%s' % (self.label, i),
                                'GeometryNodeSampleIndex')
            tmp.data_type = 'FLOAT_VECTOR'
            TransferBatoms.append(tmp)
        for i in range(4):
            gn.node_group.links.new(ObjectBatoms.outputs['Geometry'],
                                    TransferBatoms[i].inputs[0])
            gn.node_group.links.new(PositionBatoms.outputs['Position'],
                                    TransferBatoms[i].inputs[3])
            gn.node_group.links.new(GroupInput.outputs[i + 1],
                                    TransferBatoms[i].inputs['Index'])
        # ------------------------------------------------------------------
        # add positions with offsets
        # transfer offsets from object self.obj_o
        ObjectOffsets = []
        PositionOffsets = []
        TransferOffsets = []
        InputIndex = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_InputIndex' % self.label,
                                        'GeometryNodeInputIndex')
        for i in range(4):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_ObjectOffsets%s' % (self.label, i),
                                    'GeometryNodeObjectInfo')
            tmp.inputs['Object'].default_value = self.obj_o[i]
            ObjectOffsets.append(tmp)
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_PositionOffsets%s' % (self.label, i),
                                    'GeometryNodeInputPosition')
            PositionOffsets.append(tmp)
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                '%s_TransferOffsets%i' % (self.label, i),
                                'GeometryNodeSampleIndex')
            tmp.data_type = 'FLOAT_VECTOR'
            TransferOffsets.append(tmp)
            gn.node_group.links.new(ObjectOffsets[i].outputs['Geometry'],
                                    TransferOffsets[i].inputs[0])
            gn.node_group.links.new(PositionOffsets[i].outputs['Position'],
                                    TransferOffsets[i].inputs[3])
            gn.node_group.links.new(InputIndex.outputs[0],
                                    TransferOffsets[i].inputs["Index"])
        # we need five add operations
        # four: Get the positions with offset for four atoms
        # one: Get center = (positions1 + positions2)/2
        VectorAdd = []
        for i in range(5):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_VectorAdd%s' % (self.label, i),
                                    'ShaderNodeVectorMath')
            tmp.operation = 'ADD'
            VectorAdd.append(tmp)
        for i in range(4):
            gn.node_group.links.new(TransferBatoms[i].outputs[2],
                                    VectorAdd[i].inputs[0])
            gn.node_group.links.new(TransferOffsets[i].outputs[2],
                                    VectorAdd[i].inputs[1])
        #
        # divide by 2 to get the center
        VectorDivide = get_nodes_by_name(gn.node_group.nodes,
                                         'VectorDivide_%s' % self.label,
                                         'ShaderNodeVectorMath')
        VectorDivide.operation = 'DIVIDE'
        VectorDivide.inputs[1].default_value = (2, 2, 2)
        gn.node_group.links.new(
            VectorAdd[0].outputs[0], VectorAdd[4].inputs[0])
        gn.node_group.links.new(
            VectorAdd[1].outputs[0], VectorAdd[4].inputs[1])
        gn.node_group.links.new(
            VectorAdd[4].outputs[0], VectorDivide.inputs[0])
        # set center of the bond
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label,
                                        'GeometryNodeSetPosition')
        gn.node_group.links.new(
            GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(
            VectorDivide.outputs[0], SetPosition.inputs['Position'])
        # get the vector for the bond and the length
        # also the vector for the second bond
        VectorSubtract = []
        for i in range(2):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_VectorSubtract%s' % (self.label, i),
                                    'ShaderNodeVectorMath')
            tmp.operation = 'SUBTRACT'
            VectorSubtract.append(tmp)
        VectorLength = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_VectorLength' % self.label,
                                         'ShaderNodeVectorMath')
        VectorLength.operation = 'LENGTH'
        VectorCross0 = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_VectorCross0' % self.label,
                                         'ShaderNodeVectorMath')
        VectorCross0.operation = 'CROSS_PRODUCT'
        gn.node_group.links.new(
            VectorAdd[0].outputs[0], VectorSubtract[0].inputs[0])
        gn.node_group.links.new(
            VectorAdd[1].outputs[0], VectorSubtract[0].inputs[1])
        gn.node_group.links.new(
            VectorAdd[2].outputs[0], VectorSubtract[1].inputs[0])
        gn.node_group.links.new(
            VectorAdd[3].outputs[0], VectorSubtract[1].inputs[1])
        # calc the bond length, use it to build scale
        gn.node_group.links.new(
            VectorSubtract[0].outputs[0], VectorLength.inputs[0])
        #
        CombineXYZ = get_nodes_by_name(gn.node_group.nodes,
                                       '%s_CombineXYZ' % self.label,
                                       'ShaderNodeCombineXYZ')
        CombineXYZ.inputs[0].default_value = 1
        CombineXYZ.inputs[1].default_value = 1
        gn.node_group.links.new(
            VectorLength.outputs['Value'], CombineXYZ.inputs['Z'])
        # cross for rotation, for high order bond
        gn.node_group.links.new(
            VectorSubtract[0].outputs[0], VectorCross0.inputs[0])
        gn.node_group.links.new(
            VectorSubtract[1].outputs[0], VectorCross0.inputs[1])
        # get Euler for rotation
        # we need align two vectors to fix a plane
        # we build the instancer by fix bond diection to Z, and
        # high order bond shift to X
        # thus the the normal of high order bond plane is Y
        AlignEuler = []
        for i in range(2):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_AlignEuler%s' % (self.label, i),
                                    'FunctionNodeAlignEulerToVector')
            AlignEuler.append(tmp)
        AlignEuler[0].axis = 'Z'
        AlignEuler[1].axis = 'Y'
        # We should fix Z when align Y
        AlignEuler[1].pivot_axis = 'Z'
        gn.node_group.links.new(
            VectorSubtract[0].outputs[0], AlignEuler[0].inputs['Vector'])
        gn.node_group.links.new(
            AlignEuler[0].outputs[0], AlignEuler[1].inputs['Rotation'])
        gn.node_group.links.new(
            VectorCross0.outputs[0], AlignEuler[1].inputs['Vector'])
        # find bond kinds by the names of species
        for sp in self.batoms.species:
            # we need two compares for one species,
            # because we have two sockets: species_index1 and species_index2
            CompareSpecies = []
            for i in range(2):
                tmp = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_CompareSpecies_%s_%s' % (
                                            self.label, sp.name, i),
                                        compareNodeType)
                tmp.operation = 'EQUAL'
                tmp.inputs[1].default_value = string2Number(sp.name)
                CompareSpecies.append(tmp)
            gn.node_group.links.new(
                GroupInput.outputs[5], CompareSpecies[0].inputs[0])
            gn.node_group.links.new(
                GroupInput.outputs[6], CompareSpecies[1].inputs[0])
        # order
        for order in [1, 2, 3]:
            CompareOrder = get_nodes_by_name(gn.node_group.nodes,
                                             'CompareFloats_%s_%s_order' % (
                                                 self.label, order),
                                             compareNodeType)
            CompareOrder.operation = 'EQUAL'
            CompareOrder.inputs[1].default_value = order
            gn.node_group.links.new(
                GroupInput.outputs[7], CompareOrder.inputs[0])
        # style
        for style in [0, 1, 2, 3]:
            CompareStyle = get_nodes_by_name(gn.node_group.nodes,
                                             'CompareFloats_%s_%s_style' % (
                                                 self.label, style),
                                             compareNodeType)
            CompareStyle.operation = 'EQUAL'
            CompareStyle.inputs[1].default_value = style
            gn.node_group.links.new(
                GroupInput.outputs[8], CompareStyle.inputs[0])
        logger.debug('Build geometry nodes for bonds: %s' % (time() - tstart))
        #

    def update_geometry_nodes(self):
        tstart = time()
        # for sp in self.settings:
        # self.add_geometry_node(sp.as_dict())
        # only build pair in bondlists
        bondlists = self.arrays
        if len(bondlists['atoms_index1']) == 0:
            return
        pairs = np.concatenate((
            bondlists['species_index1'].reshape(-1, 1),
            bondlists['species_index2'].reshape(-1, 1),
            bondlists['order'].reshape(-1, 1),
            bondlists['style'].reshape(-1, 1)), axis=1)
        pairs = np.unique(pairs, axis=0)
        pairs = pairs.reshape(-1, 4)
        for pair in pairs:
            sp1 = number2String(pair[0])
            sp2 = number2String(pair[1])
            order = pair[2]
            style = pair[3]
            sp = self.settings['%s-%s' % (sp1, sp2)]
            self.add_geometry_node(sp.as_dict(), order, style)

        logger.debug('Update geometry nodes for bonds: %s' % (time() - tstart))

    def add_geometry_node(self, sp, order=None, style=None):
        """
        add geometry node for each bond pair
        """
        from batoms.utils.butils import get_nodes_by_name
        # tstart = time()
        if not order:
            order = sp["order"]
        if not style:
            style = int(sp["style"])
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label)
        #
        name = 'bond_%s_%s_%s_%s' % (self.label, sp["name"], order, style)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                                            'InstanceOnPoint_%s' % name,
                                            'GeometryNodeInstanceOnPoints')
        ObjectInstancer = get_nodes_by_name(gn.node_group.nodes,
                                            'ObjectInfo_%s' % name,
                                            'GeometryNodeObjectInfo')
        ObjectInstancer.inputs['Object'].default_value = \
            self.settings.instancers[sp["name"]]['%s_%s' % (order, style)]
        #
        BoolSpecies = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_BooleanMath_species' % name,
                                        'FunctionNodeBooleanMath')
        BoolOrder = get_nodes_by_name(gn.node_group.nodes,
                                      '%s_BooleanMath_order' % name,
                                      'FunctionNodeBooleanMath')
        BoolStyle = get_nodes_by_name(gn.node_group.nodes,
                                      '%s_BooleanMath_style' % name,
                                      'FunctionNodeBooleanMath')
        BoolModelStyle = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_BooleanMath_modelstyle' % name,
                                           'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes,
                                     '%s_BooleanMath_show' % name,
                                     'FunctionNodeBooleanMath')
        BoolBondLength = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_BoolBondLength' % name,
                                           'FunctionNodeBooleanMath')
        # bondlength larger than max will not show
        LessBondLength = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_LessBondLength' % name,
                                           'ShaderNodeMath')
        LessBondLength.operation = 'LESS_THAN'
        LessBondLength.inputs[1].default_value = sp['max']
        VectorLength = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_VectorLength' % self.label)
        #
        CompareSpecies0 = get_nodes_by_name(
            gn.node_group.nodes,
            '%s_CompareSpecies_%s_0' % (self.label,
                                        sp["species1"]),
            compareNodeType)
        CompareSpecies1 = get_nodes_by_name(
            gn.node_group.nodes,
            '%s_CompareSpecies_%s_1' % (self.label,
                                        sp["species2"]),
            compareNodeType)
        CompareOrder = get_nodes_by_name(
            gn.node_group.nodes,
            'CompareFloats_%s_%s_order' % (self.label,
                                           order),
            compareNodeType)
        CompareStyle = get_nodes_by_name(
            gn.node_group.nodes,
            'CompareFloats_%s_%s_style' % (self.label,
                                           style),
            compareNodeType)
        #
        CombineXYZ = get_nodes_by_name(
            gn.node_group.nodes,
            '%s_CombineXYZ' % self.label,
            'ShaderNodeCombineXYZ')
        AlignEuler1 = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_AlignEuler1' % self.label,
                                        'FunctionNodeAlignEulerToVector')
        gn.node_group.links.new(
            SetPosition.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[9], BoolShow.inputs[0])
        gn.node_group.links.new(
            GroupInput.outputs[10], BoolModelStyle.inputs[0])
        gn.node_group.links.new(
            CompareSpecies0.outputs[0], BoolSpecies.inputs[0])
        gn.node_group.links.new(
            CompareSpecies1.outputs[0], BoolSpecies.inputs[1])
        gn.node_group.links.new(BoolSpecies.outputs[0], BoolOrder.inputs[0])
        gn.node_group.links.new(CompareOrder.outputs[0], BoolOrder.inputs[1])
        gn.node_group.links.new(BoolOrder.outputs[0], BoolStyle.inputs[0])
        gn.node_group.links.new(CompareStyle.outputs[0], BoolStyle.inputs[1])
        gn.node_group.links.new(BoolStyle.outputs[0], BoolModelStyle.inputs[1])
        gn.node_group.links.new(BoolModelStyle.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(VectorLength.outputs['Value'],
                                LessBondLength.inputs[0])
        gn.node_group.links.new(BoolShow.outputs['Boolean'],
                                BoolBondLength.inputs[0])
        gn.node_group.links.new(LessBondLength.outputs[0],
                                BoolBondLength.inputs[1])
        gn.node_group.links.new(BoolBondLength.outputs['Boolean'],
                                InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(CombineXYZ.outputs[0],
                                InstanceOnPoint.inputs['Scale'])
        #
        gn.node_group.links.new(AlignEuler1.outputs[0],
                                InstanceOnPoint.inputs['Rotation'])
        #
        gn.node_group.links.new(ObjectInstancer.outputs['Geometry'],
                                InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'],
                                JoinGeometry.inputs['Geometry'])
        # print('Add geometry nodes for bonds: %s'%(time() - tstart))


    def update_geometry_node_species(self):
        # find bond kinds by the names of species
        tstart = time()
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        for sp in self.batoms.species:
            # we need two compares for one species,
            # because we have two sockets: species_index1 and species_index2
            CompareSpecies = []
            for i in range(2):
                tmp = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_CompareSpecies_%s_%s' % (
                                            self.label, sp.name, i),
                                        compareNodeType)
                tmp.operation = 'EQUAL'
                tmp.inputs[1].default_value = string2Number(sp.name)
                CompareSpecies.append(tmp)
            gn.node_group.links.new(
                GroupInput.outputs[5], CompareSpecies[0].inputs[0])
            gn.node_group.links.new(
                GroupInput.outputs[6], CompareSpecies[1].inputs[0])
        logger.debug('UPdate geometry nodes species for bonds: %s'%(time() - tstart))

    def update_geometry_node_instancer(self):
        """
        Make sure all pair has a geometry node flow
        and the instancer and material are updated.
        """
        tstart = time()
        #
        bondlists = self.arrays
        if len(bondlists['atoms_index1']) == 0:
            return
        pairs = np.concatenate((
            bondlists['species_index1'].reshape(-1, 1),
            bondlists['species_index2'].reshape(-1, 1),
            bondlists['order'].reshape(-1, 1),
            bondlists['style'].reshape(-1, 1)), axis=1)
        pairs = np.unique(pairs, axis=0)
        pairs = pairs.reshape(-1, 4)
        for pair in pairs:
            sp1 = number2String(pair[0])
            sp2 = number2String(pair[1])
            order = pair[2]
            style = pair[3]
            order_style = '%s_%s' % (order, style)
            sp = self.settings['%s-%s' % (sp1, sp2)]
            sp = sp.as_dict()
            # update geometry node
            if self.settings.instancers[sp["name"]][order_style] is None:
                self.settings.build_instancer(sp, order, style)
                self.add_geometry_node(sp)
            #always update materials
            self.settings.build_materials(sp, material_style=sp['material_style'])
            self.settings.assign_materials(sp, sp["order"], sp["style"])
            # compare radius
            if not np.isclose(sp["width"],
                    self.settings.instancers[sp["name"]][order_style].Bbond.width):
                self.settings.build_instancer(sp)
                self.settings.assign_materials(sp, sp["order"], sp["style"])
            # update  instancers
            name = 'bond_%s_%s_%s_%s' % (self.label, sp["name"],
                                         sp["order"], sp["style"])
            ObjectInstancer = get_nodes_by_name(self.gnodes.node_group.nodes,
                                                'ObjectInfo_%s' % name,
                                                'GeometryNodeObjectInfo')
            ObjectInstancer.inputs['Object'].default_value = \
                self.settings.instancers[sp["name"]][order_style]
        logger.debug('update bond instancer: %s' % (time() - tstart))

    def update(self, bondlists = None, orders = None):
        """
        Draw bonds.
        calculate bond in all farmes, and merge all bondlists
        """

        object_mode()
        # clean_coll_objects(self.coll, 'bond')
        frames = self.batoms.get_frames()
        arrays = self.batoms.arrays
        boundary_data = self.batoms.boundary.boundary_data
        show = arrays['show'].astype(bool)
        species = arrays['species'][show]
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        bond_datas = {}
        tstart = time()
        setting = self.settings.as_dict()
        if bondlists is None:
            for f in range(nframe):
                # print('update bond: ', f)
                positions = frames[f, show, :]
                # build bondlist for unit cell
                bondlist, bonddatas, peciesBondDatas, molPeciesDatas = \
                    self.build_bondlists(species, positions, self.batoms.cell,
                                        self.batoms.pbc, setting)
                # build bondlist for boundary atoms
                # for molecule with cell == [0, 0, 0], skip
                if boundary_data is not None:
                    bondlist = self.build_bondlists_with_boundary(
                        boundary_data, bondlist, bonddatas,
                        peciesBondDatas,
                        molPeciesDatas)
                    bondlist = self.check_boundary(bondlist)
                # search molecule
                # search bond
                if self.show_search:
                    self.search_bond.hide = False
                    self.search_bond.update(
                        bondlist, peciesBondDatas, molPeciesDatas,
                        arrays, self.batoms.cell)
                if f == 0:
                    bondlists = bondlist
                else:
                    bondlists = np.append(bondlists, bondlist, axis=0)
            bondlists = np.unique(bondlists, axis=0)
        bond_datas = self.calc_bond_data(species, frames[:, show, :],
                                         self.batoms.cell, bondlists,
                                         self.settings,
                                         arrays['model_style'][show])
        if orders:
            bond_datas.update({"order": orders})
        if len(bond_datas) == 0:
            return
        self.set_arrays(bond_datas)
        logger.debug('draw bond: {0:10.2f} s'.format(time() - tstart))

    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            self.build_object(default_bond_datas)
            obj = bpy.data.objects.get(self.obj_name)
        return obj

    @property
    def obj_o(self):
        return self.get_obj_o()

    def get_obj_o(self):
        objs = []
        for i in range(4):
            name = '%s_bond_offset%s' % (self.label, i)
            obj_o = bpy.data.objects.get(name)
            if obj_o is None:
                raise KeyError('%s object is not exist.' % name)
            else:
                objs.append(obj_o)
        return objs

    def get_arrays(self):
        """
        """
        object_mode()
        # tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions,
                       'offsets1': self.offsets[0],
                       'offsets2': self.offsets[1],
                       'offsets3': self.offsets[2],
                       'offsets4': self.offsets[3],
                       })
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    def set_arrays(self, arrays):
        """
        """
        attributes = self.attributes
        # same length
        dnvert = len(arrays['atoms_index1']) - len(attributes['atoms_index1'])
        if dnvert > 0:
            # add
            objs = self.obj_o
            for obj in objs:
                self.add_vertices_bmesh(dnvert, obj)
            self.add_vertices_bmesh(dnvert)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
            objs = self.obj_o
            for obj in objs:
                self.delete_vertices_bmesh(range(-dnvert), obj)
        self.set_frames(arrays)
        self.offsets = [arrays['offsets1'], arrays['offsets2'],
                        arrays['offsets3'], arrays['offsets4']]
        arrays = {key:value for key, value in arrays.items()
                if key not in ["centers", "offsets1", "offsets2",
                    "offsets3", "offsets4"]}
        self.set_attributes(arrays)
        self.update_geometry_node_species()
        self.update_geometry_node_instancer()
        self.update_geometry_nodes()

    @property
    def offsets(self):
        return self.get_offsets()

    def get_offsets(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        offsets = []
        objs = self.obj_o
        for obj in objs:
            positions = np.empty(n*3, dtype=np.float64)
            obj.data.vertices.foreach_get('co', positions)
            offsets.append(positions.reshape((n, 3)))
        return offsets

    @offsets.setter
    def offsets(self, offsets):
        self.set_offsets(offsets)

    def set_offsets(self, offsets):
        """
        Set global offsets to local vertices
        """
        object_mode()
        objs = self.obj_o
        n = len(objs[0].data.vertices)
        if len(offsets[0]) != n:
            raise ValueError('offsets has wrong shape %s != %s.' %
                             (len(offsets[0]), n))
        if n == 0:
            return
        for i in range(4):
            vertices = offsets[i].reshape((n*3, 1))
            objs[i].data.vertices.foreach_set('co', vertices)
            objs[i].data.update()

    def get_frames(self):
        """
        """
        frames = {}
        frames['centers'] = self.get_obj_frames(self.obj)
        # frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames['centers'])
        if nframe == 0:
            return
        name = '%s_bond' % (self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['centers'])
        #
        # name = '%s_bond_offset'%(self.label)
        # obj = self.obj_o
        # self.set_obj_frames(name, obj, frames['offsets'])

    @property
    def bondlists(self):
        return self.get_bondlists()

    def get_bondlists(self):
        """
        """
        object_mode()
        # tstart = time()
        arrays = self.arrays
        i = arrays['atoms_index1'].reshape(-1, 1)
        j = arrays['atoms_index2'].reshape(-1, 1)
        p = arrays['polyhedra'].reshape(-1, 1)
        offsets1 = arrays['offsets1']
        offsets2 = arrays['offsets2']
        bondlists = np.concatenate((i, j, offsets1, offsets2, p),
                                   axis=1)
        # bondlists = bondlists.astype(int)
        # print('get_arrays: %s'%(time() - tstart))
        return bondlists

    @property
    def search_bond(self):
        """search_bond object."""
        if self._search_bond is not None:
            return self._search_bond
        search_bond = SearchBond(self.label, batoms=self.batoms)
        self._search_bond = search_bond
        return search_bond

    @property
    def show_search(self):
        return self.settings.coll.Bbond.show_search

    @show_search.setter
    def show_search(self, show_search):
        self.settings.coll.Bbond.show_search = show_search
        if not show_search:
            self.search_bond.set_arrays(default_search_bond_datas.copy())
            self._search_bond = None
        else:
            self.update()

    @property
    def show_hydrogen_bond(self):
        return self.settings.coll.Bbond.show_hydrogen_bond

    @show_hydrogen_bond.setter
    def show_hydrogen_bond(self, show_hydrogen_bond):
        self.settings.coll.Bbond.show_hydrogen_bond = show_hydrogen_bond
        self.update()

    def __getitem__(self, indices):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        from batoms.bond.slicebonds import SliceBonds
        slicebonds = SliceBonds(self.label, indices, bonds=self)
        return slicebonds

    def __setitem__(self, indices, value):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[indices] = value
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
        positions = np.tile(self.positions, (M,) + (1,) *
                            (len(self.positions.shape) - 1))
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
        bbond = self.__class__(label, species,
                               self.local_positions,
                               location=self.obj.location,
                               scale=self.scale,
                               material=self.material)
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
        s = "Bonds(Total: {:6d}, {}" .format(len(self), self.arrays)
        return s

    def build_bondlists(self, species, positions, cell, pbc, setting):
        """
        build bondlist for atoms
        steps:
        1 build bondlist for atoms with pbc
        2 search connected_components (molecule),
          return pecies data and its neighbour
        3 add bondlist related with molecule

        """
        from batoms.neighborlist import bondlist_kdtree
        bondlists = np.zeros((0, 11), dtype=int)
        bonddatas = {}
        if len(setting) == 0:
            return bondlists, bonddatas, {}, {}
        #
        # tstart = time()
        # ==========================================================
        # step 1 build bondlist for atoms with pbc
        # nli: index1
        # nlj: index2
        # nlk: search bond style
        # nlp: polyhedra
        # nlt: bond type: hydrogen bond
        # nlSj: offset of atoms in nlj
        nli, nlj, nlk, nlp, nlt, nlSj = bondlist_kdtree('ijkptS', species,
                                                   positions, cell,
                                                   pbc, setting)
        nb = len(nli)
        nlSi = np.zeros((nb, 3))
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        # search type 0,
        search0 = np.where((nlk == 0) & (
            nlSj != np.array([0, 0, 0])).any(axis=1), False, True)
        # 0  1  2:5       5:8          8     9
        # i, j, offset_i, offset_j, search, search_style
        bondlists = np.concatenate((np.array([nli, nlj]).T,
                                    np.array(nlSi, dtype=int), np.array(nlSj),
                                    np.array(nlk).reshape(-1, 1),
                                    np.array(nlp).reshape(-1, 1),
                                    np.array(nlt).reshape(-1, 1)),
                                    axis=1)
        bondlists = bondlists.astype(int)
        # remove bond outside box for search0
        # not now, we need all bonds here, and then add a final check.
        # bondlists = bondlists[search0]
        # build bondatas, for each atom, save the bonds connect to it.
        argsort = bondlists[:, 0].argsort()
        bondlists = bondlists[argsort]
        u, indices = np.unique(bondlists[:, 0], return_index=True)
        indices = np.append(indices, len(bondlists))
        m = len(u)
        for i in range(m):
            bonddatas[u[i]] = bondlists[indices[i]: indices[i + 1]]
        #
        # ===================================================
        # 2 search connected_components (molecule),
        #   return pecies data and its neighbour
        peciesBondLists, molPeciesDatas = self.build_peciesBondLists(
            len(positions), bondlists)
        # build peciesBondDatas
        peciesBondDatas = {}
        for p in molPeciesDatas:
            peciesBondDatas[p] = []
        # since this is search bond type 2, bothway, first i
        argsort = peciesBondLists[:, 0].argsort()
        peciesBondLists = peciesBondLists[argsort]
        u, indices = np.unique(peciesBondLists[:, 0], return_index=True)
        indices = np.append(indices, len(peciesBondLists))
        m = len(u)
        # for each pecies, save its neighbour and offset
        for i in range(m):
            peciesBondDatas[u[i]] = peciesBondLists[indices[i]: indices[i + 1]]
        # bothway, then j
        argsort = peciesBondLists[:, 1].argsort()
        peciesBondLists = peciesBondLists[argsort]
        u, indices = np.unique(peciesBondLists[:, 1], return_index=True)
        indices = np.append(indices, len(peciesBondLists))
        m = len(u)
        # for each pecies, save its neighbour and offset
        for i in range(m):
            data = peciesBondLists[indices[i]: indices[i + 1]]
            data = data.reshape(-1, 11)
            data[:, 5:8] *= -1
            data[:, [0, 1]] = data[:, [1, 0]]
            peciesBondDatas[u[i]] = data
        # ==========================================================
        # 3 add bondlist related with molecule
        # add bondlist for molecule
        for mol in peciesBondLists:
            indices = molPeciesDatas[mol[1]]
            for i in indices:
                if i not in bonddatas:
                    continue
                data = bonddatas[i]
                n = len(data)
                bondlists = np.append(bondlists,  data, axis=0)
                bondlists[-n:, 2:5] += mol[5:8]
                bondlists[-n:, 5:8] += mol[5:8]
            indices = molPeciesDatas[mol[0]]
            for i in indices:
                if i not in bonddatas:
                    continue
                data = bonddatas[i]
                n = len(data)
                bondlists = np.append(bondlists,  data, axis=0)
                bondlists[-n:, 2:5] -= mol[5:8]
                bondlists[-n:, 5:8] -= mol[5:8]
        # print(bondlists)
        bondlists = bondlists.astype(int)
        bondlists = np.unique(bondlists, axis=0)
        self.peciesBondLists = peciesBondLists
        self.molPeciesDatas = molPeciesDatas
        self.peciesBondDatas = peciesBondDatas
        return bondlists, bonddatas, peciesBondDatas, molPeciesDatas

    def build_peciesBondLists(self, natom, bondlists):
        """
        search type 2: build molecule
        steps:
        1 search connected_components (molecules) inside atoms
        2 construct the molecules by its pecies and the coresponding offsets
        3 for
        3 return the molecules
        """
        from scipy.sparse import csgraph, csr_matrix
        molPeciesDatas = {}
        peciesBondLists = np.zeros((0, 11), dtype=int)
        # search type 2
        k = bondlists[:, 8]
        indices = np.where(k == 2)[0]
        ns2 = len(indices)
        if ns2 == 0:
            return peciesBondLists, molPeciesDatas
        bondlists1 = bondlists[indices, :]
        # ========================================================
        # 1 search connected_components (molecules) inside atoms
        # with crossed bond, entire molecule
        molDatas = {}
        ai = bondlists1[:, 0]
        aj = bondlists1[:, 1]
        data = np.ones(ns2, dtype=int)
        matrix = csr_matrix((data, (ai, aj)), shape=(natom, natom))
        n_components1, component_list1 = csgraph.connected_components(matrix)
        for i in range(n_components1):
            indices = np.where(component_list1 == i)[0]
            n = len(indices)
            if n < 2:
                continue
            molDatas[i] = {'sub': []}
            molDatas[i]['indices'] = indices
            # TODO: check here
            molDatas[i]['offsets'] = indices
        # without crossed bond, small pecies of molDatas
        mask = np.where((bondlists1[:, 2:5] != bondlists1[:, 5:8]).any(axis=1),
                        False, True)
        data = data[mask]
        ai = ai[mask]
        aj = aj[mask]
        matrix = csr_matrix((data, (ai, aj)), shape=(natom, natom))
        n_components2, component_list2 = csgraph.connected_components(matrix)
        #
        for i in range(n_components2):
            indices = np.where(component_list2 == i)[0]
            if component_list1[indices[0]] in molDatas:
                # this pecies belong to molDatas
                molDatas[component_list1[indices[0]]]['sub'].append(i)
                molPeciesDatas[i] = indices
        # pprint(molPeciesDatas)
        # cross box bond, find direct neighbour pecies
        bondlists2 = bondlists1[~mask]
        n = len(bondlists2)
        ai = bondlists2[:, 0].astype(int)
        aj = bondlists2[:, 1].astype(int)
        peciesBondLists = np.zeros((n, 11), dtype=int)
        molBondDicts = {}
        for i in range(n):
            # molDatasId = component_list1[bondlists2[i, 0]]
            i1 = component_list2[ai[i]]
            j1 = component_list2[aj[i]]
            peciesBondLists[i, 0] = i1
            peciesBondLists[i, 1] = j1
            peciesBondLists[i, 2:] = bondlists2[i, 2:]
            molBondDicts[(i1, j1)] = bondlists2[i, 2:]
        peciesBondLists = np.unique(peciesBondLists, axis=0)
        # ========================================================
        # find indirect neighbour pecies inside mol1
        # find connected path between pecies
        ai = peciesBondLists[:, 0].astype(int)
        aj = peciesBondLists[:, 1].astype(int)
        nml = len(ai)
        data = np.ones(nml, dtype=int)
        graph = csr_matrix((data, (ai, aj)), shape=(
            n_components2, n_components2))
        dist_matrix, predecessors = csgraph.shortest_path(
            csgraph=graph, return_predecessors=True)
        # print('dist_matrix: ', dist_matrix)
        # print('predecessors: ', predecessors)
        dist_matrix = dist_matrix.astype(int)
        mollist = np.zeros(11, dtype=int)
        for i, data in molDatas.items():
            indices = data['sub']
            n = len(indices)
            if n < 2:
                continue
            for j in range(n - 1):
                for k in range(j + 1, n):
                    if dist_matrix[indices[j], indices[k]] == 1:
                        continue
                    path = [indices[k]]
                    end = indices[k]
                    last = indices[k]
                    for i1 in range(dist_matrix[indices[j], end]):
                        last = predecessors[indices[j], last]
                        path = [last] + path
                    offsets = np.zeros(3)
                    for i2 in range(0, len(path) - 1):
                        offsets += molBondDicts[(path[i2], path[i2 + 1])][3:6]
                    mollist[0] = path[0]
                    mollist[1] = path[-1]
                    mollist[5:8] = offsets
                    peciesBondLists = np.append(
                        peciesBondLists, np.array([mollist]), axis=0)
        self.molDatas = molDatas
        return peciesBondLists, molPeciesDatas

    def build_bondlists_with_boundary(self, arrays, bondlists, bonddatas,
                                      peciesBondDatas, molPeciesDatas):
        """
        build extra bondlists based on boundary atoms
        """
        n = len(arrays['positions'])
        if n == 0:
            return bondlists
        # search bond type 0 and 1
        for i in range(n):
            if arrays['indices'][i] not in bonddatas:
                # print('no bonds: ', arrays['indices'][i])
                continue
            data = bonddatas[arrays['indices'][i]]
            n = len(data)
            bondlists = np.append(bondlists,  data, axis=0)
            bondlists[-n:, 2:5] += arrays['offsets'][i]
            bondlists[-n:, 5:8] += arrays['offsets'][i]
            # todo: in this case, some of the atoms overlap
            # with original atoms.
            # since it doesn't influence the 3d view, we
            # just leave it like this
            # bondlists[-n:, 8] = 1
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        bondlists = np.unique(bondlists, axis=0)
        # search bond type 2
        # divide boudanry atoms by offsets
        offsets = np.unique(arrays['offsets'], axis=0)
        n = len(offsets)*len(molPeciesDatas)
        peciesArrays = np.zeros((n, 4), dtype=int)
        n = 0
        for offset in offsets:
            indices = np.where((arrays['offsets'] == offset).all(axis=1))[0]
            indices = arrays['indices'][indices]
            for p, indices0 in molPeciesDatas.items():
                indices1 = np.intersect1d(indices0, indices)
                if len(indices1) > 0:
                    peciesArrays[n][0] = p
                    peciesArrays[n][1:4] = offset
                    n += 1
        peciesArrays = peciesArrays[:n]
        npa = len(peciesArrays)
        if npa == 0:
            return bondlists
        # search bond type 2
        for i in range(npa):
            peciesId = peciesArrays[i, 0]
            p = peciesBondDatas[peciesId]
            # repeat itself
            indices = molPeciesDatas[peciesId]
            for j in indices:
                if j not in bonddatas:
                    continue
                data = bonddatas[j]
                n = len(data)
                bondlists = np.append(bondlists,  data, axis=0)
                bondlists[-n:, 2:5] += peciesArrays[i, 1:4]
                bondlists[-n:, 5:8] += peciesArrays[i, 1:4]
            # repeat its neighbour pecies
            for pb in p:
                peciesId = pb[1]
                indices = molPeciesDatas[peciesId]
                for j in indices:
                    if j not in bonddatas:
                        continue
                    data = bonddatas[j]
                    n = len(data)
                    bondlists = np.append(bondlists,  data, axis=0)
                    bondlists[-n:, 2:5] += peciesArrays[i, 1:4] + pb[5:8]
                    bondlists[-n:, 5:8] += peciesArrays[i, 1:4] + pb[5:8]
            # todo: in this case, some of the atoms overlap with
            # original atoms.
            # since it doesn't influence the 3d view,
            # we just leave it like this
            # bondlists[-n:, 8] = 1
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        bondlists = np.unique(bondlists, axis=0)
        # search bond type 2

        return bondlists

    def check_boundary(self, bondlists, eps = 1e-6):
        """check boundary for bond search 0

        Args:
            eps (float): default 1e-6
        """
        from ase.geometry import complete_cell
        arrays = self.batoms.arrays
        positions = arrays['positions']
        cell = self.batoms.cell
        # get scaled positions
        positions = np.linalg.solve(complete_cell(cell).T,
                                positions.T).T
        npositions = positions[bondlists[:, 1]] + bondlists[:, 5:8]
        boundary = self.batoms.boundary.boundary
        # boundary condition
        mask1 = np.where((npositions[:, 0] > boundary[0][0] - eps) &
                        (npositions[:, 0] < boundary[0][1] + eps) &
                        (npositions[:, 1] > boundary[1][0] - eps) &
                        (npositions[:, 1] < boundary[1][1] + eps) &
                        (npositions[:, 2] > boundary[2][0] - eps) &
                        (npositions[:, 2] < boundary[2][1] + eps), False, True)
        mask2 = np.where(bondlists[:, 8] == 0, True, False)
        mask = ~(mask1&mask2)
        bondlists = bondlists[mask]
        return bondlists

    def search_molecule(self, natom, bondlists):
        """
        """
        from scipy.sparse import csgraph, csr_matrix
        mols = {}
        # search type
        k = bondlists[:, -1]
        indices = np.where(k == 2)[0]
        ns2 = len(indices)
        if ns2 == 0:
            return mols
        i = bondlists[indices, 0]
        j = bondlists[indices, 1]
        data = np.ones(ns2, dtype=int)
        matrix = csr_matrix((data, (i, j)), shape=(natom, natom))
        # print(matrix)
        n_components, component_list = csgraph.connected_components(matrix)
        # print(n_components)
        # print(component_list)
        for i in range(n_components):
            indices = np.where(component_list == i)[0]
            if len(indices) < 2:
                continue
            mols[i] = indices
        return mols

    def calc_bond_data(self, speciesarray, positions, cell,
                       bondlists, bondsetting,
                       model_styles):
        """
        positions: 3-demision
        todo: support frame for positions and offsets.
        """
        tstart = time()
        if not self.show_hydrogen_bond:
            bondlists = bondlists[bondlists[:, 10] != 1]
        # properties
        atoms_index1 = np.array(bondlists[:, 0], dtype=int)
        atoms_index2 = np.array(bondlists[:, 1], dtype=int)
        # second bond for high order bond
        nb = len(atoms_index1)
        second_bond = np.roll(np.arange(nb), 1)
        # atoms_index3 and atoms_index4 will be removed for Blender 3.1
        atoms_index3 = np.roll(atoms_index1, 1)
        atoms_index4 = np.roll(atoms_index2, 1)
        shows = np.ones(nb, dtype=int)
        orders = np.zeros(nb, dtype=int)  # 1, 2, 3
        styles = np.zeros(nb, dtype=int)  # 0, 1, 2
        widths = np.ones(nb, dtype=float)
        species_index1 = np.zeros(nb, dtype=int)
        species_index2 = np.zeros(nb, dtype=int)
        model_styles = model_styles[atoms_index1]
        if nb == 0:
            return default_bond_datas
        # ------------------------------------
        # offsets
        offsets1 = np.dot(bondlists[:, 2:5], cell)
        offsets2 = np.dot(bondlists[:, 5:8], cell)
        # polyhedra
        polyhedras = np.array(bondlists[:, 9], dtype=int)
        # bond center
        if len(positions.shape) == 2:
            positions = np.array([positions])
        # here we use frames directly, thus will set centers
        # to shape keys frames
        centers = (positions[:, atoms_index1] + offsets1 +
                   positions[:, atoms_index2] + offsets2)/2.0
        # ---------------------------------------------
        for b in bondsetting:
            spi = b.species1
            spj = b.species2
            indi = (speciesarray[atoms_index1] == spi)
            indj = (speciesarray[atoms_index2] == spj)
            ind = indi & indj
            orders[ind] = b.order
            styles[ind] = int(b.style)
            widths[ind] = b.width
            polyhedras[ind] = b.polyhedra
            species_index1[ind] = string2Number(spi)
            species_index2[ind] = string2Number(spj)
            if b.order > 1:
                second_bond, atoms_index3, atoms_index4 = self.secondBond(
                    b, speciesarray, second_bond,
                    atoms_index3, atoms_index4, bondlists)
        # offsets for second bond
        # this will be remove for Blender 3.1
        offsets3 = np.dot(bondlists[:, 2:5][second_bond], cell)
        offsets4 = np.dot(bondlists[:, 5:8][second_bond], cell)
        datas = {
            'atoms_index1': atoms_index1,
            'atoms_index2': atoms_index2,
            'atoms_index3': atoms_index3,
            'atoms_index4': atoms_index4,
            'second_bond': second_bond,
            'species_index1': species_index1,
            'species_index2': species_index2,
            'centers': centers,
            'offsets1': offsets1,
            'offsets2': offsets2,
            'offsets3': offsets3,
            'offsets4': offsets4,
            # 'widths': widths,
            'show': shows,
            'order': orders,
            'style': styles,
            'model_style': model_styles,
            'polyhedra': polyhedras,
        }
        # print('datas: ', datas)
        logger.debug('calc_bond_data: {0:10.2f} s'.format(time() - tstart))
        return datas

    def secondBond(self, b, speciesarray, second_bond,
                   atoms_index3, atoms_index4, bondlists):
        """
        determine the plane of high order bond.
        v1 = atoms2 - atoms1
        v2 = atoms4 - atoms3
        normal = np.cross(v1, v2)
        """
        spi = b.species1
        spj = b.species2
        # find all bonds belong to this pair spi--spj
        # here use &
        indi = (speciesarray[bondlists[:, 0]] == spi)
        indj = (speciesarray[bondlists[:, 1]] == spj)
        indices1 = np.where(indi & indj)[0]
        # find all bonds connect to this bond pair
        # here use |
        indi1 = (speciesarray[bondlists[:, 1]] == spi)
        indj2 = (speciesarray[bondlists[:, 0]] == spj)
        indices2 = np.where(indi | indj | indi1 | indj2)[0]
        bondlists2 = bondlists[indices2]
        for i in indices1:
            # find another bond for this bond: bondlists[i]
            indices3 = np.where((bondlists2[:, 0] == bondlists[i, 0]) |
                                (bondlists2[:, 1] == bondlists[i, 0]) |
                                (bondlists2[:, 0] == bondlists[i, 1]) |
                                (bondlists2[:, 1] == bondlists[i, 1]))[0]
            indices3 = indices2[indices3]
            if len(indices3) > 1:
                indices4 = np.where(indices3 != i)[0]
                second_bond[i] = indices3[indices4[0]]
                atoms_index3[i] = bondlists[second_bond[i], 0]
                atoms_index4[i] = bondlists[second_bond[i], 1]
        return second_bond, atoms_index3, atoms_index4

    def high_order_bond_plane(self, b, speciesarray,
                              positions, nvec, offsets, bondlists):
        """
        determine the plane of high order bond
        """
        spi = b.species1
        spj = b.species2
        indi = (speciesarray[bondlists[:, 0]] == spi)
        indj = (speciesarray[bondlists[:, 1]] == spj)
        bondlists1 = bondlists[indi & indj]
        nbond = len(bondlists1)
        indi1 = (speciesarray[bondlists[:, 1]] == spi)
        indj2 = (speciesarray[bondlists[:, 0]] == spj)
        bondlists2 = bondlists[indi | indj | indi1 | indj2]
        for i in range(nbond):
            # find another bond
            mask = np.logical_not(((bondlists2[:, 0] == bondlists1[i, 0])
                                   & (bondlists2[:, 1] == bondlists1[i, 1]))
                                  | ((bondlists2[:, 0] != bondlists1[i, 0])
                                     & (bondlists2[:, 0] != bondlists1[i, 1])
                                     & (bondlists2[:, 1] != bondlists1[i, 0])
                                     & (bondlists2[:, 1] != bondlists1[i, 1])))
            localbondlist = bondlists2[mask]
            if len(localbondlist) == 0:
                second_bond = np.array([0.0, 0.0, 1])
            else:
                second_bond = positions[localbondlist[0, 0]
                                        ] - positions[localbondlist[0, 1]]
            norml = np.cross(second_bond, nvec[i]) + np.array([1e-8, 0, 0])
            offset = np.cross(norml, nvec[i]) + np.array([1e-8, 0, 0])
            offsets[i] = offset/np.linalg.norm(offset)
            # print(offsets[i], offset)

        return offsets

    def bond_order_auto_set(self):
        """Set bond order by pybel
        #TODO support show, select
        """
        from openbabel import pybel
        species = self.batoms.arrays['species']

        mol = self.batoms.as_pybel(export_bond=False)
        bondlists = [[b.GetBeginAtom().GetIndex(),
                        b.GetEndAtom().GetIndex(),
                        0, 0, 0, 0, 0, 0, 0, 0, 0]
             for b in pybel.ob.OBMolBondIter(mol.OBMol)]
        #TODO detect hydrogen bond
        cutoff_dict = self.settings.cutoff_dict
        nbond = len(bondlists)
        remove = []
        bondlists = np.array(bondlists)
        for i in range(nbond):
            b = bondlists[i]
            if (species[b[0]], species[b[1]]) not in cutoff_dict:
                # swtich bond
                if (species[b[1]], species[b[0]]) in cutoff_dict:
                    bondlists[i][0], bondlists[i][1] = bondlists[i][1], bondlists[i][0]
                else:
                    # remove bond
                    remove.append(i)
        np.delete(bondlists, remove)
        orders = [b.GetBondOrder()
            for b in pybel.ob.OBMolBondIter(mol.OBMol)]
        self.update(bondlists, orders = orders)

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

    def as_dict(self):
        """
        """
        data = {
            'array': None
        }
        if len(self) > 1:
            data['array'] = dict(self.arrays)
        data['show_search'] = self.show_search
        data['show_hydrogen_bond'] = self.show_hydrogen_bond
        data['settings'] = self.settings.as_dict()
        return data
