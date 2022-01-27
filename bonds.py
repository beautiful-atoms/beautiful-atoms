"""Definition of the Bonds class.

This module defines the Bonds object in the Batoms package.

"""

from pandas import offsets
import bpy
import bmesh
from time import time
from batoms.butils import object_mode
from batoms.tools import string2Number
import numpy as np
from batoms.base import BaseObject
from batoms.bondsetting import BondSettings

default_attributes = [
            ['atoms_index1', 'INT'],
            ['atoms_index2', 'INT'],
            ['atoms_index3', 'INT'],
            ['atoms_index4', 'INT'],
            ['species_index1', 'INT'],
            ['species_index2', 'INT'],
            ['order', 'INT'],
            ['style', 'INT'],
            ['show', 'BOOLEAN'],
            ['model_style', 'INT'],
        ]

default_bond_datas = {
        'atoms_index1': np.ones(0, dtype = int),
        'atoms_index2': np.ones(0, dtype = int),
        'atoms_index3': np.ones(0, dtype = int),
        'atoms_index4': np.ones(0, dtype = int),
        'species_index1': np.ones(0, dtype = int),
        'species_index2': np.ones(0, dtype = int),
        'centers':np.zeros((0, 3)),
        # 'vectors':np.zeros((0, 3)),
        'offsets':np.zeros((0, 3)),
        # 'eulers':np.eye(3),
        # 'lengths':np.zeros((0, 3)),
        'widths': np.ones(0, dtype = float),
        'orders': np.zeros(0, dtype = int),
        'styles': np.zeros(0, dtype = int),
        'model_styles':np.ones(0, dtype = int),
        }

class Bonds(BaseObject):
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
                bond_datas = None,
                location = np.array([0, 0, 0]),
                batoms = None,
                 ):
        #
        self.batoms = batoms
        self.label = label
        obj_name = '%s_bond_center'%(self.label)
        bobj_name = 'bbond'
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        self.setting = BondSettings(self.label, batoms = batoms, bonds = self)
        if bond_datas is not None:
            self.build_object(bond_datas)
        else:
            self.load(label)
    
    def build_object(self, bond_datas, attributes = {}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        from scipy.spatial.transform import Rotation as R
        tstart = time()
        if len(bond_datas['centers'].shape) == 2:
            self._frames = (np.array([bond_datas['centers']]), 
                            # np.array([bond_datas['vectors']]),
                            np.array([bond_datas['offsets']]),
                            # np.array([bond_datas['eulers']]),
                            # np.array([bond_datas['lengths']]),
                            )
            centers = bond_datas['centers']
            # vectors = bond_datas['vectors']
            offsets = bond_datas['offsets']
            # eulers = bond_datas['eulers']
            # lengths = bond_datas['lengths']
        elif len(bond_datas['centers'].shape) == 3:
            self._frames = (bond_datas['centers'], 
                            # bond_datas['vectors'], 
                            bond_datas['offsets'], 
                            # bond_datas['eulers'], 
                            # bond_datas['lengths'],
                            )
            centers = bond_datas['centers'][0]
            # vectors = bond_datas['vectors'][0]
            offsets = bond_datas['offsets'][0]
            # eulers = bond_datas['eulers'][0]
            # lengths = bond_datas['lengths'][0]
        else:
            raise Exception('Shape of centers is wrong!')
        nbond = len(bond_datas['centers'])
        show = np.ones(nbond, dtype = int)
        attributes.update({
                            'atoms_index1': bond_datas['atoms_index1'], 
                            'atoms_index2': bond_datas['atoms_index2'], 
                            'atoms_index3': bond_datas['atoms_index3'], 
                            'atoms_index4': bond_datas['atoms_index4'], 
                            'species_index1': bond_datas['species_index1'], 
                            'species_index2': bond_datas['species_index2'], 
                            'show': show,
                            'model_style': bond_datas['model_styles'],
                            'style': bond_datas['styles'],
                            'order': bond_datas['orders'],
                            })
        name = '%s_bond_center'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(centers, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = 'POINT')
        self.setting.coll.objects.link(obj)
        #
        # name = '%s_bond_vector'%self.label
        # if name in bpy.data.objects:
        #     obj = bpy.data.objects.get(name)
        #     bpy.data.objects.remove(obj, do_unlink = True)
        # mesh = bpy.data.meshes.new(name)
        # mesh.from_pydata(vectors, [], [])
        # mesh.update()
        # obj = bpy.data.objects.new(name, mesh)
        # self.setting.coll.objects.link(obj)
        # obj.hide_set(True)
        #
        name = '%s_bond_offset'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.setting.coll.objects.link(obj)
        obj.hide_set(True)
        # # calc euler angles
        # name = '%s_bond_rotation'%self.label
        # if name in bpy.data.objects:
        #     obj = bpy.data.objects.get(name)
        #     bpy.data.objects.remove(obj, do_unlink = True)
        # mesh = bpy.data.meshes.new(name)
        # mesh.from_pydata(eulers, [], [])
        # mesh.update()
        # obj = bpy.data.objects.new(name, mesh)
        # self.setting.coll.objects.link(obj)
        # obj.hide_set(True)
        #
        # scales = np.concatenate((np.ones((nbond, 1)), 
        #                 np.ones((nbond, 1)),
        #                 lengths.reshape(-1, 1)), axis = 1)
        # name = '%s_bond_scale'%self.label
        # if name in bpy.data.objects:
        #     obj = bpy.data.objects.get(name)
        #     bpy.data.objects.remove(obj, do_unlink = True)
        # mesh = bpy.data.meshes.new(name)
        # mesh.from_pydata(scales, [], [])
        # mesh.update()
        # obj = bpy.data.objects.new(name, mesh)
        # self.setting.coll.objects.link(obj)
        # obj.hide_set(True)
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis = True)
        print('bonds: build_object: {0:10.2f} s'.format(time() - tstart))
    
    def load(self, label):
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
        #------------------------------------------------------------------
        # select attributes
        GroupInput = modifier.node_group.nodes.get('Group Input')
        # add new output sockets
        for att in default_attributes:
            GroupInput.outputs.new(type = att[1], name = att[0])
        # the above codes not works. maybe bug in blender, 
        # we add this, maybe deleted in the future
        for i in range(1, 11):
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
        TransferBatoms = []
        for i in range(4):
            tmp = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_TransferBatoms%s'%(self.label, i),
                        'GeometryNodeAttributeTransfer')
            tmp.mapping = 'INDEX'
            tmp.data_type = 'FLOAT_VECTOR'
            TransferBatoms.append(tmp)
        for i in range(4):
            gn.node_group.links.new(ObjectBatoms.outputs['Geometry'], TransferBatoms[i].inputs['Target'])
            gn.node_group.links.new(PositionBatoms.outputs['Position'], TransferBatoms[i].inputs['Attribute'])
            gn.node_group.links.new(GroupInput.outputs[i + 1], TransferBatoms[i].inputs['Index'])
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
        # we need three add operations
        # two: Get the positions with offset for atoms2, and atom4
        # one: Get center = (positions1 + positions2)/2
        VectorAdd = []
        for i in range(3):
            tmp = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_VectorAdd%s'%(self.label, i),
                        'ShaderNodeVectorMath')
            tmp.operation = 'ADD'
            VectorAdd.append(tmp)
        gn.node_group.links.new(TransferBatoms[1].outputs[0], VectorAdd[0].inputs[0])
        gn.node_group.links.new(TransferOffsets.outputs[0], VectorAdd[0].inputs[1])
        gn.node_group.links.new(TransferBatoms[3].outputs[0], VectorAdd[1].inputs[0])
        gn.node_group.links.new(TransferOffsets.outputs[0], VectorAdd[1].inputs[1])
        #
        # divide by 2 to get the center
        VectorDivide = get_nodes_by_name(gn.node_group.nodes, 
                    'VectorDivide_%s'%self.label,
                    'ShaderNodeVectorMath')
        VectorDivide.operation = 'DIVIDE'
        VectorDivide.inputs[1].default_value = (2, 2, 2)
        gn.node_group.links.new(TransferBatoms[0].outputs[0], VectorAdd[2].inputs[0])
        gn.node_group.links.new(VectorAdd[0].outputs[0], VectorAdd[2].inputs[1])
        gn.node_group.links.new(VectorAdd[2].outputs[0], VectorDivide.inputs[0])
        # set center of the bond
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label, 
                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(VectorDivide.outputs[0], SetPosition.inputs['Position'])
        # get the vector for the bond and the length
        # also the vector for the second bond
        VectorSubtract = []
        for i in range(2):
            tmp = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_VectorSubtract%s'%(self.label, i),
                        'ShaderNodeVectorMath')
            tmp.operation = 'SUBTRACT'
            VectorSubtract.append(tmp)
        VectorLength = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorLength'%self.label,
                    'ShaderNodeVectorMath')
        VectorLength.operation = 'LENGTH'
        VectorCross0 = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorCross0'%self.label,
                    'ShaderNodeVectorMath')
        VectorCross0.operation = 'CROSS_PRODUCT'
        gn.node_group.links.new(TransferBatoms[0].outputs[0], VectorSubtract[0].inputs[0])
        gn.node_group.links.new(VectorAdd[0].outputs[0], VectorSubtract[0].inputs[1])
        gn.node_group.links.new(TransferBatoms[2].outputs[0], VectorSubtract[1].inputs[0])
        gn.node_group.links.new(VectorAdd[1].outputs[0], VectorSubtract[1].inputs[1])
        # calc the bond length, use it to build scale
        gn.node_group.links.new(VectorSubtract[0].outputs[0], VectorLength.inputs[0])
        #
        CombineXYZ = get_nodes_by_name(gn.node_group.nodes,
                        '%s_CombineXYZ'%self.label, 
                        'ShaderNodeCombineXYZ')
        CombineXYZ.inputs[0].default_value = 1
        CombineXYZ.inputs[1].default_value = 1
        gn.node_group.links.new(VectorLength.outputs['Value'], CombineXYZ.inputs['Z'])
        # cross for rotation, for high order bond
        gn.node_group.links.new(VectorSubtract[0].outputs[0], VectorCross0.inputs[0])
        gn.node_group.links.new(VectorSubtract[1].outputs[0], VectorCross0.inputs[1])
        # get Euler for rotation
        # we need align two vectors to fix a plane
        # we build the instancer by fix bond diection to Z, and high order bond shift to X
        # thus the the normal of high order bond plane is Y
        AlignEuler = []
        for i in range(2):
            tmp = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_AlignEuler%s'%(self.label, i),
                        'FunctionNodeAlignEulerToVector')
            AlignEuler.append(tmp)
        AlignEuler[0].axis = 'Z'
        AlignEuler[1].axis = 'Y'
        # We should fix Z when align Y
        AlignEuler[1].pivot_axis = 'Z'
        gn.node_group.links.new(VectorSubtract[0].outputs[0], AlignEuler[0].inputs['Vector'])
        gn.node_group.links.new(AlignEuler[0].outputs[0], AlignEuler[1].inputs['Rotation'])
        gn.node_group.links.new(VectorCross0.outputs[0], AlignEuler[1].inputs['Vector'])
        # find bond kinds by the names of species
        for sp in self.batoms.species:
            # we need two compares for one species,
            # because we have two sockets: species_index1 and species_index2
            CompareSpecies = []
            for i in range(2):
                tmp = get_nodes_by_name(gn.node_group.nodes, 
                            '%s_CompareSpecies_%s_%s'%(self.label, sp.name, i),
                            'FunctionNodeCompareFloats')
                tmp.operation = 'EQUAL'
                tmp.inputs[1].default_value = string2Number(sp.name)
                CompareSpecies.append(tmp)
            gn.node_group.links.new(GroupInput.outputs[5], CompareSpecies[0].inputs[0])
            gn.node_group.links.new(GroupInput.outputs[6], CompareSpecies[1].inputs[0])
        # order 
        for order in [1, 2, 3]:
            CompareOrder = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareFloats_%s_%s_order'%(self.label, order),
                        'FunctionNodeCompareFloats')
            CompareOrder.operation = 'EQUAL'
            CompareOrder.inputs[1].default_value = order
            gn.node_group.links.new(GroupInput.outputs[7], CompareOrder.inputs[0])
        # style 
        for style in [0, 1, 2]:
            CompareStyle = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareFloats_%s_%s_style'%(self.label, style),
                        'FunctionNodeCompareFloats')
            CompareStyle.operation = 'EQUAL'
            CompareStyle.inputs[1].default_value = style
            gn.node_group.links.new(GroupInput.outputs[8], CompareStyle.inputs[0])
        #
        for sp in self.setting:
            self.add_geometry_node(sp.as_dict())
        
        print('Build geometry nodes for bonds: %s'%(time() - tstart))

    def add_geometry_node(self, sp):
        from batoms.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label)
        #
        order = sp['order']
        style = int(sp['style'])
        name = '%s_%s_%s_%s'%(self.label, sp["name"], order, style)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                    'InstanceOnPoint_%s'%name,
                    'GeometryNodeInstanceOnPoints')
        ObjectInstancer = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s'%name,
                    'GeometryNodeObjectInfo')
        ObjectInstancer.inputs['Object'].default_value = \
                        self.setting.instancers[sp["name"]]['%s_%s'%(order, style)]
        #
        BoolSpecies = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_BooleanMath_species'%name,
                        'FunctionNodeBooleanMath')
        BoolOrder = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_BooleanMath_order'%name,
                        'FunctionNodeBooleanMath')
        BoolStyle = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_BooleanMath_style'%name,
                        'FunctionNodeBooleanMath')
        BoolModelStyle = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_BooleanMath_modelstyle'%name,
                        'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_BooleanMath_show'%name,
                    'FunctionNodeBooleanMath')
        BoolBondLength = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_BoolBondLength'%name,
                    'FunctionNodeBooleanMath')
        #bondlength larger than max will not show
        LessBondLength = get_nodes_by_name(gn.node_group.nodes,
                        '%s_LessBondLength'%name, 
                        'ShaderNodeMath')
        LessBondLength.operation = 'LESS_THAN'
        LessBondLength.inputs[1].default_value = sp['max']
        VectorLength = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorLength'%self.label)
        #
        CompareSpecies0 = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_CompareSpecies_%s_0'%(self.label, sp["species1"]))
        CompareSpecies1 = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_CompareSpecies_%s_1'%(self.label, sp["species2"]))
        CompareOrder = get_nodes_by_name(gn.node_group.nodes, 
                'CompareFloats_%s_%s_order'%(self.label, order))
        CompareStyle = get_nodes_by_name(gn.node_group.nodes, 
                'CompareFloats_%s_%s_style'%(self.label, style))
        #
        CombineXYZ = get_nodes_by_name(gn.node_group.nodes,
                        '%s_CombineXYZ'%self.label, 
                        'ShaderNodeCombineXYZ')
        AlignEuler1 = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_AlignEuler1'%self.label,
                    'FunctionNodeAlignEulerToVector')
        gn.node_group.links.new(SetPosition.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[9], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[10], BoolModelStyle.inputs[0])
        gn.node_group.links.new(CompareSpecies0.outputs[0], BoolSpecies.inputs[0])
        gn.node_group.links.new(CompareSpecies1.outputs[0], BoolSpecies.inputs[1])
        gn.node_group.links.new(BoolSpecies.outputs[0], BoolOrder.inputs[0])
        gn.node_group.links.new(CompareOrder.outputs[0], BoolOrder.inputs[1])
        gn.node_group.links.new(BoolOrder.outputs[0], BoolStyle.inputs[0])
        gn.node_group.links.new(CompareStyle.outputs[0], BoolStyle.inputs[1])
        gn.node_group.links.new(BoolStyle.outputs[0], BoolModelStyle.inputs[1])
        gn.node_group.links.new(BoolModelStyle.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(VectorLength.outputs['Value'], LessBondLength.inputs[0])
        gn.node_group.links.new(BoolShow.outputs['Boolean'], BoolBondLength.inputs[0])
        gn.node_group.links.new(LessBondLength.outputs[0], BoolBondLength.inputs[1])
        gn.node_group.links.new(BoolBondLength.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(CombineXYZ.outputs[0], InstanceOnPoint.inputs['Scale'])
        #
        gn.node_group.links.new(AlignEuler1.outputs[0], InstanceOnPoint.inputs['Rotation'])
        #
        gn.node_group.links.new(ObjectInstancer.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])

    def update_bonds(self, ):
        """Draw bonds.
        calculate bond in all farmes, and save
        get the max number of bonds for each pair
        draw the bonds
        add shape key
        the extra bonds use the find bond data.
        """
        
        from batoms.butils import clean_coll_objects
        # if not self.bondlist:
        object_mode()
        # clean_coll_objects(self.coll, 'bond')
        frames = self.batoms.get_frames()
        arrays = self.batoms.arrays
        size = arrays['radius']*arrays['scale']
        species = arrays['species']
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        bond_datas = {}
        tstart = time()
        for f in range(nframe):
            print('update bond: ', f)
            positions = frames[f]
            # if len(frames_boundary) > 0:
            #     positions_boundary = frames_boundary[f]
            #     positions = positions + positions_boundary
            # if len(frames_search) > 0:
            #     positions_search = frames_search[f]
            #     positions = positions + positions_search
            bondlist = build_bondlists(species, positions, 
                        self.batoms.cell, self.batoms.pbc, self.setting.cutoff_dict)
            bond_kinds = calc_bond_data(species, positions, 
                        self.batoms.cell, size, bondlist, self.setting,
                        arrays['model_style'])
            if f == 0:
                bond_datas = bond_kinds
            else:
                for kind, bond_data in bond_kinds.items():
                    if kind not in bond_datas:
                        bond_datas[kind] = bond_kinds[kind]
                    else:
                        bond_datas[kind]['positions'].extend(bond_data['positions'])
                        if len(bond_data['positions']) > 0:
                            if bond_datas[kind]['nposition'] < len(bond_data['positions'][0]):
                                bond_datas[kind]['nposition'] = len(bond_data['positions'][0])
        # print('calc bond: {0:10.2f} s'.format(time() - tstart))
            # nframe = len(bond_datas['positions'])
            # if nframe == 0: continue
            # change to same length
            # find max
            # nbs = [bond_data['positions'][i].shape[0] for i in range(nframe)]
            # nb_max = max(nbs)
            # frames_bond = np.zeros((nframe, nb_max, 7))
            # for i in range(nframe):
            #     frames_bond[i, 0:nbs[i], :] = bond_data['positions'][i]
            #     frames_bond[i, nbs[i]:, 0:3] = frames[i][0]
        if len(bond_datas) == 0:
            return
        self.set_bond_datas(bond_datas)
        # self.coll.objects.link(bb.obj)
        # bpy.data.collections['Collection'].objects.unlink(bb.obj)
        # bb.set_frames()
        # bpy.context.scene.frame_set(self.batoms.nframe)
        print('draw bond: {0:10.2f} s'.format(time() - tstart))

    @property
    def obj_v(self):
        return self.get_obj_v()
    
    def get_obj_v(self):
        name = '%s_bond_vector'%self.label
        obj_v = bpy.data.objects.get(name)
        if obj_v is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_v

    @property
    def obj_o(self):
        return self.get_obj_o()
    
    def get_obj_o(self):
        name = '%s_bond_offset'%self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_o

    @property
    def obj_r(self):
        return self.get_obj_r()
    
    def get_obj_r(self):
        name = '%s_bond_rotation'%self.label
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

    def get_arrays(self):
        """
        """
        object_mode()
        tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions[0],
                        # 'vectors': self.positions[1],
                        'offsets': self.positions[2],
                        # 'scales': self.positions[3],
                        })
        # print('get_arrays: %s'%(time() - tstart))
        return arrays
    
    @property
    def bondlists(self):
        return self.get_bondlists()
    
    def get_bondlists(self):
        """
        """
        object_mode()
        bondlists = self.arrays
        return bondlists

    @property
    def bond_datas(self):
        return self.get_bond_datas()
    
    @bond_datas.setter    
    def bond_datas(self, bond_datas):
        self.set_bond_datas(bond_datas)
    
    def get_bond_datas(self):
        return np.array(self.mesh.bond_datas)
    
    def set_bond_datas(self, bond_datas):
        """
        """
        if len(bond_datas['centers']) == 0:
            return
        attributes = self.attributes
        # same length
        if len(bond_datas['centers']) == len(attributes['show']):
            self.set_positions(bond_datas['centers'],
                                # bond_datas['vectors'],
                                bond_datas['offsets'],
                                # bond_datas['eulers'],
                                # bond_datas['lengths'],
                                )
            self.set_attributes({'atoms_index1': bond_datas['atoms_index1']})
            self.set_attributes({'atoms_index2': bond_datas['atoms_index2']})
            self.set_attributes({'atoms_index3': bond_datas['atoms_index3']})
            self.set_attributes({'atoms_index4': bond_datas['atoms_index4']})
            self.set_attributes({'species_index1': bond_datas['species_index1']})
            self.set_attributes({'species_index2': bond_datas['species_index2']})
            self.set_attributes({'order': bond_datas['orders']})
            self.set_attributes({'style': bond_datas['styles']})
        else:
            # add or remove vertices
            self.build_object(bond_datas)

    @property
    def local_positions(self):
        return self.get_local_positions()
    
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = []
        for obj in [self.obj, self.obj_o]:
            co = np.empty(n*3, dtype=np.float64)
            obj.data.vertices.foreach_get('co', co)  
            local_positions.append(co.reshape((n, 3)))
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
        positions = []
        objs = [self.obj, self.obj_o]
        for i in range(2):
            obj = objs[i]
            positions.append(local2global(self.local_positions[i], 
                np.array(obj.matrix_world)))
        return positions
    
    def set_positions(self, centers, offsets):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        n = len(self.obj.data.vertices)
        # scales = np.concatenate((np.ones((n, 1)), 
                        # np.ones((n, 1)),
                        # lengths.reshape(-1, 1)), axis = 1)
        for obj, positions in zip([self.obj, self.obj_o], 
                    [centers, offsets]):
            if len(positions) != n:
                raise ValueError('positions has wrong shape %s != %s.' %
                                    (len(positions), n))
            positions = local2global(positions, 
                    np.array(obj.matrix_world), reversed = True)
            positions = positions.reshape((n*3, 1))
            obj.data.shape_keys.key_blocks[0].data.foreach_set('co', positions)
            obj.data.update()
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode = 'EDIT')
            bpy.ops.object.mode_set(mode = 'OBJECT')
    
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
        # print('set_attributes: %s'%(time() - tstart))

    def set_attribute_with_indices(self, name, indices, data):
        data0 = self.attributes[name]
        data0[indices] = data
        self.set_attributes({name: data0})
    
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
        arrays.update({'positions': self.positions[0],
                        })
        print('get_arrays: %s'%(time() - tstart))
        return arrays
    
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
        centers = frames
        nframe = len(centers)
        if nframe == 0 : return
        nbond = len(centers[0])
        # scales = np.concatenate((np.ones((nframe, nbond, 1)), 
                        # np.ones((nframe, nbond, 1)),
                        # lengths.reshape(nframe, nbond, 1)), axis = 2)
        for sp, frames in zip(['center'],
                                [centers]):
            name = '%s_bond_%s'%(self.label, sp)
            obj = bpy.data.objects.get(name)
            base_name = 'Basis_%s_bond_%s'%(self.label, sp)
            if obj.data.shape_keys is None:
                obj.shape_key_add(name = base_name)
            elif base_name not in obj.data.shape_keys.key_blocks:
                obj.shape_key_add(name = base_name)
        if only_basis:
            return
            """
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

            """
    
    def __len__(self):
        return len(self.obj.data.vertices)
    
    
    def __getitem__(self, index):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system
        
        """
        from batoms.bond import Bond
        if isinstance(index, int):
            bond = Bond(self.label, index, bonds=self)
            # bpy.ops.object.mode_set(mode=mode)
            return bond
        else:
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
        s = "Bonds(Total: {:6d}, {}" .format(len(self), self.arrays)
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


def build_bondlists(species, positions, cell, pbc, cutoff):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from batoms.neighborlist import neighbor_kdtree
    if len(cutoff) == 0: return {}
    #
    tstart = time()
    # nli, nlj, nlS = primitive_neighbor_list('ijS', pbc,
    #                                cell,
    #                                positions, cutoff, species=species,
    #                                self_interaction=False,
    #                                max_nbins=1e6)
    nli, nlj, nlS = neighbor_kdtree('ijS', species, 
                positions, cell, pbc,
            cutoff)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    bondlists1 = bondlists.copy()
    bondlists1[:, [0, 1]] = bondlists1[:, [1, 0]]
    bondlists2 = np.concatenate((bondlists, bondlists1), axis=0)
    np.unique(bondlists2)
    return bondlists

def calc_bond_data(speciesarray, positions, cell, radii, 
            bondlists, bondsetting,
            model_styles):
    """
    """
    from ase.data import chemical_symbols
    from batoms.tools import calc_euler_angle

    tstart = time()
    chemical_symbols = np.array(chemical_symbols)
    # properties
    atoms_index1 = np.array(bondlists[:, 0])
    atoms_index2 = np.array(bondlists[:, 1])
    atoms_index3 = np.roll(bondlists[:, 0], 1)
    atoms_index4 = np.roll(bondlists[:, 1], 1)
    nb =len(bondlists)
    orders = np.zeros(nb, dtype = int) # 1, 2, 3
    styles = np.zeros(nb, dtype = int) # 0, 1, 2
    widths = np.ones(nb, dtype = float) 
    species_index1 = np.zeros(nb, dtype = int)
    species_index2 = np.zeros(nb, dtype = int)
    model_styles = model_styles[bondlists[:, 0]]
    if nb == 0:
        return default_bond_datas
    #------------------------------------
    # offsets
    offsets = np.dot(bondlists[:, 2:5], cell)
    # bond vectors and lengths
    centers = (positions[bondlists[:, 0]] + 
            positions[bondlists[:, 1]] + offsets)/2.0
    #---------------------------------------------
    for b in bondsetting:
        spi = b.species1
        spj = b.species2
        indi = (speciesarray[bondlists[:, 0]] == spi)
        indj = (speciesarray[bondlists[:, 1]] == spj)
        ind = indi & indj
        orders[ind] = b.order
        styles[ind] = int(b.style)
        widths[ind] = b.width
        species_index1[ind] = string2Number(spi)
        species_index2[ind] = string2Number(spj)
        if b.order >1:
            atoms_index3, atoms_index4 = secondBond(b, speciesarray, atoms_index3, atoms_index4, bondlists)
            # offsets = high_order_bond_plane(b, speciesarray, positions, nvecs, offsets, bondlists)
    # eulers = calc_euler_angle(offsets, nvecs)
    datas = {
        'atoms_index1': atoms_index1,
        'atoms_index2': atoms_index2,
        'atoms_index3': atoms_index3,
        'atoms_index4': atoms_index4,
        'species_index1': species_index1,
        'species_index2': species_index2,
        'centers':centers,
        # 'vectors':nvecs,
        'offsets':offsets,
        # 'eulers':eulers,
        # 'lengths':lengths,
        'widths':widths,
        'orders':orders,
        'styles':styles,
        'model_styles':model_styles,
    }
    # print('datas: ', datas)
    print('calc_bond_data: {0:10.2f} s'.format(time() - tstart))
    return datas


def secondBond(b, speciesarray, atoms_index3, atoms_index4, bondlists):
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
            atoms_index3[i] = 0
            atoms_index4[i] = 1
        else:
            atoms_index3[i] = localbondlist[0, 0]
            atoms_index4[i] = localbondlist[0, 1]

    return atoms_index3, atoms_index4

def high_order_bond_plane(b, speciesarray, positions, nvec, offsets, bondlists):
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
            second_bond = positions[localbondlist[0, 0]] - positions[localbondlist[0, 1]]
        norml = np.cross(second_bond, nvec[i]) + np.array([1e-8, 0, 0])
        offset = np.cross(norml, nvec[i]) + np.array([1e-8, 0, 0])
        offsets[i] = offset/np.linalg.norm(offset)
        # print(offsets[i], offset)

    return offsets