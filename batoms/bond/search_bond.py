import bpy
import numpy as np
from time import time
from batoms.attribute import Attributes
from batoms.base.object import ObjectGN
from batoms.utils.butils import object_mode, compareNodeType
from batoms.utils import number2String, string2Number
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

default_attributes = [
    {"name": 'atoms_index', "data_type": 'INT', "dimension": 0},
    {"name": 'species_index', "data_type": 'INT', "dimension": 0},
    {"name": 'show', "data_type": 'INT', "dimension": 0},
    {"name": 'select', "data_type": 'INT', "dimension": 0},
    {"name": 'model_style', "data_type": 'INT', "dimension": 0},
    {"name": 'scale', "data_type": 'FLOAT', "dimension": 0},
]

default_GroupInput = [
    ['atoms_index', 'NodeSocketInt'],
    ['species_index', 'NodeSocketInt'],
    ['show', 'NodeSocketBool'],
    ['select', 'NodeSocketInt'],
    ['model_style', 'NodeSocketInt'],
    ['scale', 'NodeSocketFloat'],
]


default_search_bond_datas = {
    'atoms_index': np.ones(0, dtype=int),
    'species_index': np.ones(0, dtype=int),
    'species': np.ones(0, dtype='U4'),
    'positions': np.zeros((0, 3)),
    'scales': np.zeros(0),
    'offsets': np.zeros((0, 3)),
    'model_styles': np.ones(0, dtype=int),
    'shows': np.ones(0, dtype=int),
    'selects': np.ones(0, dtype=int),
}


class SearchBond(ObjectGN):

    def __init__(self,
                 label=None,
                 search_bond_datas=None,
                 batoms=None,
                 ):
        """SearchBond Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            search_bond_datas (_type_, optional): _description_. Defaults to None.
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = 'search_bond'
        ObjectGN.__init__(self, label, name)
        if search_bond_datas is not None:
            self.build_object(search_bond_datas)
        elif self.loadable():
            self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        else:
            self.build_object(default_search_bond_datas)

    def loadable(self):
        """Check loadable or not
        """
        # object exist
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            return False
        # batoms exist, and flag is True
        # coll = bpy.data.collections.get(self.label)
        # if coll is None:
            # return False
        # return coll.Bbond.flag
        return True

    def build_object(self, datas, attributes={}):
        """
        build child object and add it to main objects.
        """
        # tstart = time()
        if len(datas['positions'].shape) == 2:
            self._frames = {'positions': np.array([datas['positions']]),
                            'offsets': np.array([datas['offsets']]),
                            }
            positions = datas['positions']
            offsets = datas['offsets']
        elif len(datas['positions'].shape) == 3:
            self._frames = {'positions': datas['positions'],
                            'offsets': datas['offsets'],
                            }
            positions = datas['positions'][0]
            offsets = datas['offsets'][0]
        else:
            raise Exception('Shape of positions is wrong!')
        #
        attributes.update({
            'atoms_index': datas['atoms_index'],
            'species_index': datas['species_index'],
            'show': datas['shows'],
            'model_style': datas['model_styles'],
            'select': datas['selects'],
            'scale': datas['scales'],
        })
        name = '%s_search_bond' % self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(positions, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        # Add attributes
        self._attributes.add(default_attributes)
        self.batoms.coll.objects.link(obj)
        obj.batoms.type = 'BOND'
        obj.batoms.label = self.batoms.label
        obj.Bbond.label = self.batoms.label
        obj.parent = self.batoms.obj
        #
        name = '%s_search_bond_offset' % self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.batoms.coll.objects.link(obj)
        obj.hide_set(True)
        obj.parent = self.obj
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis=True)
        # print('boundary: build_object: {0:10.2f} s'.format(time() - tstart))

    def update(self, bondlist, mollists, moldatas, arrays, cell):
        self.hide = False
        search_bond_data = self.calc_search_bond_data(
            bondlist, mollists, moldatas, arrays, cell)
        self.set_arrays(search_bond_data)

    def build_geometry_node(self):
        """
        """
        from batoms.utils.butils import get_nodes_by_name, build_modifier
        name = 'GeometryNodes_%s_search_bond' % self.label
        modifier = build_modifier(self.obj, name)
        # ------------------------------------------------------------------
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = modifier.node_group.nodes[1]
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
        gn.node_group.links.new(
            GroupInput.outputs['Geometry'], JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(
            JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])
        # ------------------------------------------------------------------
        # transform postions of batoms to boundary
        ObjectBatoms = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_ObjectBatoms' % self.label,
                                         'GeometryNodeObjectInfo')
        ObjectBatoms.inputs['Object'].default_value = self.batoms.obj
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_PositionBatoms' % (self.label),
                                           'GeometryNodeInputPosition')
        TransferBatoms = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_TransferBatoms' % (self.label),
                                           'GeometryNodeSampleIndex')
        TransferBatoms.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectBatoms.outputs['Geometry'],
                                TransferBatoms.inputs[0])
        gn.node_group.links.new(PositionBatoms.outputs['Position'],
                                TransferBatoms.inputs[1])
        gn.node_group.links.new(GroupInput.outputs[1],
                                TransferBatoms.inputs['Index'])
        # ------------------------------------------------------------------
        # add positions with offsets
        # transfer offsets from object self.obj_o
        ObjectOffsets = get_nodes_by_name(gn.node_group.nodes,
                                          '%s_ObjectOffsets' % (self.label),
                                          'GeometryNodeObjectInfo')
        ObjectOffsets.inputs['Object'].default_value = self.obj_o
        PositionOffsets = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_PositionOffsets' % (
                                                self.label),
                                            'GeometryNodeInputPosition')
        TransferOffsets = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_TransferOffsets' % self.label,
                                            'GeometryNodeSampleIndex')
        TransferOffsets.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectOffsets.outputs['Geometry'],
                                TransferOffsets.inputs[0])
        gn.node_group.links.new(PositionOffsets.outputs['Position'],
                                TransferOffsets.inputs[1])
        OffsetNode = self.vectorDotMatrix(
            gn, TransferOffsets.outputs[2], self.batoms.cell, '')
        # we need one add operation to get the positions with offset
        VectorAdd = get_nodes_by_name(gn.node_group.nodes,
                                      '%s_VectorAdd' % (self.label),
                                      'ShaderNodeVectorMath')
        VectorAdd.operation = 'ADD'
        gn.node_group.links.new(TransferBatoms.outputs[2], VectorAdd.inputs[0])
        gn.node_group.links.new(OffsetNode.outputs[0], VectorAdd.inputs[1])
        # set positions
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label,
                                        'GeometryNodeSetPosition')
        gn.node_group.links.new(
            GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(
            VectorAdd.outputs[0], SetPosition.inputs['Position'])
        #
        # ------------------------------------------------------------------
        # transform scale of batoms to boundary
        if bpy.app.version_string >= '3.2.0':
            ScaleBatoms = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_ScaleBatoms' % (self.label),
                                            'GeometryNodeInputNamedAttribute')
            # need to be "FLOAT_VECTOR", because scale is "FLOAT_VECTOR"
            ScaleBatoms.data_type = "FLOAT_VECTOR"
            ScaleBatoms.inputs[0].default_value = "scale"
            TransferScale = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_TransferScale' % (self.label),
                                            'GeometryNodeSampleIndex')
            TransferScale.data_type = 'FLOAT_VECTOR'
            gn.node_group.links.new(ObjectBatoms.outputs['Geometry'],
                                    TransferScale.inputs[0])
            gn.node_group.links.new(ScaleBatoms.outputs['Attribute'],
                                    TransferScale.inputs[1])
            gn.node_group.links.new(GroupInput.outputs[1],
                                    TransferScale.inputs['Index'])

    def add_geometry_node(self, spname):
        """
        """
        from batoms.utils.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes,
                                           'CompareFloats_%s_%s' % (
                                               self.label, spname),
                                           compareNodeType)
        CompareSpecies.operation = 'EQUAL'
        # CompareSpecies.data_type = 'INT'
        CompareSpecies.inputs[1].default_value = string2Number(spname)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                                            'InstanceOnPoint_%s_%s' % (
                                                self.label, spname),
                                            'GeometryNodeInstanceOnPoints')
        ObjectInfo = get_nodes_by_name(gn.node_group.nodes,
                                       'ObjectInfo_%s_%s' % (
                                           self.label, spname),
                                       'GeometryNodeObjectInfo')
        ObjectInfo.inputs['Object'].default_value = \
            self.batoms.species.instancers[spname]
        BoolShow = get_nodes_by_name(gn.node_group.nodes,
                                     'BooleanMath_%s_%s_1' % (
                                         self.label, spname),
                                     'FunctionNodeBooleanMath')
        #
        gn.node_group.links.new(SetPosition.outputs['Geometry'],
                                InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[2],
                                CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3],
                                BoolShow.inputs[0])
        # transfer scale
        if bpy.app.version_string >= '3.2.0':
            TransferScale = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_TransferScale' % (self.label),
                                            'GeometryNodeSampleIndex')
            gn.node_group.links.new(
                TransferScale.outputs[2], InstanceOnPoint.inputs['Scale'])
        else:
            gn.node_group.links.new(GroupInput.outputs[6],
                                InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSpecies.outputs[0],
                                BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'],
                                InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'],
                                InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'],
                                JoinGeometry.inputs['Geometry'])

    def update_geometry_node_instancer(self, spname, instancer):
        """When instances are re-build, we need also update
        the geometry node.

        Args:
            spname (str): name of the species
        """
        from batoms.utils.butils import get_nodes_by_name
        # update  instancers
        ObjectInfo = get_nodes_by_name(self.gnodes.node_group.nodes,
                                       'ObjectInfo_%s_%s' % (
                                           self.label, spname),
                                       'GeometryNodeObjectInfo')
        ObjectInfo.inputs['Object'].default_value = instancer
        logger.debug('update boundary instancer: {}'.format(spname))

    @property
    def obj_o(self):
        return self.get_obj_o()

    def get_obj_o(self):
        name = '%s_search_bond_offset' % self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.' % name)
        return obj_o

    def get_search_bond(self):
        boundary = np.array(self.batoms.obj.batoms.boundary)
        return boundary.reshape(3, -1)

    def set_arrays(self, arrays):
        """
        """
        attributes = self.attributes
        # same length
        dnvert = len(arrays['species_index']) - \
            len(attributes['species_index'])
        if dnvert > 0:
            self.add_vertices_bmesh(dnvert)
            self.add_vertices_bmesh(dnvert, self.obj_o)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
            self.delete_vertices_bmesh(range(-dnvert), self.obj_o)
        if len(self) == 0:
            self.update_mesh()
            return
        self.positions = arrays["positions"][0]
        self.offsets = arrays["offsets"][0]
        self.set_frames(arrays)
        species_index = [string2Number(sp) for sp in arrays['species']]
        self.set_attributes({
                            'atoms_index': arrays["atoms_index"],
                            'species_index': species_index,
                            'scale': arrays['scales'],
                            'show': arrays['shows'],
                            })
        species = np.unique(arrays['species'])
        for sp in species:
            self.add_geometry_node(sp)

    def get_arrays(self):
        """
        """
        object_mode()
        # tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        arrays.update({'offsets': self.offsets})
        # radius
        radius = self.batoms.radius
        arrays.update({'radius': np.zeros(len(self))})
        species = np.array([number2String(i)
                           for i in arrays['species_index']], dtype='U20')
        arrays['species'] = species
        for sp, value in radius.items():
            mask = np.where(arrays['species'] == sp)
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

    def get_search_bond_data(self, include_batoms=False):
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
        offsets = np.empty(n*3, dtype=int)
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
        n = len(self.obj_o.data.vertices)
        if len(offsets) != n:
            raise ValueError('offsets has wrong shape %s != %s.' %
                             (len(offsets), n))
        if n == 0:
            return
        offsets = offsets.reshape((n*3, 1))
        if self.obj_o.data.shape_keys is None and len(self) > 0:
            base_name = "Basis_%s"%self.obj_o.name
            self.obj_o.shape_key_add(name=base_name)
        self.obj_o.data.shape_keys.key_blocks[0].data.foreach_set(
            'co', offsets)
        self.obj_o.data.update()
        self.update_mesh(self.obj_o)

    @property
    def bondlists(self):
        return self.get_bondlists()

    def get_bondlists(self):
        bondlists = self.batoms.bond.arrays
        return bondlists

    def get_frames(self):
        """
        """
        frames = {}
        frames['positions'] = self.get_obj_frames(self.obj)
        frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0:
            return
        name = '%s_search_bond' % (self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['positions'])
        #
        name = '%s_search_bond_offset' % (self.label)
        obj = self.obj_o
        self.set_obj_frames(name, obj, frames['offsets'])

    def calc_search_bond_data(self, bondlists,
                              mollists, moldatas, arrays, cell):
        """
        """
        # tstart = time()
        # 9th column of bondlists is search bond type 1
        # search type 1: atoms search by bond, the one with offset
        indices1 = np.where(
            (bondlists[:, 2:5] != np.array([0, 0, 0])).any(axis=1))[0]
        indices2 = np.where(
            (bondlists[:, 5:8] != np.array([0, 0, 0])).any(axis=1))[0]
        n = len(indices1) + len(indices2)
        if n == 0:
            return default_search_bond_datas
        bondlists1 = bondlists[indices1]
        bondlists1 = bondlists1[:, [0, 2, 3, 4]]
        bondlists1 = np.unique(bondlists1, axis=0)
        indices1 = bondlists1[:, 0].astype(int)
        model_styles1 = arrays['model_style'][indices1]
        shows1 = arrays['show'][indices1]
        selects1 = arrays['select'][indices1]
        scales1 = arrays['scale'][indices1]
        species_indexs1 = arrays['species_index'][indices1]
        species1 = arrays['species'][indices1]
        offset_vectors1 = bondlists1[:, 1:4]
        offsets1 = np.dot(offset_vectors1, cell)
        positions1 = arrays['positions'][indices1] + offsets1
        # ------------------------------------
        #
        bondlists2 = bondlists[indices2]
        bondlists2 = bondlists2[:, [1, 5, 6, 7]]
        bondlists2 = np.unique(bondlists2, axis=0)
        indices2 = bondlists2[:, 0].astype(int)
        model_styles2 = arrays['model_style'][indices2]
        shows2 = arrays['show'][indices2]
        selects2 = arrays['select'][indices2]
        scales2 = arrays['scale'][indices2]
        species_indexs2 = arrays['species_index'][indices2]
        species2 = arrays['species'][indices2]
        offset_vectors2 = bondlists2[:, 1:4]
        offsets2 = np.dot(offset_vectors2, cell)
        positions2 = arrays['positions'][indices2] + offsets2
        #
        indices = np.append(indices1, indices2)
        species_indexs = np.append(species_indexs1, species_indexs2)
        species = np.append(species1, species2)
        positions = np.append(positions1, positions2, axis=0)
        offset_vectors = np.append(offset_vectors1, offset_vectors2, axis=0)
        model_styles = np.append(model_styles1, model_styles2)
        selects = np.append(selects1, selects2)
        shows = np.append(shows1, shows2)
        scales = np.append(scales1, scales2)
        datas = {
            'atoms_index': np.array(indices),
            'species_index': species_indexs,
            'species': species,
            'positions': np.array([positions]),
            # 'offsets':offsets,
            'offsets': np.array([offset_vectors]),
            'model_styles': model_styles,
            'shows': shows,
            'selects': selects,
            'scales': scales,
        }
        # =========================
        # # search molecule
        # atoms_index = []
        # for b in mollists:
        #     indices = moldatas[b[0]]
        #     datas['atoms_index'] = np.append(datas['atoms_index'], indices)
        #     datas['offsets'] = np.append(datas['offsets'], b[5:8])

        # print('datas: ', datas)
        # print('calc_search_bond_data: {0:10.2f} s'.format(time() - tstart))
        return datas
