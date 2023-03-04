"""Definition of the polyhedras class.

This module defines the polyhedras object in the Batoms package.

"""

import bpy
from time import time
from batoms.attribute import Attributes
from batoms.utils.butils import object_mode, compareNodeType, get_nodes_by_name
from batoms.utils import string2Number
import numpy as np
from batoms.base.object import ObjectGN
from .setting import PolyhedraSettings
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


default_attributes = [
    {"name":'atoms_index1', "data_type": 'INT', "dimension": 0},
    {"name":'atoms_index2', "data_type": 'INT', "dimension": 0},
    {"name":'species_index', "data_type": 'INT', "dimension": 0},
    {"name":'face_species_index', "data_type": 'INT', "dimension": 0, "domain": 'FACE'},
    {"name":'show', "data_type": 'INT', "dimension": 0},
    {"name":'style', "data_type": 'INT', "dimension": 0},
]


default_GroupInput = [
    ['atoms_index1', 'NodeSocketInt', 'POINT'],
    ['atoms_index2', 'NodeSocketInt', 'POINT'],
    ['species_index', 'NodeSocketInt', 'POINT'],
    ['face_species_index', 'NodeSocketInt', 'FACE'],
    ['show', 'NodeSocketBool', 'POINT'],
    ['style', 'NodeSocketInt', 'POINT'],
]

default_polyhedra_datas = {
    'atoms_index1': np.ones(0, dtype=int),
    'atoms_index2': np.ones(0, dtype=int),
    'species_index': np.ones(0, dtype=int),
    'face_species_index': np.ones(0, dtype=int),
    'vertices': np.zeros((0, 3)),
    'offsets': np.zeros((0, 3)),
    'edges': [],
    'faces': [],
    'widths': np.ones(0, dtype=float),
    'shows': np.zeros(0, dtype=int),
}


class Polyhedra(ObjectGN):
    """Polyhedras Class

    A Polyhedras object is linked to this main collection in Blender.

    Parameters:

    label: str
        Name of the Bpolyhedras.
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
                 label=None,
                 polyhedra_datas=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 ):
        #
        self.batoms = batoms
        self.label = label
        name = 'polyhedra'
        ObjectGN.__init__(self, label, name)
        self.settings = PolyhedraSettings(
            self.label, batoms=batoms, polyhedras=self)
        flag = self.load()
        if not flag and polyhedra_datas is not None:
            self.build_object(polyhedra_datas)
        else:
            self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)


    def build_object(self, datas, attributes={}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        tstart = time()
        if len(datas['vertices'].shape) == 2:
            self._frames = {'vertices': np.array([datas['vertices']]),
                            'offsets': np.array([datas['offsets']]),
                            }
            vertices = datas['vertices']
            offsets = datas['offsets']
        elif len(datas['vertices'].shape) == 3:
            self._frames = {'vertices': datas['vertices'],
                            'offsets': datas['offsets'],
                            }
            vertices = datas['vertices'][0]
            offsets = datas['offsets'][0]
        else:
            raise Exception('Shape of vertices is wrong!')

        attributes.update({
            'atoms_index1': datas['atoms_index1'],
            'atoms_index2': datas['atoms_index2'],
            'species_index': datas['species_index'],
            'face_species_index': datas['face_species_index'],
            'show': datas['shows'],
            # 'model_style': datas['model_styles'],
        })
        name = self.obj_name
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(
            vertices, datas['edges'], datas['faces'])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        # Add attributes
        self._attributes.add(default_attributes)
        self.settings.coll.objects.link(obj)
        obj.batoms.type = 'POLYHEDRA'
        obj.batoms.label = self.label
        obj.parent = self.batoms.obj
        #
        name = '%s_polyhedra_offset' % self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.settings.coll.objects.link(obj)
        obj.hide_set(True)
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis=True)
        # self.assign_materials()
        self.update_geometry_node_material()
        obj.parent = self.obj
        logger.debug('polyhedras: build_object: {0:10.2f} s'.format(time() - tstart))

    def assign_materials(self):
        # sort element by occu
        me = self.obj.data
        me.materials.clear()
        for sp in self.settings:
            me.materials.append(self.settings.materials[sp.species][0])

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        from batoms.utils.butils import get_nodes_by_name, build_modifier
        tstart = time()
        name = 'GeometryNodes_%s_polyhedra' % self.label
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
        # ------------------------------------------------------------------
        # calculate bond vector, length, rotation based on the index
        # Get four positions from batoms, bond and the second bond
        #  for high order bond plane
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
                                TransferBatoms.inputs[3])
        gn.node_group.links.new(GroupInput.outputs[2],
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
        InputIndex = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_InputIndex' % self.label,
                                        'GeometryNodeInputIndex')
        gn.node_group.links.new(ObjectOffsets.outputs['Geometry'],
                                TransferOffsets.inputs[0])
        gn.node_group.links.new(PositionOffsets.outputs['Position'],
                                TransferOffsets.inputs[3])
        gn.node_group.links.new(InputIndex.outputs[0],
                                TransferOffsets.inputs["Index"])
        # we need one add operation to get the positions with offset
        VectorAdd = get_nodes_by_name(gn.node_group.nodes,
                                      '%s_VectorAdd' % (self.label),
                                      'ShaderNodeVectorMath')
        VectorAdd.operation = 'ADD'
        gn.node_group.links.new(TransferBatoms.outputs[2],
                                VectorAdd.inputs[0])
        gn.node_group.links.new(TransferOffsets.outputs[2],
                                VectorAdd.inputs[1])
        # set positions
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label,
                                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                SetPosition.inputs['Geometry'])
        gn.node_group.links.new(VectorAdd.outputs[0],
                                SetPosition.inputs['Position'])
        gn.node_group.links.new(SetPosition.outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])

        # find kinds by the names of species
        for sp in self.settings:
            self.add_geometry_node(sp.as_dict())
        logger.debug('Build geometry nodes for polyhedras: %s' % (time() - tstart))

    def add_geometry_node(self, sp):
        from batoms.utils.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        GroupOutput = gn.node_group.nodes[1]
        previousNode = GroupOutput.inputs['Geometry'].links[0].from_socket
        # print(previousNode)
        # we need two compares for one species,
        # vertices
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_CompareSpecies_%s_vertex' % (
                                               self.label, sp["species"]),
                                           compareNodeType)
        CompareSpecies.operation = 'EQUAL'
        CompareSpecies.inputs[1].default_value = string2Number(sp["species"])
        gn.node_group.links.new(
            GroupInput.outputs[3], CompareSpecies.inputs[0])
        # face
        CompareSpeciesFace = get_nodes_by_name(gn.node_group.nodes,
                                               '%s_CompareSpecies_%s_face' % (
                                                   self.label, sp["species"]),
                                               compareNodeType)
        CompareSpeciesFace.operation = 'EQUAL'
        CompareSpeciesFace.inputs[1].default_value = string2Number(
            sp["species"])
        gn.node_group.links.new(
            GroupInput.outputs[4], CompareSpeciesFace.inputs[0])
        #
        setMaterialIndex = get_nodes_by_name(gn.node_group.nodes,
                                             '%s_setMaterialIndex_%s' % (
                                                 self.label, sp["species"]),
                                             'GeometryNodeSetMaterialIndex')
        setMaterialIndex.inputs[2].default_value = \
            self.settings.materials[sp["species"]][1]
        gn.node_group.links.new(previousNode,
                                setMaterialIndex.inputs['Geometry'])
        gn.node_group.links.new(CompareSpeciesFace.outputs[0],
                                setMaterialIndex.inputs['Selection'])
        gn.node_group.links.new(setMaterialIndex.outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])

    def set_realize_instances(self, realize_instances):
        """ Make instancing object real
        # TODO: add make real to geometry node
        """
        #
        i = len(self.settings.bpy_setting) - 1
        sp = self.settings.bpy_setting[i]
        setMaterialIndex = get_nodes_by_name(self.gnodes.node_group.nodes,
                                            '%s_setMaterialIndex_%s' % (
                                                self.label, sp["species"]),
                                            'GeometryNodeSetMaterialIndex')

        nodes = self.gnodes.node_group.nodes
        RealizeInstances = get_nodes_by_name(self.gnodes.node_group.nodes,
                                             '%s_RealizeInstances' % self.label,
                                             'GeometryNodeRealizeInstances')
        if not realize_instances:
            # switch off
            if len(RealizeInstances.outputs[0].links) > 0:
                link = RealizeInstances.outputs[0].links[0]
                self.gnodes.node_group.links.remove(link)
            self.gnodes.node_group.links.new(setMaterialIndex.outputs[0],
                nodes[1].inputs[0])
        else:
            self.gnodes.node_group.links.new(setMaterialIndex.outputs[0],
                RealizeInstances.inputs[0])
            self.gnodes.node_group.links.new(
                RealizeInstances.outputs[0],
                nodes[1].inputs[0])
        self.gnodes.node_group.update_tag()

    def update(self, ):
        """Draw polyhedras.
        calculate bond in all farmes, and save
        get the max number of polyhedras for each pair
        draw the polyhedras
        add shape key
        the extra polyhedras use the find bond data.
        """
        # if not self.bondlist:
        object_mode()
        # clean_coll_objects(self.coll, 'bond')
        frames = self.batoms.get_frames()
        arrays = self.batoms.arrays
        show = arrays['show'].astype(bool)
        species = arrays['species']
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        polyhedra_datas = {}
        tstart = time()
        for f in range(nframe):
            logger.debug('update polyhedra: {}'.format(f))
            positions = frames[f, show, :]
            bondlists = self.bondlists
            if len(bondlists) == 0:
                self.batoms.bond.update()
                bondlists = self.batoms.bond.bondlists
            polyhedra_kinds = \
                self.calc_polyhedra_data(bondlists,
                                         species,
                                         positions,
                                         arrays['model_style'][show])
            if f == 0:
                polyhedra_datas = polyhedra_kinds

        if len(polyhedra_datas) == 0:
            return
        self.set_arrays(polyhedra_datas)
        # self.coll.objects.link(bb.obj)
        # bpy.data.collections['Collection'].objects.unlink(bb.obj)
        # bb.set_frames()
        # bpy.context.scene.frame_set(self.batoms.nframe)
        logger.debug('draw polyhedra: {0:10.2f} s'.format(time() - tstart))

    def update_geometry_node_material(self):
        """
        Make sure all species has a geometry node flow
        and the material are updated.
        """
        tstart = time()
        for sp in self.settings:
            sp = sp.as_dict()
            # update material in geometry node
            self.settings.build_materials(sp)
        self.assign_materials()
        for i in range(len(self.settings.bpy_setting)):
            sp = self.settings.bpy_setting[i]
            setMaterialIndex = get_nodes_by_name(self.gnodes.node_group.nodes,
                                                '%s_setMaterialIndex_%s' % (
                                                    self.label, sp["species"]),
                                                'GeometryNodeSetMaterialIndex')
            setMaterialIndex.inputs[2].default_value = i
        logger.debug('update polyhedra materials: %s'%(time() - tstart))

    @property
    def obj_o(self):
        return self.get_obj_o()

    def get_obj_o(self):
        name = '%s_polyhedra_offset' % self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.' % name)
        return obj_o

    def get_arrays(self):
        """
        """
        object_mode()
        # tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions[0],
                       'offsets': self.positions[1],
                       })
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    @property
    def bondlists(self):
        return self.get_bondlists()

    def get_bondlists(self):
        bondlists = self.batoms.bond.bondlists
        return bondlists

    def set_arrays(self, arrays):
        """
        """
        attributes = self.attributes
        # same length
        if len(arrays['vertices']) == len(attributes['show']):
            self.positions = arrays['vertices']
            self.offsets = arrays['offsets']
            self.set_attributes({'atoms_index1': arrays['atoms_index1']})
            self.set_attributes({'atoms_index2': arrays['atoms_index2']})
            self.set_attributes({'species_index': arrays['species_index']})
            self.set_attributes(
                {'face_species_index': arrays['face_species_index']})
            self.update_geometry_node_material()
        else:
            # add or remove vertices, rebuild the object
            self.build_object(arrays)


    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            self.build_object(default_polyhedra_datas)
        return obj

    @property
    def offsets(self):
        return self.get_offsets()

    def get_offsets(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        offsets = np.empty(n*3, dtype=np.float64)
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
        self.obj_o.data.shape_keys.key_blocks[0].data.foreach_set(
            'co', offsets)
        self.obj_o.data.update()

    def get_frames(self):
        """
        """
        frames = {}
        frames['vertices'] = self.get_obj_frames(self.obj)
        frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0:
            return
        name = '%s_polyhedra' % (self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['vertices'])
        #
        name = '%s_polyhedra_offset' % (self.label)
        obj = self.obj_o
        self.set_obj_frames(name, obj, frames['offsets'])

    def __getitem__(self, index):
        """Return a subset of the polyhedras.
        # TODO
        """
        pass

    def __setitem__(self, index, value):
        """set a subset of the polyhedras.
        # TODO
        """
        pass

    def calc_polyhedra_data(self, bondlists, species, positions,
                            model_styles):
        """
        """
        from ase.data import chemical_symbols
        from scipy.spatial import ConvexHull

        tstart = time()
        chemical_symbols = np.array(chemical_symbols)
        # find bonds contribute to polyhedra
        indices = bondlists[:, 8].astype(bool)
        bondlists = bondlists[indices]
        # maxinum number of poly
        npoly = len(bondlists)
        if npoly == 0:
            return default_polyhedra_datas
        # properties
        npoly = npoly*12
        atoms_index1 = np.zeros(npoly, dtype=int)  # 0, 1, 2
        atoms_index2 = np.zeros(npoly, dtype=int)  # 0, 1, 2
        face_species_index = np.zeros(npoly, dtype=int)  # 0, 1, 2
        widths = np.ones(npoly, dtype=float)
        species_index = np.zeros(npoly, dtype=int)
        # ------------------------------------
        # offsets10 = bondlists[:, 2:5]
        offsets20 = bondlists[:, 5:8]
        indices10 = bondlists[:, 0].astype(int)
        indices20 = bondlists[:, 1].astype(int)
        positions20 = positions[indices20] + offsets20
        # ---------------------------------------------
        vertices = []
        offsets = []
        edges = []
        faces = []
        nv = 0
        nf = 0
        speciesarray = species[indices10]
        for poly in self.settings:
            # find center atoms == species
            spis = np.where(speciesarray == poly.species)[0]
            indices1 = indices10[spis]
            bondlist1 = bondlists[spis]
            # center atoms is define by i and the offset
            atom_centers = bondlist1[:, [0, 2, 3, 4]]
            positions21 = positions20[spis]
            # find indices of center atoms
            u, indices = np.unique(atom_centers, axis=0, return_inverse=True)
            indices = np.append(indices, len(indices1))
            n = len(u)
            # print(n)
            for i in range(n):
                if model_styles[int(u[i, 0])] != 2:
                    continue
                mask = np.where(indices == i)[0]
                vertices1 = positions21[mask]
                dnv = len(vertices1)
                if dnv >= 4:
                    # search face for polyhedra
                    hull = ConvexHull(vertices1)
                    face = hull.simplices
                    # tri = Delaunay(vertices1)
                    # face = tri.simplices
                    dnf = len(face)
                    face = face + nv
                    edge = []
                    for f in face:
                        edge.append([f[0], f[1]])
                        edge.append([f[0], f[2]])
                        edge.append([f[1], f[2]])
                    edges = edges + list(edge)
                    faces = faces + list(face)
                    vertices.extend(vertices1)
                    offsets.extend(bondlist1[mask, 5:8])
                    nv1 = nv + dnv
                    nf1 = nf + dnf
                    atoms_index1[nv:nv1] = bondlist1[mask, 0]
                    atoms_index2[nv:nv1] = bondlist1[mask, 1]
                    widths[nv:nv1] = poly.width
                    species_index[nv:nv1] = string2Number(poly.species)
                    face_species_index[nf:nf1] = string2Number(poly.species)
                    nv = nv1
                    nf = nf1
        n = len(vertices)
        if n == 0:
            datas = default_polyhedra_datas
        else:
            nf = len(faces)
            shows = np.ones(n, dtype=int)
            datas = {
                'atoms_index1': atoms_index1[0:n],
                'atoms_index2': atoms_index2[0:n],
                'species_index': species_index[0:n],
                'face_species_index': face_species_index[0:nf],
                'vertices': np.array(vertices),
                'offsets': np.array(offsets),
                'widths': widths[0:n],
                'edges': edges,
                'faces': faces,
                'shows': shows,
                'model_styles': model_styles,
            }
        # print('datas: ', datas)
        logger.debug('calc_polyhedra_data: {0:10.2f} s'.format(time() - tstart))
        return datas

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
        data['settings'] = self.settings.as_dict()
        return data
