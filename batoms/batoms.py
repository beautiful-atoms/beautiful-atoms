"""Definition of the Batoms class in the batoms package.

# TODO: get evaluated positiosn
# TODO: add location for boundary, bond and all child objects
# TODO: add feature: cavity
"""
import bpy
from batoms.bspecies import Bspecies
from batoms.volumetric_data import VolumetricData
from batoms.attribute import Attributes
from batoms.cell import Bcell
from batoms.bselect import Selects
from batoms.base.collection import BaseCollection
from batoms.base.object import ObjectGN
from batoms.ribbon.ribbon import Ribbon
from batoms.utils.butils import object_mode, show_index, \
    get_nodes_by_name
from batoms.utils import string2Number, read_from_others, deprecated
from batoms.plugins import plugin_info

import numpy as np
from time import time

import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

default_attributes = [
    {"name": 'select', "data_type": 'INT', "dimension": 0},
    {"name": 'species_index', "data_type": 'INT', "dimension": 0},
    {"name": 'species', "data_type": 'STRING', "dimension": 0},
    {"name": 'show', "data_type": 'INT', "dimension": 0},
    {"name": 'scale', "data_type": 'FLOAT', "dimension": 0},
    {"name": 'model_style', "data_type": 'INT', "dimension": 0},
]


default_GroupInput = [
    ['select', 'NodeSocketInt'],
    ['species_index', 'NodeSocketInt'],
    ['show', 'NodeSocketBool'],
    ['scale', 'NodeSocketFloat'],
]

subcollections = ['instancer', 'surface', 'ribbon', 'plane']


class Batoms(BaseCollection, ObjectGN):
    """_summary_

    Args:
        BaseCollection (_type_): _description_
        ObjectGN (_type_): _description_
    """

    def __init__(self,
                 label='batoms',
                 species=None,
                 positions=None,
                 attributes=None,
                 species_props=None,
                 pbc=False, cell=np.array([0, 0, 0]),
                 location=np.array([0, 0, 0]),
                 volume=None,
                 info=None,
                 show_unit_cell=True,
                 scale=1.0,
                 model_style=0,
                 radius_style='0',
                 color_style='0',
                 from_ase=None,
                 from_pymatgen=None,
                 from_pybel=None,
                 movie=True,
                 segments=None,
                 ):
        """Batoms Class
        The Batoms object is a interface to a batoms object in Blender.

        Args:
            label (str, optional):
                Name for the object in Blender. Defaults to 'batoms'.
            species (list, optional):
                List of species, e.g. ['O', 'H', 'H'] or ['Fe_u', 'Fe_d', 'O'].
                Defaults to [].
            positions (list, optional):
                Positions of atoms. Defaults to [].
            attributes (dict, optional):
                Properties for each atoms. Defaults to {}.
            species_props (dict, optional):
                Properties for each species. Defaults to {}.
            info (dict, optional):
                Other properties, e.g. secondary structure of protein.
                Defaults to {}.
            pbc (bool or list of bool, optional):
                Periodic boundary conditions. Defaults to False.
            cell (3x3 matrix or length 3 or 6 vector, optional):
                Unit cell. Defaults to None.
            location (array, optional):
                The location of whole Batoms object, the same as
                the origin of cell.
                Defaults to np.array([0, 0, 0]).
            show_unit_cell (bool, optional):
                Show unit cell during rendering. Defaults to True.
            volume (_type_, optional):
                Volumetric data. e.g. electron density.
                Defaults to None.
            scale (float, optional):
                Scale factor for all atoms. Defaults to 1.0.
            model_style (int, optional):
                Enum in [0, 1, 2, 3]. Defaults to 0.
            from_ase (ASE Atoms, optional):
                Import structure from ASE. Defaults to None.
            from_pymatgen (Pymatgen structure, optional):
                Import structure from pymatgen. Defaults to None.
            metaball (bool, optional):
                Show atoms as metaball. Defaults to False.
            movie (bool, optional):
                Load all frames. Defaults to True.
            segments (_type_, optional):
                Resolution of the sphere. Defaults to None.

        Examples:

        >>> from batoms import Batoms
        >>> h2o = Batom('h2o', species = ['O', 'H', 'H'],
                    positions= [[0, 0, 0.40],
                                [0, -0.76, -0.2],
                                [0, 0.76, -0.2]])
        >>> h2o.pbc = True
        """
        #
        self.obj_name = label
        ObjectGN.__init__(self, label)
        BaseCollection.__init__(self, coll_name=label)
        self._render = None
        self._bond = None
        self._polyhedra = None
        self._boundary = None
        self._volumetric_data = None
        #
        if from_ase or from_pymatgen or from_pybel:
            species, positions, attributes, cell, pbc, info = read_from_others(from_ase,
                                                                               from_pymatgen,
                                                                               from_pybel)
        if species is None and self.check_batoms(label):
            self.from_batoms(label)
        else:
            if species is None:
                species = []
                positions = []
            self.set_collection(label, color_style=color_style, radius_style=radius_style)
            self._cell = Bcell(label, cell, batoms=self)
            positions = np.array(positions)
            if len(positions.shape) == 3:
                self._frames = {"positions": positions}
                positions = self._frames["positions"][0]
            else:
                self._frames = {"positions": np.array([positions])}
            # init attributes
            natom = len(positions)
            self.coll.batoms.scale = scale
            self.coll.batoms.model_style = str(model_style)
            arrays = {
                      'species': species,
                      'species_index': [string2Number(sp) for sp in species],
                      'scale': np.ones(natom)*scale,
                      'show': np.ones(natom, dtype=int),
                      'model_style': np.ones(natom, dtype=int)*model_style,
                      'select': np.zeros(natom, dtype=int),
                      }
            if attributes is not None:
                arrays.update(attributes)
            self.build_object(label, positions, arrays, location)
            self.selects = Selects(label, self)
            if species_props is None:
                species_props = species
            self._species = Bspecies(
                label, species_props, self, segments=segments)
            self.selects.add('all', np.arange(len(self)))
            self._volumetric_data = VolumetricData(label, volume, self)
            self.set_pbc(pbc)
            # self.label = label
            if movie:
                self.set_frames()
            self.show_unit_cell = show_unit_cell

        self.ribbon = Ribbon(self.label, batoms=self, datas=info, update=True)
        show_index()
        self.hideOneLevel()

    def set_collection(self, label, color_style='0',
                        radius_style='0'):
        """Build main collection and its child collections.

        Args:
            label (str): name of Batoms object

        Raises:
            Exception: The label is already in use.
        """
        if bpy.data.collections.get(label):
            raise Exception("Failed, the name %s already in use!" % label)
        coll = bpy.data.collections.new(label)
        bpy.context.scene.collection.children.link(coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s' % (label, sub_name))
            coll.children.link(subcoll)
        coll.batoms.type = 'BATOMS'
        coll.batoms.label = label
        # set internal data
        self.coll.batoms.color_style = str(color_style)
        self.coll.batoms.radius_style = str(radius_style)

    def hideOneLevel(self):
        """Hide one level of collecitons in the outline in Blender
        """
        from batoms.utils.butils import hideOneLevel
        # hideOneLevel()
        pass

    def build_object(self, label, positions, arrays, location=[0, 0, 0]):
        """Build the main Batoms object

        Args:
            label (str):
                Name of the object
            arrays (array):
                arrays of properties for each atoms
            location (list, optional):
                Location of the object. Defaults to [0, 0, 0].
        """
        self.delete_obj(label)
        mesh = bpy.data.meshes.new(label)
        obj = bpy.data.objects.new(label, mesh)
        obj.data.from_pydata(positions, [], [])
        obj.location = location
        obj.batoms.type = 'BATOMS'
        obj.batoms.label = label
        self.coll.objects.link(obj)
        # add shape_keys
        if self.obj.data.shape_keys is None:
            self.obj.shape_key_add(name="Basis_{}".format(self.label))
        #
        self._attributes = Attributes(label=label, parent=self, obj_name=self.obj_name)
        # Add attributes
        self._attributes.add(default_attributes)
        # add cell object as its child
        self.cell.obj.parent = self.obj
        self.set_attributes(arrays)
        self.build_geometry_node()
        logger.info("Add object {}".format(self.label))

    def build_geometry_node(self):
        """Geometry node for instancing sphere on vertices!
        """
        from batoms.utils.butils import build_modifier
        name = 'GeometryNodes_%s' % self.label
        modifier = build_modifier(self.obj, name)
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = modifier.node_group.nodes[1]
        for att in default_GroupInput:
            inputs.new(att[1], att[0])
            id = inputs[att[0]].identifier
            modifier['%s_use_attribute' % id] = True
            modifier['%s_attribute_name' % id] = att[0]
        gn = modifier
        # print(gn.name)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        SeparateGeometry = \
            get_nodes_by_name(gn.node_group.nodes,
                              '%s_SeparateGeometry' % self.label,
                              'GeometryNodeSeparateGeometry')
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                SeparateGeometry.inputs['Geometry'])
        gn.node_group.links.new(GroupInput.outputs[3],
                                SeparateGeometry.inputs['Selection'])
        gn.node_group.links.new(SeparateGeometry.outputs[0],
                                JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])
        # set positions
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_PositionBatoms' % (self.label),
                                           'GeometryNodeInputPosition')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label,
                                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                SetPosition.inputs['Geometry'])
        # TODO: use cell object directly.
        cell = self.cell.array
        if np.isclose(np.linalg.det(self.cell.array), 0):
            cell = np.eye(3)
        icell = np.linalg.inv(cell)
        scaledPositionsNode = self.vectorDotMatrix(
            gn, PositionBatoms.outputs[0], icell, 'scaled')
        VectorWrap = get_nodes_by_name(gn.node_group.nodes,
                                       '%s_VectorWrap' % (self.label),
                                       'ShaderNodeVectorMath')
        VectorWrap.operation = 'WRAP'
        VectorWrap.inputs[1].default_value = [1, 1, 1]
        VectorWrap.inputs[2].default_value = [0, 0, 0]
        gn.node_group.links.new(
            scaledPositionsNode.outputs[0], VectorWrap.inputs['Vector'])
        PositionsNode = self.vectorDotMatrix(gn, VectorWrap.outputs[0], cell, '')
        gn.node_group.links.new(
            PositionsNode.outputs[0], SetPosition.inputs['Position'])
        self.wrap = self.pbc

    def add_geometry_node(self, spname, instancer):
        """Add geometry node for each species

        Args:
            spname (str):
                Name of the species
            instancer (bpy.data.object):
                Object to be instanced
        """
        from batoms.utils.butils import compareNodeType
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                        '%s_SetPosition' % self.label)
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
        ObjectInfo.inputs['Object'].default_value = instancer
        BoolShow = get_nodes_by_name(gn.node_group.nodes,
                                     'BooleanMath_%s_%s_1' % (
                                         self.label, spname),
                                     'FunctionNodeBooleanMath')
        # set materials
        SetMaterial = get_nodes_by_name(gn.node_group.nodes,
                                            'SetMaterial_%s_%s' % (
                                                self.label, spname),
                                            'GeometryNodeSetMaterial')
        SetMaterial.inputs[2].default_value = instancer.data.materials[0]
        #
        gn.node_group.links.new(SetPosition.outputs['Geometry'],
                                InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[2],
                                CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[4],
                                InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSpecies.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'],
                                InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'],
                                InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'],
                                SetMaterial.inputs['Geometry'])
        gn.node_group.links.new(SetMaterial.outputs['Geometry'],
                                JoinGeometry.inputs['Geometry'])

    def check_batoms(self, label):
        """Check batoms exist or not

        Args:
            label (str):
                Name of the Batoms

        Returns:
            bool: Batoms label exist or not
        """
        flag = True
        if label not in bpy.data.collections:
            flag = False
        elif bpy.data.collections[label].batoms.type == 'OTHER':
            flag = False
        if label not in bpy.data.objects:
            flag = False
        elif bpy.data.objects[label].batoms.type == 'OTHER':
            flag = False
        return flag

    def from_batoms(self, label):
        """Load a Batoms object

        Args:
            label (str):
                Name of the Batoms
        """
        self.coll_name = label
        self.obj_name = label
        self._cell = Bcell(label=label, batoms = self)
        self._species = Bspecies(label, {}, self)
        self._volumetric_data = VolumetricData(label, None, self)
        self._attributes = Attributes(label=label, parent=self, obj_name=self.obj_name)
        self.selects = Selects(label, self)

    @classmethod
    def from_dict(cls, data):
        """Build batoms object from a dictionary.

        Args:
            data (dict): _description_


        Returns:
            Batoms: _description_

        """
        pass

    def as_dict(self):
        """
        """
        data = {
            'batoms': {'label': self.label},
            'cell': None,
            'bond': {},
            'polyhedra': {},
            'boundary': {},
        }
        data['batoms']['array'] = dict(self.arrays)
        data['batoms']['species'] = self.species.as_dict()
        data['batoms'].update(self.coll.batoms.as_dict())
        data['cell'] = self.cell.as_dict()
        data['bond'] = self.bond.as_dict()
        data['polyhedra'] = self.polyhedra.as_dict()
        data['boundary'] = self.boundary.as_dict()
        # plugins
        for key, info in plugin_info.items():
            _name = '_{}'.format(key)
            if len(info) == 0:
                continue
            if getattr(getattr(self.coll, info[2]), 'active'):
                plugin = getattr(self, key)
                data[key] = plugin.as_dict()
        return data

    def set_arrays(self, arrays):
        """
        """
        # if len(arrays['positions']) == 0:
        #     return
        attributes = self.attributes
        # same length
        dnvert = len(arrays['species_index']) - \
            len(attributes['species_index'])
        if dnvert > 0:
            # self.obj.data.vertices.add(dnvert)
            self.add_vertices_bmesh(dnvert)
            # self.udpate_mesh(self.obj)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
        self.set_frames(arrays)
        for key, value in arrays.items():
            if key == "positions": continue
            self.set_attributes({key: value})
        self.update_mesh()

    @property
    def label(self):
        return self.get_label()

    @label.setter
    def label(self, label):
        self.set_label(label)

    def get_label(self):
        return self.coll.batoms.label

    def set_label(self, label):
        self.coll.batoms.label = label

    @property
    def scale(self):
        return self.get_scale()

    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)

    def get_scale(self):
        # scale = self.attributes['scale']
        scale = self.coll.batoms.scale
        return scale

    def set_scale(self, scale):
        """
        """
        self.coll.batoms.scale = scale
        scale_array = np.ones(len(self))*scale
        self.set_scale_array(scale_array)

    def set_scale_array(self, scale_array):
        self.set_attributes({'scale': scale_array})


    @property
    def species(self):
        return self.get_species()

    def get_species(self):
        return self._species

    @species.setter
    def species(self, species):
        self.set_species(species)

    def set_species(self, species):
        for key, data in species.items():
            self._species[key] = data

    @property
    def volumetric_data(self):
        return self.get_volumetric_data()

    def get_volumetric_data(self):
        return self._volumetric_data

    @volumetric_data.setter
    def volumetric_data(self, value):
        self.set_volumetric_data(value)

    def set_volumetric_data(self, value):
        for key, data in value.items():
            self._volumetric_data[key] = data

    @property
    def elements(self):
        return self.get_elements()

    def get_elements(self):
        # main elements
        arrays = self.arrays
        main_elements = self.species.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        return elements

    @property
    def model_style(self):
        return self.get_model_style()

    @model_style.setter
    def model_style(self, model_style):
        self.set_model_style(model_style)

    def get_model_style(self):
        return int(self.coll.batoms.model_style)

    def set_model_style(self, model_style):
        if int(model_style) == 0:
            self.scale = 1
        elif int(model_style) == 1:
            self.scale = 0.4
        elif int(model_style) == 2:
            self.scale = 0.4
        elif int(model_style) == 3:
            self.scale = 0.0001
        scale = np.ones(len(self))*self.scale
        self.set_attribute('scale', scale)
        self.coll.batoms.model_style = str(model_style)
        model_style_array = np.ones(len(self), dtype=int)*int(model_style)
        self.set_model_style_array(model_style_array)

    @property
    def model_style_array(self):
        return self.get_attribute('model_style')

    @model_style_array.setter
    def model_style_array(self, model_style_array):
        self.set_model_style_array(model_style_array)

    def set_model_style_array(self, model_style_array):
        model_style = {'model_style': model_style_array}
        self.set_attributes(model_style)
        self.draw()
        if self._boundary is not None:
            self.boundary.update()

    @property
    def polyhedra_style(self):
        return int(self.coll.batoms.polyhedra_style)

    @polyhedra_style.setter
    def polyhedra_style(self, polyhedra_style):
        self.coll.batoms.polyhedra_style = str(polyhedra_style)
        self.draw()

    @property
    def show_unit_cell(self):
        return self.coll.batoms.show_unit_cell

    @show_unit_cell.setter
    def show_unit_cell(self, show_unit_cell):
        if show_unit_cell:
            self.coll.batoms.show_unit_cell = show_unit_cell
            self.cell.hide = False
        else:
            self.cell.hide = True

    @property
    def radius(self):
        return self.get_radius()

    def get_radius(self):
        radius = {}
        instancers = self.species.instancers
        for sp in self.species:
            radius[sp.name] = instancers[sp.name].batoms.atom.radius
        return radius

    @property
    def radius_style(self):
        return self.coll.batoms.radius_style

    @radius_style.setter
    def radius_style(self, radius_style):
        from batoms.utils import get_default_species_data
        self.coll.batoms.radius_style = str(radius_style)
        for name, sp in self.species.items():
            data = sp.as_dict()
            props = get_default_species_data(data['elements'],
                                             radius_style=radius_style,
                                             color_style=self.color_style)
            data.update(props)
            sp.update(data)

    @property
    def color_style(self):
        return self.coll.batoms.color_style

    @color_style.setter
    def color_style(self, color_style):
        from batoms.utils import get_default_species_data
        self.coll.batoms.color_style = str(color_style)
        for name, sp in self.species.items():
            data = sp.as_dict()
            props = get_default_species_data(data['elements'],
                                             radius_style=self.radius_style,
                                             color_style=color_style)
            data.update(props)
            self.species.bpy_setting[name].color = props["elements"][sp.main_element]["color"]
            sp.update(data)

    @property
    def radii_vdw(self):
        from ase.data import vdw_radii, chemical_symbols
        # object_mode()
        radii = []
        elements = self.arrays['elements']
        for element in elements:
            if element == 'X':
                continue
            number = chemical_symbols.index(element)
            radii.append(vdw_radii[number])
        return np.array(radii)

    @property
    def size(self):
        return self.get_size()

    @size.setter
    def size(self, size):
        self.set_size(size)

    def get_size(self):
        return self.arrays['size']

    def set_size(self, size):
        scale = {}
        radius = self.radius
        for sp in self.species:
            scale[sp] = [size[sp]/radius[sp]]*3
        self.scale = scale

    def get_scaled_positions(self, cell=None):
        """
        Get array of scaled_positions.
        """
        from ase.cell import Cell
        if not cell:
            cell = self.cell
        cell = Cell.new(cell)
        scaled_positions = cell.scaled_positions(self.local_positions)
        return scaled_positions

    def get_arrays(self, batoms=None, local=False, X=False, sort=True):
        """
        """
        # object_mode()
        # tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        # radius
        radius = self.radius
        arrays.update({'radius': np.zeros(len(self))})
        for sp, value in radius.items():
            mask = np.where(arrays['species'] == sp)
            arrays['radius'][mask] = value
        # size
        arrays['size'] = arrays['radius']*arrays['scale']
        main_elements = self.species.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        arrays.update({'elements': np.array(elements, dtype='U20')})
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    @property
    def cell(self):
        return self._cell

    @cell.setter
    def cell(self, cell):
        from ase.cell import Cell
        cell = Cell.ascell(cell)
        self._cell[:] = cell

    def set_cell(self, cell, scale_atoms=False):
        """Set unit cell vectors.

        Args:
            cell (array):
                New cell array
            scale_atoms (bool, optional):
                Scale all atoms or not. Defaults to False.
        """
        from ase.cell import Cell
        cell = Cell.new(cell)
        oldcell = Cell(self.cell)
        self.cell = cell
        if scale_atoms:
            M = np.linalg.solve(oldcell.complete(), cell.complete())
            self.positions = np.dot(self.positions, M)

    @property
    def pbc(self):
        return self.get_pbc()

    @pbc.setter
    def pbc(self, pbc):
        self.set_pbc(pbc)

    def get_pbc(self):
        return list(self.cell.obj.batoms.cell.pbc)

    def set_pbc(self, pbc):
        if isinstance(pbc, bool):
            pbc = [pbc]*3
        self.cell.obj.batoms.cell.pbc = pbc

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
        index = np.zeros(nvert, dtype=int)
        layer = me.vertex_layers_int.get('index')
        layer.data.foreach_get("value", index)
        return index

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
        instancer = self.build_instancer(radius=radius,
                                         scale=scale,
                                         subdivisions=subdivisions,
                                         shape='ICO_SPHERE')
        instancer.parent = self.obj

    @property
    def shape(self):
        return self.get_shape()

    @shape.setter
    def shape(self, shape):
        self.set_shape(shape)

    def get_shape(self):
        """
        TODO
        """
        # nverts = len(self.instancer.data.vertices)
        return 'to do'

    def set_shape(self, shape):
        """
        "UV_SPHERE", "ICO_SPHERE", "CUBE"
        """
        shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]
        scale = self.scale
        if shape not in [0, 1, 2]:
            raise Exception('Shape %s is not supported!' % shape)
        scale = self.scale
        radius = self.radius
        self.clean_batom_objects(self.instancer_name)
        instancer = self.build_instancer(radius=radius,
                                         scale=scale,
                                         shape=shapes[shape])
        instancer.parent = self.obj

    def delete(self, index=[]):
        """Delete atoms by index.

        index: list
            index of atoms to be delete

        Example:

        >>> h2o.delete([1])

        Please note that index start from 0.

        """
        if isinstance(index[0], (bool, np.bool_)):
            index = np.where(index)[0]
        if isinstance(index, int):
            index = [index]
        self.delete_vertices_bmesh(index)

    def __delitem__(self, index):
        """
        """
        self.delete(index)

    def draw_constraints(self):
        """Draw label for constraints
        # TODO
        """
        pass

    def get_frames(self, local=True):
        """
        """
        frames = self.get_obj_frames(self.obj, local=local)
        return frames

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames['positions'])
        if nframe == 0:
            return
        name = self.label
        obj = self.obj
        self.set_obj_frames(name, obj, frames['positions'])

    def __getitem__(self, indices):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        from batoms.slicebatoms import SliceBatoms
        if isinstance(indices, str):
            bspecies = self.species[indices]
            return bspecies
        else:
            slicebatoms = SliceBatoms(self.label, indices, batoms=self)
            # bpy.ops.object.mode_set(mode=mode)
            return slicebatoms


    def __setitem__(self, indices, value):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[indices] = value
        self.set_positions(positions)

    def __imul__(self, m):
        """
        In-place repeat of atoms.

        >>> from batoms.batoms import Batoms
        >>> co = Batoms('co', ['C', 'O'], [[0, 0, 0], [1.2, 0, 0]])
        >>> co.cell = [3, 3, 3]
        >>> co *= [2, 2, 1]
        """
        cell = self.cell
        if isinstance(m, int):
            m = (m, m, m)
        for x, vec in zip(m, cell):
            if x != 1 and not vec.any():
                raise ValueError('Cannot repeat along undefined lattice '
                                 'vector')
        M = np.product(m)
        n = len(self)
        frames = self.frames
        nframe = len(frames)
        attributes = self.attributes
        for key, data in attributes.items():
            attributes[key] = np.tile(
                data, (M,) + (1,) * (len(data.shape) - 1))
        # positions will be updated in frames
        # repeat frames
        frames_new = []
        for i in range(0, nframe):
            positions = np.tile(
                frames[i], (M,) + (1,) * (len(frames[i].shape) - 1))
            i0 = 0
            for m0 in range(m[0]):
                for m1 in range(m[1]):
                    for m2 in range(m[2]):
                        i1 = i0 + n
                        positions[i0:i1] += np.dot((m0, m1, m2), cell)
                        i0 = i1
            frames_new.append(positions)
        attributes["positions"] = np.array(frames_new)
        self.set_arrays(attributes)
        # self.set_frames(frames_new)
        self.cell.repeat(m)
        self.update_gn_cell()
        if self.volumetric_data is not None:
            self._volumetric_data *= m
        self.species.update_geometry_node()
        if self._boundary is not None:
            self.boundary.update()
        self.draw()
        return self

    def repeat(self, m):
        """
        """
        self *= m
        return self

    def __mul__(self, m):
        self.repeat(m)
        return self

    def copy(self, label, displacement=[2, 2, 2]):
        """
        Return a copy.

        label: str
            The label of the copy.

        For example, copy H species:

        >>> h_new = h.copy(label = 'h_new', species = 'H')
        # TODO Support copy of other prroperties: materials, info ...
        """
        # object_mode()
        # copy object first
        arrays = self.arrays
        batoms = self.__class__(label=label,
                                location=self.location + np.array(displacement),
                                species=arrays['species'],
                                positions=self.local_positions,
                                pbc=self.pbc,
                                cell=self.cell.array,
                                )
        return batoms

    def extend(self, other):
        """
        Extend batom object by appending batoms from *other*.

        >>> slab = au111 + co
        todo: merge bond setting
        """
        # could also use self.add_arrays(other.positions)
        # object_mode()
        # merge bondetting
        self.bond.settings.extend(other.bond.settings)
        # creat new selections
        n1 = len(self)
        n2 = len(other)
        indices1 = list(range(n1))
        indices2 = list(range(n1, n1 + n2))
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        other.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.join()
        # update species and species_index
        self._species.extend(other._species)
        self.selects.add(self.label, indices1)
        self.selects.add(other.label, indices2)
        # remove shape key from mol
        sp = self.obj.data.shape_keys.key_blocks.get('Basis_%s' % other.label)
        self.obj.shape_key_remove(sp)
        # remove old
        bpy.ops.batoms.delete(label=other.label)

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

    def separate(self):
        """
        Separate batoms object based on selects.

        >>> slab = au111 + co
        >>> slab.separate()
        >>> au111=Batoms('au111')
        >>> co=Batoms('co')
        """
        # could also use self.add_arrays(other.positions)
        # object_mode()
        arrays = self.arrays
        selects = self.selects.selects
        self_indices = []
        for name, sel in selects.items():
            if name == "all": continue
            if name == self.label:
                self_indices = sel.indices
                continue
            indices = sel.indices
            sel = self.__class__(name, species = arrays['species'][indices],
                                positions = arrays['positions'][indices])
            new_arrays = {}
            for key, array in arrays.items():
                if key in ["positions", "species"]: continue
                new_arrays[key] = array[indices]
            sel.set_attributes(new_arrays)
        if len(self_indices) == 0:
            # remove old
            bpy.ops.batoms.delete(label=self.label)
            bpy.context.view_layer.objects.active = sel.obj
        else:
            remove_indices = list(set(range(len(self))) - set(self_indices))
            del self[remove_indices]
            bpy.context.view_layer.objects.active = self.obj




    def __iter__(self):
        batom = self.obj
        for i in range(len(self)):
            yield batom.matrix_world @ batom.data.vertices[i].co

    def __repr__(self) -> str:
        text = []
        text.append('label={0}, '.format(self.label))
        text.append('species=%s, ' % (list(self.species.species)))
        text.append('cell={0}, '.format(self.cell))
        text.append('pbc={0}'.format(self.pbc))
        # text.append('positions={0}'.format(self.positions))
        text = "".join(text)
        text = "Batoms(%s)" % text
        return text

    def replace(self, indices, species):
        """Replace species.
        Parameters:

        indices: list
            indices of atoms will be replaced.
        species: str
            atoms will be changed to this species.

        >>> h2o.replace([0], 'O_1')
        >>> h2o.replace([0], 'S')
        >>> h2o.replace(range(2), 'N')

        # TODO remove species which is completely replaced.
        """
        from batoms.utils import get_default_species_data
        # if kind exists, merger, otherwise build a new kind and add.
        mode = self.obj.mode
        bpy.context.view_layer.objects.active = self.obj
        if isinstance(species, str):
            ele = species.split('_')[0]
            species = [species, {'elements': {ele: {"occupancy": 1.0}}}]
            props = get_default_species_data(species[1]['elements'],
                                             radius_style=self.radius_style,
                                             color_style=self.color_style)
            species[1].update(props)
            species[1]["elements"].update(props["elements"])
        if species[0] not in self.species:
            if len(species[0]) > 4:
                print("Warning: the name of the species: {}, should not be longer than four characters.".format(species[0]))
                logger.warning("The name of the species: {}, should not be longer than four characters.".format(species[0]))
            self.species[species[0]] = species[1]
            # add geometry node
            self.add_geometry_node(
                species[0], self.species.instancers[species[0]])
            #
        species_index = self.attributes['species_index']
        species_array = self.attributes['species']
        #
        species_index[indices] = string2Number(species[0])
        species_array[indices] = species[0]
        self.set_attributes({'species_index': species_index})
        self.set_attributes({'species': species_array})
        bpy.context.view_layer.objects.active = self.obj
        # print(mode)
        # print(self.species)
        # for sp in self.species:
        self.bond.settings.add_species(species[0])
        self.polyhedra.settings.add_species(species[0])

    def auto_build_species(self, tol = 1e-5):
        """Auto build species based on equivalent_atoms.

        Args:
            tol (float, optional): _description_. Defaults to 1e-5.
        """
        from batoms.data import jmol_colors
        from batoms.utils import get_equivalent_atoms
        equivalent_atoms = get_equivalent_atoms(self.as_ase(), tol)
        species = self.species
        species_array = self.get_attribute('species')
        icolor = 1
        for name, sp in species.items():
            n = 0
            atoms_sp = equivalent_atoms[sp.indices]
            unique_species = np.unique(atoms_sp)
            for usp in unique_species:
                indices = np.where(equivalent_atoms==usp)[0]
                if len(indices) > 0:
                    if n > 0:
                        name = '{}_{}'.format(species_array[indices[0]], n)
                        self.replace(indices, name)
                        # self.species[name].color = jmol_colors[icolor]
                    n += 1
                    icolor += 1
            if n > 10:
                logger.critical("Name of species longer than 4 characters.")
            logger.info("{} kinds of species are found for {}".format(n, name))

    def add_atoms(self, arrays):
        """ Used to add small number of atoms
        Todo: find a fast way.
        """
        from batoms.utils import local2global
        import bmesh
        # object_mode()
        n0 = len(self)
        positions = arrays.pop('positions')
        n1 = n0 + len(positions)
        # add positions
        positions = local2global(positions,
                                 np.array(self.obj.matrix_world),
                                 reversed=True)
        bm = bmesh.new()
        bm.from_mesh(self.obj.data)
        bm.verts.ensure_lookup_table()
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.obj.data)
        bm.clear()
        # add species
        self.species.add(list(set(arrays['species'])))
        self.set_attribute_with_indices(
            'species', range(n0, n1), arrays['species'])
        if 'species_index' not in arrays:
            species_index = [string2Number(sp) for sp in arrays['species']]
            self.set_attribute_with_indices('species_index',
                                            range(n0, n1),
                                            species_index)
        # add arrays
        if 'show' not in arrays:
            show = np.ones(n1 - n0)
            self.set_attribute_with_indices('show', range(n0, n1), show)
        if 'scale' not in arrays:
            scale = np.ones(n1 - n0)
            self.set_attribute_with_indices('scale', range(n0, n1), scale)
        # add bond setting
        bond_dicts = self.bond.settings.as_dict()
        species0 = self.species.keys()
        species1 = np.unique(arrays['species'])
        for sp1 in species1:
            for sp0 in species0:
                pair1 = (sp1, sp0)
                pair2 = (sp0, sp1)
                if pair1 not in bond_dicts and pair2 not in bond_dicts:
                    self.bond.settings.add(pair1)
                    sp3 = self.bond.settings.find(pair1)
                    sp4 = self.bond.settings.find(pair2)
                    if sp3:
                        self.bond.update_geometry_node_instancer(
                            sp3.as_dict())
                    if sp4:
                        self.bond.update_geometry_node_instancer(
                            sp4.as_dict())

    def get_cell(self):
        if self.label not in bpy.data.collections:
            return None
        bcell = bpy.data.collections['%s_cell' % self.label]
        cell = np.array(
            [bcell.matrix_world @ bcell.data.vertices[i].co for i in range(3)])
        return cell

    def get_distances(self, i, indices, mic=False):
        """
        Return distances of atom No.i with a list of atoms.

        Use mic=True to use the Minimum Image Convention.

        >>> h2o.get_distances(0, 1)
        >>> h2o.get_distances(0, [1, 2])
        """
        from ase.geometry import get_distances
        positions = self.positions
        p1 = positions[i]
        p2 = positions[indices]
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)
        D_len.shape = (-1,)
        return D_len

    def get_angle(self, i1, i2, i3, mic=False):
        """
        Get angle in degrees between the vectors i2->i1 and
        i2->i3.
        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.

        >>> h2o.get_angle(0, 1, 2)

        """
        from ase.geometry import get_angles
        positions = self.positions
        p1 = positions[i1]
        p2 = positions[i2]
        p3 = positions[i3]
        v12 = p1 - p2
        v32 = p3 - p2
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        return get_angles([v12], [v32], cell=cell, pbc=pbc)

    def get_dihedral(self, i1, i2, i3, i4, mic=False):
        """

        >>> h2o2.get_dihedral(0, 1, 2, 3)

        """
        from ase.geometry import get_dihedrals
        positions = self.positions
        p1 = positions[i1]
        p2 = positions[i2]
        p3 = positions[i3]
        p4 = positions[i4]
        v21 = p2 - p1
        v32 = p3 - p2
        v43 = p4 - p3
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        return get_dihedrals([v21], [v32], [v43], cell=cell, pbc=pbc)

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned.
        """
        return self.as_ase().get_center_of_mass(scaled=scaled)

    def get_center_of_geometry(self, coll=None):
        """
        """
        vertices = self.get_all_vertices(coll=coll, cell=self.show_unit_cell)
        canvas = np.zeros([2, 3])
        canvas[0] = np.min(vertices, axis=0)
        canvas[1] = np.max(vertices, axis=0)
        center = np.mean(canvas, axis=0)
        return center

    @property
    def show(self):
        return self.get_show()

    @show.setter
    def show(self, state):
        self.set_show(state)

    def get_show(self):
        return self.attributes['show']

    def set_show(self, show, only_atoms=False):
        #
        if not only_atoms:
            for obj in self.obj.children:
                obj.hide_render = not show
                obj.hide_set(not show)
        #
        if isinstance(show, (bool, int)):
            show = np.ones(len(self), dtype=bool)*show
        self.set_attributes({'show': show})
        #TODO: /home/xing/apps/beautiful-atoms/batoms/batoms.py:1269:
        # DeprecationWarning: In future, it will be an error for 'np.bool_'
        #  scalars to be interpreted as an index
        self.coll.batoms.show = show[0]

    @property
    def wrap(self):
        return self.get_wrap()

    @wrap.setter
    def wrap(self, state):
        self.set_wrap(state)

    def get_wrap(self):
        return list(self.coll.batoms.wrap)

    def set_wrap(self, wrap):
        #
        nodes = self.gnodes.node_group.nodes
        self.update_gn_cell()
        if isinstance(wrap, bool):
            wrap = [wrap]*3
        self.coll.batoms.wrap = list(wrap)
        wrap = np.array(wrap)
        if not wrap.any():
            # switch off
            n = len(nodes[
                '%s_CombineXYZ_' % self.label].outputs['Vector'].links)
            if n > 0:
                link = nodes['%s_CombineXYZ_' %
                             self.label].outputs['Vector'].links[0]
                self.gnodes.node_group.links.remove(link)
        else:
            self.gnodes.node_group.links.new(
                nodes['%s_CombineXYZ_' % self.label].outputs[0],
                nodes['%s_SetPosition' % self.label].inputs['Position'])
        self.gnodes.node_group.update_tag()
        # TODO: support selection

    def update_gn_cell(self):
        # update cell
        cell = self.cell.array
        if np.isclose(np.linalg.det(self.cell.array), 0):
            cell = np.eye(3)
        icell = np.linalg.inv(cell)
        # set positions
        gn = self.gnodes
        for i in range(3):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_VectorDot%s_%s' % (self.label, i, ''),
                                    'ShaderNodeVectorMath')
            tmp.operation = 'DOT_PRODUCT'
            tmp.inputs[1].default_value = cell[:, i]
            # icell
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_VectorDot%s_%s' % (
                                        self.label, i, 'scaled'),
                                    'ShaderNodeVectorMath')
            tmp.operation = 'DOT_PRODUCT'
            tmp.inputs[1].default_value = icell[:, i]

    def get_spacegroup_number(self, symprec=1e-5):
        """
        """
        try:
            import spglib
        except ImportError:
            return 1
        sg = spglib.get_spacegroup((self.cell, self.get_scaled_positions(),
                                    self.arrays['numbers']),
                                   symprec=symprec)
        if sg is None:
            return None
        no = int(sg[sg.find('(') + 1:sg.find(')')])
        return no

    def find_primitive(self, ):
        """
        """
        from ase import Atoms
        import spglib
        atoms = self.atoms
        lattice = atoms.cell
        points = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, points, numbers)
        lattice, points,  numbers = spglib.find_primitive(cell)
        atoms = Atoms(numbers=numbers, scaled_positions=points, cell=lattice)
        return atoms

    def get_all_vertices(self, coll=None, cell=True):
        """
        Get position of all vertices from all mesh in batoms.
        Used for plane boundary and calc_camera_data
        """
        positions = self.positions
        # isosurface, plane
        if coll is None:
            coll = self.coll
        for obj in coll.all_objects:
            if obj.type != 'MESH':
                continue
            if obj.batoms.type == 'INSTANCER' or obj.batoms.type == 'VOLUME':
                continue
            if obj.hide_get():
                continue
            n = len(obj.data.vertices)
            vertices = np.empty(n*3, dtype=np.float64)
            obj.data.vertices.foreach_get('co', vertices)
            vertices = vertices.reshape((n, 3))
            vertices = np.append(vertices, np.ones((n, 1)), axis=1)
            mat = np.array(obj.matrix_world)
            vertices = mat.dot(vertices.T).T
            # (natom, 4) back to (natom, 3)
            vertices = vertices[:, :3]
            positions = np.concatenate((positions, vertices), axis=0)
        return positions

    def get_canvas_box(self, direction=[0, 0, 1], padding=None, coll=None):
        """
        Calculate the canvas box from [0, 0, 1] and other direction.

        """
        from batoms.utils import get_canvas
        vertices = self.get_all_vertices(coll=coll, cell=self.show_unit_cell)
        canvas = get_canvas(vertices, direction=direction, padding=padding)
        width = canvas[1][0] - canvas[0][0]
        height = canvas[1][1] - canvas[0][1]
        depth = canvas[1][2] - canvas[0][2]
        return width, height, depth

    def lock_to_camera(self, obj):
        from batoms.utils.butils import lock_to
        for sp, instancer in self.species.instancers.items():
            lock_to(instancer, obj, location=False, rotation=True)

    @property
    def bonds(self):
        """bonds object."""
        deprecated('"bonds" will be deprecated in the furture, please use "bond".')
        return self.bond

    @property
    def bond(self):
        """bond object."""
        from batoms.bond.bond import Bond, default_bond_datas
        if self._bond is not None:
            return self._bond
        bond = Bond(self.label, batoms=self)
        self.bond = bond
        return bond

    @bond.setter
    def bond(self, bond):
        self._bond = bond

    @property
    def polyhedras(self):
        """polyhedras object."""
        deprecated('"polyhedras" will be deprecated in the furture, please use "polyhedra".')
        return self.polyhedra

    @property
    def polyhedra(self):
        """polyhedra object."""
        from batoms.polyhedra.polyhedra import Polyhedra, default_polyhedra_datas
        if self._polyhedra is not None:
            return self._polyhedra
        polyhedra = Polyhedra(self.label,
                                polyhedra_datas=default_polyhedra_datas,
                                batoms=self)
        self.polyhedra = polyhedra
        return polyhedra

    @polyhedra.setter
    def polyhedra(self, polyhedra):
        self._polyhedra = polyhedra

    @property
    def boundary(self):
        """boundary object."""
        from batoms.boundary import Boundary, default_boundary_datas
        if self._boundary is not None:
            return self._boundary
        boundary = Boundary(self.label, batoms=self,
                            location=self.location)
        self._boundary = boundary
        return boundary

    @boundary.setter
    def boundary(self, boundary):
        if isinstance(boundary, (int, float)):
            boundary = np.array([[-boundary, 1 + boundary]]*3)
        elif len(boundary) == 3:
            if isinstance(boundary[0], (int, float)):
                boundary = np.array([[-boundary[0], 1 + boundary[0]],
                                    [-boundary[1], 1 + boundary[1]],
                                    [-boundary[2], 1 + boundary[2]]])
            elif len(boundary[0]) == 2:
                boundary = np.array(boundary)
        self.update_gn_cell()
        self.boundary[:] = boundary

    @property
    def isosurfaces(self):
        """isosurfaces object."""
        deprecated('"isosurfaces" will be deprecated in the furture, please use "isosurface".')
        return self.isosurface



    @property
    def ms(self):
        """ms object."""
        deprecated('"ms" will be deprecated in the furture, please use "molecular_surface".')
        return self.molecular_surface


    @property
    def render(self):
        """Render object."""
        from batoms.render.render import Render
        if self._render is not None:
            return self._render
        render = Render()
        self.render = render
        return render

    @render.setter
    def render(self, render):
        render.batoms = self
        self._render = render
        self.lock_to_camera(render.camera.obj)

    def get_image(self, viewport=None, engine=None,
                  frame=1,
                  animation=False,
                  output=None,
                  center=None,
                  padding=None,
                  canvas=None,
                  gpu=False):
        """Rendering the model.

        Ask the attached render to rendering the model.

        frame: int
        animation: bool
        output: str
        center: array
        padding: float
        canvas: array of 3

        """
        if output is None:
            output = '%s.png' % self.label
        if self.render is None:
            raise RuntimeError('Batoms object has no render.')
        self.render.run_render = True
        if viewport is not None:
            self.render.viewport = viewport
        if engine is not None:
            self.render.engine = engine
        self.render.gpu = gpu
        image = self.render.run(self, frame=frame,
                                animation=animation,
                                output=output,
                                center=center,
                                padding=padding,
                                canvas=canvas,
                                )
        return image

    def draw(self, model_style=None):
        """
        Draw atoms, bond, polyhedra, .

        Parameters:

        model_style: str
        """
        if model_style is not None:
            self.model_style = model_style
        # self.draw_cell()
        self.draw_space_filling()
        self.draw_ball_and_stick()
        self.draw_polyhedra()
        self.draw_wireframe()

    def draw_space_filling(self):
        # mask = np.where(self.model_style_array == 0, True, False)
        # self.set_attribute_with_indices('scale', mask, self.scale)
        self.boundary.update()


    def draw_ball_and_stick(self):
        mask = np.where(self.model_style_array >= 1, True, False)
        if not mask.any():
            from batoms.bond.bond import default_bond_datas
            self.bond.set_arrays(default_bond_datas.copy())
            return
        # self.set_attribute_with_indices('scale', mask, self.scale)
        self.bond.hide = False
        self.boundary.hide = False
        self.bond.update()

    def draw_polyhedra(self):
        mask = np.where(self.model_style_array == 2, True, False)
        if not mask.any():
            from batoms.polyhedra.polyhedra import default_polyhedra_datas
            self.polyhedra.set_arrays(default_polyhedra_datas)
            return
        self.polyhedra.update()
        self.set_attribute_with_indices('show', mask, True)
        if self.polyhedra_style == 0:
            # self.set_attribute_with_indices('scale', mask, self.scale)
            self.boundary.hide = False
            self.bond.hide = False
            self.bond.search_bond.hide = False
            # self.bond.update()
        if self.polyhedra_style == 1:
            # only hide bond
            # self.set_attribute_with_indices('scale', mask, self.scale)
            self.bond.hide = True
            # show atoms
            if self.bond._search_bond is not None:
                self.bond.search_bond.hide = False
            self.boundary.hide = False
        elif self.polyhedra_style == 2:
            # hide bond and atoms, except center atoms
            for b in self.bond.settings:
                if b.polyhedra:
                    mask1 = np.where(
                        self.attributes['species'] == b.species1, True, False)
                    self.set_attribute_with_indices('scale', mask1, self.scale)
                    mask[mask1] = False
            scale = 0
            # self.set_attribute_with_indices('scale', mask, self.scale)
            # self.boundary.hide = True
            self.bond.hide = True
            self.bond.search_bond.hide = True
            # self.set_attribute_with_indices('show', mask, False)
        elif self.polyhedra_style == 3:
            # hide all bond and atoms
            scale = 0.0001
            self.set_attribute_with_indices('scale', mask, scale)
            self.boundary.hide = True
            self.bond.hide = True
            self.bond.search_bond.hide = True
            # self.set_attribute_with_indices('show', mask, False)

    def draw_wireframe(self):
        mask = np.where(self.model_style_array == 3, True, False)
        # self.set_attribute_with_indices('show', mask, 0)
        self.set_attribute_with_indices('scale', mask, self.scale)
        # self.update(mask)

    def as_ase(self, local=True, with_attribute=True):
        """
        local: bool
            if True, use the origin of uint cell as the origin
        """
        from ase import Atoms
        frames = self.get_frames()
        arrays = self.arrays
        nframe = len(frames)
        images = []
        for f in range(nframe):
            positions = frames[f]
            if not local:
                positions += self.cell.origin
            atoms = Atoms(symbols=arrays['elements'],
                          positions=positions,
                          cell=self.cell, pbc=self.pbc)
            if with_attribute:
                for name, array in arrays.items():
                    if name in ['elements', 'positions', 'numbers']:
                        continue
                    atoms.set_array(name, np.array(array))
            images.append(atoms)
        if nframe == 1:
            return images[0]
        else:
            return images

    def as_pybel(self, export_bond=False):
        """Convert an Batoms object to an OBMol object.

        Returns:
            OBMOL: OBMOL
        """
        from openbabel import pybel
        from ase.data import atomic_numbers
        mol = pybel.ob.OBMol()
        arrays = self.arrays
        natom = len(self)
        for i in range(natom):
            a = mol.NewAtom()
            a.SetAtomicNum(atomic_numbers[arrays['elements'][i]])
            a.SetVector(arrays['positions'][i][0],
                        arrays['positions'][i][1],
                        arrays['positions'][i][2])
        if export_bond:
            bond_arrays = self.bond.arrays
            nbond = len(bond_arrays)
            for i in range(nbond):
                mol.AddBond(int(bond_arrays['atoms_index1'][i]) + 1,
                            int(bond_arrays['atoms_index2'][i]) + 1,
                            bond_arrays['order'][i])
        else:
            mol.ConnectTheDots()
            mol.PerceiveBondOrders()
        mol = pybel.Molecule(mol)
        return mol

    def write(self, filename, local=True):
        """
        Save batoms to structure file.
        >>> h2o.write('h2o.xyz')
        """
        self.as_ase(local).write(filename)

    def transform(self, matrix=None):
        """
        Transformation matrix
        make sure np.linalg.det(P) > 0, otherwise,
        a error raised by ASE last column [5, 5, 0] to
        translate new atoms, thus not overlap with old one.

        """
        from ase.build.supercells import make_supercell
        if matrix is not None:
            rotation = np.array([matrix[0][:3], matrix[1][:3], matrix[2][:3]])
            translation = np.array([matrix[0][3], matrix[1][3], matrix[2][3]])
            atoms = self.as_ase()
            atoms = make_supercell(atoms, rotation)
            batoms = self.__class__(label='%s_transform' % self.label,
                                    from_ase=atoms,
                                    model_style=self.model_style)
        else:
            return
        batoms.translate(translation)
        # self.hide = True
        return batoms

    @property
    def draw_text(self):
        """draw_text object.
        Shoud be a global variable.
        """
        from batoms.draw.draw_screen import DrawText
        dns = bpy.app.driver_namespace
        if self.label in dns:
            return dns[self.label]
        context = bpy.context
        dns[self.label] = DrawText(context, self.label)
        return dns[self.label]

    @property
    def show_label(self):
        return self.coll.batoms.show_label

    @show_label.setter
    def show_label(self, label="species"):
        self.coll.batoms.show_label = str(label)
        arrays = self.arrays
        positions = arrays["positions"]
        self.draw_text.remove_handle()
        if label is None or label=="":
            return
        if label.lower() == "index":
            n = len(positions)
            texts = [str(i) for i in range(n)]
        elif label.lower() in arrays:
            texts = arrays[label]
        else:
            raise Exception("%s does not have %s attribute" %
                            (self.label, label))
        self.draw_text.add_handle(positions, texts)

    def optmize(self, forcefield="mmff94", steps=500):
        """_summary_

        Args:
            forcefield (str, optional): _description_. Defaults to "mmff94".
            steps (int, optional): _description_. Defaults to 500.
        """
        from openbabel import openbabel as ob
        from openbabel import pybel
        # from batoms.utils import read_from_pybel
        mol = self.as_pybel()
        mol.localopt(forcefield, steps)
        positions = []
        for atom in ob.OBMolAtomIter(mol.OBMol):
            positions.append([atom.GetX(), atom.GetY(), atom.GetZ()])
        # species, positions, arrays, cell, pbc, info = read_from_pybel(mol)
        self.positions = positions

    def export_mesh(self, filename, with_cell=False,
                    with_polyhedra=False,
                    with_bond=False,
                    with_boundary=False,
                    with_search_bond=False):
        """_summary_
        ctmviewer

        Args:
            filename (_type_): _description_

        """
        # remove shape key
        self.realize_instances = True
        objs = [self.obj]
        if with_cell:
            self.cell.obj.select_set(True)
        if with_bond:
            self.bond.realize_instances = True
            objs.insert(0, self.bond.obj)
        if with_boundary:
            self.boundary.realize_instances = True
            objs.insert(0, self.boundary.obj)
        if with_search_bond:
            self.bond.search_bond.realize_instances = True
            objs.insert(0, self.bond.search_bond.obj)
        if with_polyhedra:
            self.polyhedra.realize_instances = True
            objs.insert(0, self.polyhedra.obj)
        for obj in objs:
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.shape_key_remove(all=True)
            bpy.ops.object.modifier_apply(modifier=obj.modifiers[0].name)
            obj.select_set(True)
        if filename.endswith('.x3d'):
            bpy.ops.export_scene.x3d(filepath=filename, use_selection=True)
        elif filename.endswith('.obj') or filename.endswith('.mtl'):
            bpy.ops.export_scene.obj(filepath=filename, use_selection=True)
        elif filename.endswith('.fbx'):
            bpy.ops.export_scene.fbx(filepath=filename, use_selection=True)
        else:
            raise('File format %s is not supported.' % filename)

    @property
    def crystal_view(self):
        return self.get_crystal_view()

    @crystal_view.setter
    def crystal_view(self, state):
        self.set_crystal_view(state)

    def get_crystal_view(self):
        return self.coll.batoms.crystal_view

    def set_crystal_view(self, crystal_view):
        #
        self.coll.batoms.crystal_view = crystal_view
        if crystal_view:
            self.boundary = 0.01
            self.bond.show_search = True
            self.bond.update()
            self.polyhedra.update()
        else:
            self.boundary = 0
            self.bond.show_search = False
            self.bond.update()
            self.polyhedra.update()

    @property
    def segments(self):
        return self.get_segments()

    @segments.setter
    def segments(self, segments):
        self.set_segments(segments)

    def get_segments(self):
        return self.coll.batoms.segments[:]

    def set_segments(self, segments):
        if isinstance(segments, int):
            segments = [segments, segments]
            # raise Exception('Segments should be int!')
        species = self._species
        self.coll.batoms.segments = segments
        for name, sp in species.items():
            sp.data.segments = segments
            sp.update(sp.data.as_dict())

    def auto_assign_charge(self, force_field="mmff94"):
        """Estimates atomic partial charges using openbabel.

        MMFF94 Partial Charges are assigned using a four stage algorithm.

        Args:
            force_field (str, optional): force field. Defaults to "mmff94".
        """
        tstart = time()
        ob = self.as_pybel()
        charges = np.array(ob.calccharges(model = force_field))
        self.set_attributes({'charges':charges})
        logger.debug("time: {:1.3f}".format(time() - tstart))
        return charges

    def calc_electrostatic_potential(self, points):
        """Calc electrostatic potential on given points
            based on partial charges.

        Args:
            points (array): points

        Returns:
            array: potentials at given points
        """
        from ase.geometry import get_distances
        tstart = time()
        positions = self.positions
        charges = self.get_attribute('charges')
        # parameter
        eps = 8.854e-12 # Permittivity of free space (F.m-1)
        k = 1/4/np.pi/eps
        # calculate distance
        D, D_len = get_distances(positions, points)
        # matrix
        pot_mat = k/D_len*charges[:, None]
        # sum over all charge (atoms)
        potentials = np.sum(pot_mat, axis = 0)
        logger.debug("time: {:1.3f}".format(time() - tstart))
        return potentials

    def set_hide(self, state):
        self.obj.hide_render = state
        self.obj.hide_set(state)
        self.cell.hide = True
        self.bond.hide = True


def hook_plugins(cls, plugin_info):
    tstart = time()
    def create_getter(name, data):
        def plugin_getter(self):
            """plugin object."""
            import importlib
            _name = '_{}'.format(name)
            if hasattr(self, _name):
                plugin = getattr(self, _name)
                logger.debug('has plugin: {}'.format(name))
                return plugin
            else:
                logger.debug('Does not have plugin: {}'.format(name))
                plugin_module = getattr(importlib.import_module("batoms.plugins.{}".format(data[0])), data[1])
                plugin = plugin_module(self.label, batoms=self)
                setattr(self, _name, plugin)
                return plugin
        return plugin_getter

    def create_setter(name):
        def plugin_setter(self, plugin):
            """Get plugin
            Args:
                plugin (_type_): _description_
            """
            logger.debug('set plugin: {}'.format(name))
            _name = '_{}'.format(name)
            setattr(self, _name, plugin)
        return plugin_setter

    for name, data in plugin_info.items():
        logger.debug("Add plugin: {}".format(name))
        setattr(cls, name, property(
                    fget = create_getter(name, data),
                    fset = create_setter(name),
                    doc = 'test',
                ))

    logger.debug("Hook plugins: {:1.2f}".format(time() - tstart))



hook_plugins(Batoms, plugin_info)
