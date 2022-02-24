"""Definition of the Batoms class.

This module defines the Batoms object in the batoms package.

# TODO: get evaluated positiosn
# TODO: add location for boundary, bonds and all child objects
# TODO: add feature: cavity
"""
import bpy
from batoms.batom import Batom
from batoms.bond import Bonds, default_bond_datas
from batoms.polyhedra import Polyhedras, default_polyhedra_datas
from batoms.boundary import Boundary, default_boundary_datas
from batoms.bspecies import Bspecies
from batoms.cell import Bcell
from batoms.render.render import Render
from batoms.bselect import Selects
from batoms.utils import string2Number, read_from_ase, read_from_pymatgen
from batoms.base.collection import BaseCollection
from batoms.base.object import ObjectGN
from batoms.isosurfacesetting import IsosurfaceSetting
from batoms.planesetting import PlaneSetting
from batoms.mssetting import MSsetting
from batoms.ribbon.ribbon import Ribbon
from batoms.utils.butils import object_mode, show_index, \
    get_nodes_by_name, compareNodeType
import numpy as np
# from time import time

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]

default_attributes = [
    ['select', 'INT'],
    ['species_index', 'INT'],
    ['species', 'STRING'],
    ['show', 'BOOLEAN'],
    ['scale', 'FLOAT'],
    ['model_style', 'INT'],
]

default_GroupInput = [
    ['select', 'NodeSocketInt'],
    ['species_index', 'NodeSocketInt'],
    ['show', 'NodeSocketBool'],
    ['scale', 'NodeSocketFloat'],
]

subcollections = ['instancer', 'surface', 'ribbon', 'plane']


class Batoms(BaseCollection, ObjectGN):
    def __init__(self,
                 label='batoms',
                 species=[],
                 positions=[],
                 attributes={},
                 species_props={},
                 info={},
                 pbc=False, cell=None,
                 location=np.array([0, 0, 0]),
                 show_unit_cell=True,
                 volume=None,
                 scale=1.0,
                 model_style=0,
                 polyhedra_style=0,
                 from_ase=None,
                 from_pymatgen=None,
                 metaball=False,
                 movie=True,
                 segments=None,
                 ):
        """Batom Class
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
            polyhedra_style (int, optional):
                Enum in [0, 1, 2]. Defaults to 0.
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
        if from_ase is not None:
            species, positions, attributes, cell, pbc, info = \
                read_from_ase(from_ase)
        if from_pymatgen is not None:
            species, positions, attributes, cell, pbc, info = \
                read_from_pymatgen(from_pymatgen)
        if len(species) == 0 and self.check_batoms(label):
            self.from_batoms(label)
        else:
            self.set_collection(label)
            self._cell = Bcell(label, cell, batoms=self)
            positions = np.array(positions)
            if len(positions.shape) == 3:
                self._frames = positions
                positions = self._frames[0]
            else:
                self._frames = np.array([positions])
            #
            natom = len(positions)
            self.build_object(label, positions, location)
            self.selects = Selects(label, self)
            if not species_props:
                species_props = {sp: {'elements': {sp.split('_')[0]: 1.0}}
                                 for sp in species}
            self.build_geometry_node()
            self._species = Bspecies(
                label, label, species_props, self, segments=segments)
            self.selects.add('all', np.arange(len(self)))
            if isinstance(scale, (int, float)):
                scale = np.ones(natom)*scale
            show = np.ones(natom, dtype=int)
            species_index = [string2Number(sp) for sp in species]
            attributes.update({'species': species,
                               'species_index': species_index,
                               'scale': scale,
                               'show': show,
                               })
            self.set_attributes(attributes)
            if volume is not None:
                self.build_volume(volume)
            self.set_pbc(pbc)
            # self.label = label
            if movie:
                self.set_frames()
        self.isosurfacesetting = IsosurfaceSetting(self.label, batoms=self)
        self.planesetting = PlaneSetting(self.label, batoms=self)
        self.mssetting = MSsetting(self.label, probe=1.4, batoms=self)
        self.ribbon = Ribbon(self.label, batoms=self, datas=info, update=True)
        self._render = None
        self._bonds = None
        self._polyhedras = None
        self._boundary = None
        show_index()
        self.hideOneLevel()

    def set_collection(self, label):
        """Build main collection and its child collections.

        Args:
            label (str): name of Batoms object

        Raises:
            Exception: The label is already in use.
        """
        if bpy.data.collections.get(label):
            raise Exception("Failed, the name %s already in use!" % label)
        coll = bpy.data.collections.new(label)
        bpy.data.scenes['Scene'].collection.children.link(coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s' % (label, sub_name))
            coll.children.link(subcoll)
        coll.batoms.flag = True
        coll.batoms.label = label

    def hideOneLevel(self):
        """Hide one level of collecitons in the outline in Blender
        """
        from batoms.utils.butils import hideOneLevel
        hideOneLevel()

    def build_object(self, label, positions, location=[0, 0, 0]):
        """Build the main Batoms object

        Args:
            label (str):
                Name of the object
            positions (array):
                Positions of each atoms
            location (list, optional):
                Location of the object. Defaults to [0, 0, 0].
        """
        self.delete_obj(label)
        mesh = bpy.data.meshes.new(label)
        # Add attributes
        for attribute in default_attributes:
            mesh.attributes.new(
                name=attribute[0], type=attribute[1], domain='POINT')
        obj = bpy.data.objects.new(label, mesh)
        obj.data.from_pydata(positions, [], [])
        obj.location = location
        obj.batoms.flag = True
        obj.batoms.label = label
        self.coll.objects.link(obj)
        # add cell object as its child
        self.cell.obj.parent = self.obj

    def build_geometry_node(self):
        """Geometry node for instancing sphere on vertices!
        """
        name = 'GeometryNodes_%s' % self.label
        modifier = self.obj.modifiers.new(name=name, type='NODES')
        modifier.node_group.name = name
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes.get('Group Input')
        GroupOutput = modifier.node_group.nodes.get('Group Output')
        # add new output sockets
        for att in default_GroupInput:
            GroupInput.outputs.new(type=att[1], name=att[0])
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
            gn, PositionBatoms, icell, 'scaled')
        VectorWrap = get_nodes_by_name(gn.node_group.nodes,
                                       '%s_VectorWrap' % (self.label),
                                       'ShaderNodeVectorMath')
        VectorWrap.operation = 'WRAP'
        VectorWrap.inputs[1].default_value = [1, 1, 1]
        VectorWrap.inputs[2].default_value = [0, 0, 0]
        gn.node_group.links.new(
            scaledPositionsNode.outputs[0], VectorWrap.inputs['Vector'])
        PositionsNode = self.vectorDotMatrix(gn, VectorWrap, cell, '')
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
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
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
                                JoinGeometry.inputs['Geometry'])

    def build_volume(self, volume):
        """Save volumetric data as a mesh

        Args:
            volume (array):
                volumetric data, e.g. electron density
        """
        # remove old volume point
        # tstart = time()
        if volume is None:
            return
        name = "%s_volume" % self.label
        if name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[name], do_unlink=True)
        shape = volume.shape
        volume = volume.reshape(-1, 1)
        npoint = len(volume)
        dn = 3 - npoint % 3
        verts = np.append(volume, np.zeros((dn, 1)), axis=0)
        verts = verts.reshape(-1, 3)
        mesh = bpy.data.meshes.new("%s_volume" % self.label)
        mesh.from_pydata(verts, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        obj.data = mesh
        obj.batoms.bvolume.is_bvolume = True
        obj.batoms.bvolume.shape = shape
        self.coll.objects.link(obj)
        obj.hide_set(True)
        obj.hide_render = True
        # print('Draw volume: {0:1.2f}'.format(time() - tstart))

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
        elif not bpy.data.collections[label].batoms.flag:
            flag = False
        if label not in bpy.data.objects:
            flag = False
        elif not bpy.data.objects[label].batoms.flag:
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
        self._cell = Bcell(label=label)
        self._species = Bspecies(label, label, {}, self)
        self.selects = Selects(label, self)

    @property
    def volume(self):
        """Retrieve volume data from a mesh

        Returns:
            array: Volumetric data
        """
        # tstart = time()
        obj = bpy.data.objects.get('%s_volume' % self.label)
        if obj is None:
            return None
        n = len(obj.data.vertices)
        volume = np.empty(n*3, dtype=np.float64)
        obj.data.vertices.foreach_get('co', volume)
        volume = volume.reshape(-1, 1)
        shape = self.volumeShape
        npoint = np.product(shape)
        volume = volume[:npoint]
        volume = volume.reshape(shape)
        # print('Read volume: {0:1.2f}'.format(time() - tstart))
        return volume

    @volume.setter
    def volume(self, volume):
        self.build_volume(volume)

    @property
    def volumeShape(self):
        if "%s_volume" % self.label not in bpy.data.objects:
            return 0
        return bpy.data.objects["%s_volume" % self.label].batoms.bvolume.shape

    @volumeShape.setter
    def volumeShape(self, volumeShape):
        bpy.data.objects["%s_volume" %
                         self.label].batoms.bvolume.shape = volumeShape

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
        scale = self.attributes['scale']
        return scale

    def set_scale(self, scale):
        """
        """
        if isinstance(scale, (int, float)):
            scale = np.ones(len(self))*scale
        # for species
        elif isinstance(scale, dict):
            species = self.attributes['species']
            scale0 = self.scale
            for key, value in scale.items():
                scale0[np.where(species == key)] = value
            scale = scale0
        self.set_attributes({'scale': scale})

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
    def model_style(self):
        return self.get_model_style()

    @model_style.setter
    def model_style(self, model_style):
        self.set_model_style(model_style)

    def get_model_style(self):
        return self.attributes['model_style']

    def set_model_style(self, model_style):
        model_style = {'model_style': np.ones(
            len(self), dtype=int)*int(model_style)}
        self.set_attributes(model_style)
        self.draw(model_style['model_style'], draw_isosurface=False)
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
        self.coll.batoms.show_unit_cell = show_unit_cell
        self.cell.draw_cell()

    @property
    def radius(self):
        return self.get_radius()

    def get_radius(self):
        radius = {}
        instancers = self.species.instancers
        for sp in self.species:
            radius[sp.name] = instancers[sp.name].batoms.batom.radius
        return radius

    @property
    def radius_style(self):
        radius_style = {}
        for name, sp in self.species.items():
            radius_style[name] = sp.radius_style
        return radius_style

    @radius_style.setter
    def radius_style(self, radius_style):
        for name, sp in self.species.items():
            sp.radius_style = radius_style

    @property
    def color_style(self):
        color_style = {}
        for name, sp in self.species.items():
            color_style[name] = sp.color_style
        return color_style

    @color_style.setter
    def color_style(self, color_style):
        for name, sp in self.species.items():
            sp.color_style = color_style

    @property
    def radii_vdw(self):
        from ase.data import vdw_radii, chemical_symbols
        object_mode()
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
        object_mode()
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
        # main elements
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
        return list(self.cell.obj.batoms.bcell.pbc)

    def set_pbc(self, pbc):
        if isinstance(pbc, bool):
            pbc = [pbc]*3
        self.cell.obj.batoms.bcell.pbc = pbc

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

    def delete_verts(self, index=[]):
        """
        delete verts
        """
        import bmesh
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
        self.delete_verts(index)

    def __delitem__(self, index):
        """
        """
        self.delete(index)

    def draw_constraints(self):
        """Draw label for constraints
        # TODO
        """
        pass

    def get_frames(self):
        """
        """
        frames = self.get_obj_frames(self.obj)
        return frames

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0:
            return
        name = self.label
        obj = self.obj
        self.set_obj_frames(name, obj, frames)

    def __getitem__(self, index):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """

        if isinstance(index, int):
            batom = Batom(self.label, index, batoms=self)
            # bpy.ops.object.mode_set(mode=mode)
            return batom
        if isinstance(index, str):
            bspecies = self.species[index]
            return bspecies
        else:
            return self.positions[index]

    def __setitem__(self, index, value):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[index] = value
        self.set_positions(positions)

    def __imul__(self, m):
        """
        In-place repeat of atoms.

        >>> from batoms.batom import Batom
        >>> c = Batom('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
        >>> c.repeat([3, 3, 3], np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]))
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
        positions = self.positions
        positions = np.tile(positions, (M,) + (1,) *
                            (len(positions.shape) - 1))
        attributes = self.attributes
        for key, data in attributes.items():
            attributes[key] = np.tile(
                data, (M,) + (1,) * (len(data.shape) - 1))
        i0 = 0
        n1 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
                    n1 += 1
        self.add_vertices(positions[n:])
        self.set_attributes(attributes)
        self.cell.repeat(m)
        self.update_gn_cell()
        # if self.volume is not None:
        # self.volume = np.tile(self.volume, m)
        # repeat frames
        frames_new = []
        if self.nframe > 1:
            for i in range(0, self.nframe):
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
        self.set_frames(frames_new)
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

    def copy(self, label):
        """
        Return a copy.

        label: str
            The label of the copy.

        For example, copy H species:

        >>> h_new = h.copy(label = 'h_new', species = 'H')
        # TODO Support copy of other prroperties: materials, info ...
        """
        object_mode()
        # copy object first
        arrays = self.arrays
        batoms = self.__class__(label=label,
                                species=arrays['species'],
                                positions=arrays['positions'],
                                pbc=self.pbc,
                                cell=self.cell.array,
                                )
        batoms.translate([2, 2, 2])
        return batoms

    def extend(self, other):
        """
        Extend batom object by appending batoms from *other*.

        >>> slab = au111 + co
        todo: merge bonds setting
        """
        # could also use self.add_vertices(other.positions)
        object_mode()
        # merge bondsetting
        self.bonds.setting.extend(other.bonds.setting)
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
        # update species and species_indextsa
        self._species.extend(other._species)
        self.selects.add(self.label[0:min(4, len(self.label))], indices1)
        self.selects.add(other.label[0:min(4, len(other.label))], indices2)
        # remove shape key from mol
        sp = self.obj.data.shape_keys.key_blocks.get('Basis_%s' % other.label)
        self.obj.shape_key_remove(sp)

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
        # if kind exists, merger, otherwise build a new kind and add.
        object_mode()
        if isinstance(species, str):
            species = [species, {'elements': {species.split('_')[0]: 1.0}}]
        if species[0] not in self.species:
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
        # print(self.species)
        # for sp in self.species:
        # self.bonds.setting.add([(species[0], sp.name)])
        # self.polyhedrasetting.add([species[0]])

    def add_vertices(self, positions):
        """
        Todo: find a fast way.
        """
        import bmesh
        object_mode()
        positions = positions - self.obj.location
        bm = bmesh.new()
        bm.from_mesh(self.obj.data)
        bm.verts.ensure_lookup_table()
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.obj.data)
        bm.clear()

    def get_cell(self):
        if self.label not in bpy.data.collections:
            return None
        bcell = bpy.data.collections['%s_cell' % self.label]
        cell = np.array(
            [bcell.matrix_world @ bcell.data.vertices[i].co for i in range(3)])
        return cell

    def make_real(self):
        """ Make instancing object real
        # TODO: add make real to geometry node
        """
        self.select = True
        bpy.ops.object.duplicates_make_real()

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

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned.
        """
        return self.as_ase().get_center_of_mass(scaled=scaled)

    def get_center_of_geometry(self, colls=None):
        """
        """
        vertices = self.get_all_vertices(colls=colls, cell=self.show_unit_cell)
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

    def set_show(self, show, only_atoms=True):
        #
        if not only_atoms:
            names = self.coll.all_objects.keys()
            for name in names:
                obj = bpy.data.objects.get(name)
                obj.hide_render = not show
                obj.hide_set(not show)
        #
        if isinstance(show, (int, float)):
            show = np.ones(len(self))*show
        self.set_attributes({'show': show})

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

    def get_all_vertices(self, colls=None, cell=True):
        """
        Get position of all vertices from all mesh in batoms.
        Used for plane boundary and calc_camera_data
        """
        positions = self.positions
        # isosurface, plane
        if colls is None:
            colls = subcollections
        for coll in colls:
            if not cell and coll == 'cell':
                continue
            if 'atom' == coll:
                continue
            if 'instancer' in coll:
                continue
            if 'render' in coll:
                continue
            for obj in self.coll.children['%s_%s' % (self.label,
                                                     coll)].all_objects:
                if obj.type != 'MESH':
                    continue
                if 'volume' in obj.name:
                    continue
                if 'instancer' in obj.name:
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

    def get_canvas_box(self, direction=[0, 0, 1], padding=None, colls=None):
        """
        Calculate the canvas box from [0, 0, 1] and other direction.

        """
        from batoms.utils import get_canvas
        vertices = self.get_all_vertices(colls=colls, cell=self.show_unit_cell)
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
        if self._bonds is not None:
            return self._bonds
        bonds = Bonds(self.label, bond_datas=default_bond_datas,
                      batoms=self)
        self.bonds = bonds
        return bonds

    @bonds.setter
    def bonds(self, bonds):
        self._bonds = bonds

    @property
    def polyhedras(self):
        """polyhedras object."""
        if self._polyhedras is not None:
            return self._polyhedras
        polyhedras = Polyhedras(self.label,
                                polyhedra_datas=default_polyhedra_datas,
                                batoms=self)
        self.polyhedras = polyhedras
        return polyhedras

    @polyhedras.setter
    def polyhedras(self, polyhedras):
        self._polyhedras = polyhedras

    @property
    def boundary(self):
        """boundary object."""
        if self._boundary is not None:
            return self._boundary
        boundary = Boundary(self.label, boundary_datas=default_boundary_datas,
                            batoms=self,
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
        self.boundary[:] = boundary

    def get_arrays_with_boundary(self):
        """
        get arrays with boundary atoms
        """
        # arrays = self.arrays
        arrays_b = None
        if self._boundary is not None:
            # arrays_b = {
            #         'indices': np.arange(natom),
            #         'species': arrays['species'],
            #         'positions': arrays['positions'],
            #         'offsets': np.zeros((natom, 3)),
            #         }
            arrays_b = self.boundary.boundary_data
            # arrays_b['positions'] = np.append(arrays_b['positions'],
            #             boundary_data['positions'], axis = 0)
            # arrays_b['indices'] = np.append(arrays_b['indices'],
            #             boundary_data['indices'])
            # arrays_b['species'] = np.append(arrays_b['species'],
            #             boundary_data['species'])
            # arrays_b['offsets'] = np.append(arrays_b['offsets'],
            #             boundary_data['offsets'], axis = 0)
        else:
            arrays_b = None
        return arrays_b

    @property
    def render(self):
        """Render object."""
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

    def draw(self, model_style=None, draw_isosurface=True):
        """
        Draw atoms, bonds, polyhedra, .

        Parameters:

        model_style: str
        draw_isosurface: bool
        """
        # self.draw_cell()
        self.draw_space_filling()
        self.draw_ball_and_stick()
        self.draw_polyhedra()
        self.draw_wireframe()

    def draw_space_filling(self, scale=1.0):
        mask = np.where(self.model_style == 0, True, False)
        self.set_attribute_with_indices('scale', mask, scale)

    def draw_ball_and_stick(self, scale=0.4):
        mask = np.where(self.model_style == 1, True, False)
        self.set_attribute_with_indices('scale', mask, scale)
        self.bonds.update()

    def draw_polyhedra(self, scale=0.4):
        mask = np.where(self.model_style == 2, True, False)
        self.polyhedras.update()
        self.set_attribute_with_indices('show', mask, True)
        if self.polyhedra_style == 0:
            self.set_attribute_with_indices('scale', mask, scale)
            # self.bonds.update()
        if self.polyhedra_style == 1:
            self.set_attribute_with_indices('scale', mask, scale)
        elif self.polyhedra_style == 2:
            for b in self.bonds.setting:
                if b.polyhedra:
                    mask1 = np.where(
                        self.attributes['species'] == b.species1, True, False)
                    self.set_attribute_with_indices('scale', mask1, scale)
                    mask[mask1] = False
            scale = 0
            self.set_attribute_with_indices('scale', mask, scale)
            # self.set_attribute_with_indices('show', mask, False)
        elif self.polyhedra_style == 3:
            scale = 0
            self.set_attribute_with_indices('scale', mask, scale)
            # self.set_attribute_with_indices('show', mask, False)

    def draw_wireframe(self):
        mask = np.where(self.model_style == 3, True, False)
        self.set_attribute_with_indices('show', mask, 0)
        self.set_attribute_with_indices('scale', mask, 0.0001)
        # self.update(mask)

    def as_ase(self, local=True):
        """
        local: bool
            if True, use the origin of uint cell as the origin
        """
        from ase import Atoms
        arrays = self.arrays
        positions = arrays['positions']
        if local:
            positions -= self.cell.origin
        atoms = Atoms(symbols=arrays['elements'],
                      positions=positions,
                      cell=self.cell, pbc=self.pbc)
        for name, array in arrays.items():
            if name in ['elements', 'positions']:
                continue
            atoms.set_array(name, np.array(array))

        return atoms

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
