"""Definition of the Batoms class.

This module defines the Batoms object in the batoms package.

"""
import bpy
import bmesh
from batoms.batom import Batom
from batoms.bonds import Bonds, default_bond_datas
from batoms.polyhedras import Polyhedras, default_polyhedra_datas
from batoms.boundary import Boundary, default_boundary_datas
from batoms.bspecies import Bspecies
from batoms.cell import Bcell
from batoms.render.render import Render
from batoms.bselect import Selects
from batoms.tools import get_default_species_data, number2String, string2Number, read_from_ase, read_from_pymatgen
from batoms.base import BaseCollection
from batoms.isosurfacesetting import IsosurfaceSetting
from batoms.planesetting import PlaneSetting
from batoms.mssetting import MSsetting
from batoms.ribbon import Ribbon
from batoms.bdraw import draw_cylinder
from batoms.butils import object_mode, show_index, get_nodes_by_name
import numpy as np
from time import time

shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]

default_attributes = [
        ['select', 'INT'],
        ['species_index', 'INT'], 
        ['species', 'STRING'], 
        ['show', 'BOOLEAN'], 
        ['scale', 'FLOAT'], 
        ['model_style', 'INT'],
        ['radius_style', 'INT'],
        ]
    
default_GroupInput = [
        ['select', 'INT'],
        ['species_index', 'INT'], 
        ['show', 'BOOLEAN'], 
        ['scale', 'FLOAT'], 
        ]

subcollections = ['instancer', 'bond', 'polyhedra', 'surface', 'ribbon', 'plane']

class Batoms(BaseCollection):
    """Batom Class
    
    The Batoms object is a interface to a batoms object in Blender.

    Parameters:

    label: str
        Name for the object in Blender.
    species: list of str
        ['O', 'H', 'H']
        ['Fe_u', 'Fe_d', 'O']
    positions: array
        positions
    attributes: array
        eg. index, chianid
    locations: array
        The object's origin location in global coordinates.
    elements: dict
        elements for each species, includes fractional Occupancy
        str: 'O'
        dict:{'Al': {'Al': 0.887, 'Si': 0.333}
              'O' : {'O': 1.0}}
    pbc: Bool or three Bool
        Periodic boundary conditions. Examples: True,
        False, (True, False, False).  Default value: False.
    cell: 3x3 matrix or length 3 or 6 vector
        Unit cell.
    segments: list of 2 Int
        Value should be int, and in [3, 100000]
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
    model_style: int
        enum in [0, 1, 2, 3], Default value: 0
    atoms: ase.atoms.Atoms object or a list of ase.atoms.Atoms object
        or pymatgen structure object
    boundary:  list 
        search atoms at the boundary
    
    Examples:

    >>> from batoms import Batoms
    >>> h2o = Batom('h2o', species = ['O', 'H', 'H'], 
                elements={'O':{'O':0.8, 'N': 0.2}, 'H':{'H': 0.8}}, 
                positions= [[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]])

    """
    def __init__(self, 
                label = 'batoms',
                species = [],
                positions = [],
                attributes = {},
                species_props = {},
                info = {},
                pbc = False, cell = None,
                location = np.array([0, 0, 0]),
                boundary = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                show_unit_cell = True,
                volume = None,
                scale = 1.0, 
                model_style = 0, 
                polyhedra_style = 0, 
                from_ase = None, 
                from_pymatgen = None, 
                metaball = False,
                movie = True,
                segments = None,
                 ):
        #
        BaseCollection.__init__(self, coll_name = label)
        if from_ase is not None:
            species, positions, attributes, cell, pbc, info = read_from_ase(from_ase)
        if from_pymatgen is not None:
            species, positions, attributes, cell, pbc, info = read_from_pymatgen(from_pymatgen)            
        if len(species) == 0 and self.check_batoms(label):
            self.from_batoms(label)
        else:
            self.set_collection(label, boundary)
            self._cell = Bcell(label, cell, batoms = self)
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
                species_props = {sp: {'elements':{sp.split('_')[0]:1.0}} for sp in species}
            self._species = Bspecies(label, label, species_props, self, segments = segments)
            self.build_geometry_node()
            self.selects.add('sel0', np.arange(len(self)))
            if isinstance(scale, (int, float)):
                scale = np.ones(natom)*scale
            show = np.ones(natom, dtype = int)
            # elements = self.check_elements(elements)
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
        self.isosurfacesetting = IsosurfaceSetting(self.label, batoms = self)
        self.planesetting = PlaneSetting(self.label, batoms = self)
        self.mssetting = MSsetting(self.label, probe = 1.4, batoms = self)
        self.ribbon = Ribbon(self.label, batoms = self, datas = info, update = True)
        self._render = None
        self._bonds = None
        self._polyhedras = None
        self._boundary = None
        show_index()
    
    def set_collection(self, label, boundary = [0, 0, 0]):
        """
        build main collection and its child collections.
        """
        if bpy.data.collections.get(label):
                raise Exception("Failed, the name %s already in use!"%label)
        coll = bpy.data.collections.new(label)
        bpy.data.scenes['Scene'].collection.children.link(coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s'%(label, sub_name))
            coll.children.link(subcoll)
        coll.batoms.flag = True
        coll.batoms.label = label
        coll.batoms.boundary = boundary
        
    def build_object(self, label, positions, location = [0, 0, 0]):
        """
        build child object and add it to main objects.
        """
        self.obj_name = label
        if label not in bpy.data.objects:
            mesh = bpy.data.meshes.new(label)
            # Add attributes
            for attribute in default_attributes:
                mesh.attributes.new(name = attribute[0], type = attribute[1], domain = 'POINT')
            obj = bpy.data.objects.new(label, mesh)
            obj.data.from_pydata(positions, [], [])
            obj.location = location
            obj.batoms.batom.flag = True
            obj.batoms.batom.label = label
            self.coll.objects.link(obj)
        elif bpy.data.objects[label].batoms.batom.flag:
            obj = bpy.data.objects[label]
            self.coll.objects.link(obj)
            print("Object %s already exists, Load it!"%label)
        else:
            raise Exception("Failed, the name %s already in use and is not Batom object!"%label)

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        name = 'GeometryNodes_%s'%self.label
        modifier = self.obj.modifiers.new(name = name, type = 'NODES')
        modifier.node_group.name = name
        GroupInput = modifier.node_group.nodes.get('Group Input')
        # add new output sockets
        for att in default_GroupInput:
            GroupInput.outputs.new(type = att[1], name = att[0])
        # the above codes not works. maybe bug in blender, 
        # we add this, maybe deleted in the future
        for i in range(1, 5):
            test = get_nodes_by_name(modifier.node_group.nodes, 
                            'BooleanMath_%s'%i,
                            'FunctionNodeCompareFloats')
            modifier.node_group.links.new(GroupInput.outputs[i], test.inputs[0])
        #
        i = 2
        for att in default_GroupInput:
            modifier['Input_%s_use_attribute'%i] = 1
            modifier['Input_%s_attribute_name'%i] = att[0]
            i += 1
        gn = modifier
        # print(gn.name)
        # print(GroupInput.outputs[:])
        GroupOutput = gn.node_group.nodes.get('Group Output')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])
        # set positions
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_PositionBatoms'%(self.label),
                        'GeometryNodeInputPosition')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label, 
                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        # gn.node_group.links.new(VectorAdd.outputs[0], SetPosition.inputs['Position'])
        # todo: use cell object directly.
        cell = self.cell.array
        if np.isclose(np.linalg.det(self.cell.array), 0):
            cell = np.eye(3)
        icell = np.linalg.inv(cell)
        scaledPositionsNode = self.vectorDotMatrix(gn, PositionBatoms, icell, 'scaled')
        VectorWrap = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_VectorWrap'%(self.label),
                        'ShaderNodeVectorMath')
        VectorWrap.operation = 'WRAP'
        VectorWrap.inputs[1].default_value = [1, 1, 1]
        VectorWrap.inputs[2].default_value = [0, 0, 0]
        gn.node_group.links.new(scaledPositionsNode.outputs[0], VectorWrap.inputs['Vector'])
        PositionsNode = self.vectorDotMatrix(gn, VectorWrap, cell, '')
        gn.node_group.links.new(PositionsNode.outputs[0], SetPosition.inputs['Position'])
        self.wrap = self.pbc
    def vectorDotMatrix(self, gn, vectorNode, matrix, name):
        """
        """
        CombineXYZ = get_nodes_by_name(gn.node_group.nodes,
                        '%s_CombineXYZ_%s'%(self.label, name),
                        'ShaderNodeCombineXYZ')
        #
        VectorDot = []
        for i in range(3):
            tmp = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_VectorDot%s_%s'%(self.label, i, name),
                        'ShaderNodeVectorMath')
            tmp.operation = 'DOT_PRODUCT'
            VectorDot.append(tmp)
            tmp.inputs[1].default_value = matrix[i]
            gn.node_group.links.new(vectorNode.outputs[0], tmp.inputs[0])
            gn.node_group.links.new(tmp.outputs['Value'], CombineXYZ.inputs[i])
        return CombineXYZ
        
    def add_geometry_node(self, spname, selname):
        """
        """
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label)
        CompareSelect = get_nodes_by_name(gn.node_group.nodes, 
                    'select_%s_%s'%(self.label, selname),
                    'FunctionNodeCompareFloats')
        CompareSelect.operation = 'EQUAL'
        # CompareSelect.data_type = 'INT'
        CompareSelect.inputs[1].default_value = string2Number(selname)
        gn.node_group.links.new(GroupInput.outputs[1], CompareSelect.inputs[0])
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
        ObjectInfo.inputs['Object'].default_value = self.species.instancers[selname][spname]
        #
        BoolSelectSpecies = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_%s_%s_0'%(self.label, selname, spname),
                        'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes, 
                    'BooleanMath_%s_%s_%s_1'%(self.label, selname, spname),
                    'FunctionNodeBooleanMath')
        #
        # gn.node_group.links.new(GroupInput.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(SetPosition.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[1], CompareSelect.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[2], CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[4], InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSelect.outputs[0], BoolSelectSpecies.inputs[0])
        gn.node_group.links.new(CompareSpecies.outputs[0], BoolSelectSpecies.inputs[1])
        gn.node_group.links.new(BoolSelectSpecies.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])        

    def build_volume(self, volume):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        # remove old volume point
        tstart = time()
        if volume is None: return
        name = "%s_volume"%self.label
        if name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[name], do_unlink = True)
        shape = volume.shape
        volume = volume.reshape(-1, 1)
        npoint = len(volume)
        dn = 3 - npoint % 3
        verts = np.append(volume, np.zeros((dn, 1)), axis = 0)
        verts = verts.reshape(-1, 3)
        mesh = bpy.data.meshes.new("%s_volume"%self.label)
        mesh.from_pydata(verts, [], [])  
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        obj.data = mesh
        obj.batoms.bvolume.is_bvolume = True
        obj.batoms.bvolume.shape = shape
        self.coll.objects.link(obj)
        obj.hide_set(True)
        obj.hide_render = True
        print('Draw volume: {0:1.2f}'.format(time() - tstart))

    def check_batoms(self, label):
        flag = True
        if label not in bpy.data.collections:
            flag = False
        elif not bpy.data.collections[label].batoms.flag:
            flag = False
        if label not in bpy.data.objects:
            flag = False
        elif not bpy.data.objects[label].batoms.batom.flag:
            flag = False
        return flag

    def from_batoms(self, label):
        
        self.coll_name = label
        self.obj_name = label
        self._cell = Bcell(label = label)
        self._species = Bspecies(label, label, {}, self)
        self.selects = Selects(label, self)
    
    @property
    def obj(self):
        return self.get_obj()
    
    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.'%self.obj_name)
        return obj

    @property    
    def volume(self):
        return self.get_volume()
    
    def get_volume(self):
        tstart = time()
        obj = bpy.data.objects.get('%s_volume'%self.label)
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
        if "%s_volume"%self.label not in bpy.data.objects:
            return 0
        return bpy.data.objects["%s_volume"%self.label].batoms.bvolume.shape
    
    @volumeShape.setter    
    def volumeShape(self, volumeShape):
        bpy.data.objects["%s_volume"%self.label].batoms.bvolume.shape = volumeShape

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
        scale = np.ones(len(self))*scale
        # for species
        if isinstance(scale, dict):
            scale_dict = {}
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
        model_style = {'model_style': np.ones(len(self))*int(model_style)}
        self.set_attributes(model_style)
        self.draw(model_style, draw_isosurface = False)
    
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
        for sel in self.selects:
            radius[sel.name] = {}
            for sp in self.species:
                radius[sel.name][sp.name] = instancers[sel.name][sp.name].batoms.batom.radius
        return radius
    
    @property
    def radius_style(self):
        return self.get_radius_style()
    
    @radius_style.setter
    def radius_style(self, radius_style):
        self.set_radius_style(radius_style)
    
    def get_radius_style(self):
        return self.attributes['radius_style']
    
    def set_radius_style(self, radius_style):
        self.coll.batoms.radius_style = str(radius_style)
        scale = self.scale
        species_props = get_default_species_data(self.species,
                                radius_style = radius_style)
        # print(species_props)
        self.species

    @property
    def radii_vdw(self):
        from ase.data import vdw_radii, chemical_symbols
        object_mode()
        radii = []
        elements = self.arrays['elements']
        for element in elements:
            if element == 'X': continue
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
    
    def get_scaled_positions(self, cell = None):
        """
        Get array of scaled_positions.
        """
        from ase.cell import Cell
        if not cell:
            cell = self.cell
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
        arrays.update({'positions': self.positions})
        # radius
        radius = self.radius
        arrays.update({'radius': np.zeros(len(self))})
        for sel, data in radius.items():
            for sp, value in data.items():
                mask = np.where((arrays['species'] == sp) & (arrays['select'] == string2Number(sel)))
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

        Parameters:

        cell: 

        Examples:

        """
        from ase.cell import Cell
        from ase.geometry.cell import complete_cell

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
    def gnodes(self):
        return self.get_gnodes()
    
    @gnodes.setter
    def gnodes(self, gnodes):
        self.set_gnodes(gnodes)
    
    def get_gnodes(self):
        name = 'GeometryNodes_%s'%self.label
        gnodes = self.obj.modifiers.get(name)
        return gnodes
    
    def set_gnodes(self, gnodes):
        pass
    
    @property
    def mnode(self):
        return self.get_mnode()
    
    @mnode.setter
    def mnode(self, mnode):
        self.set_mnode(mnode)
    
    def get_mnode(self):
        return self.materials[self.main_elements].node_tree.nodes
    
    def set_mnode(self, mnode):
        for key, value in mnode.items():
            self.materials[self.main_elements].node_tree.nodes['Principled BSDF'].inputs[key].default_value = value
    
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
        bpy.context.view_layer.update()
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0 : return
        obj = self.obj
        base_name = 'Basis_%s'%self.label
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
        
        if isinstance(index, int):
            batom = Batom(self.label, index, batoms=self)
            # bpy.ops.object.mode_set(mode=mode)
            return batom
        if isinstance(index, str):
            bspecies = Bspecies(self.label, index, batoms=self)
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
        positions = np.tile(positions, (M,) + (1,) * (len(positions.shape) - 1))
        attributes = self.attributes
        for key, data in attributes.items():
            attributes[key] = np.tile(data, (M,) + (1,) * (len(data.shape) - 1))
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
        # if self.volume is not None:
            # self.volume = np.tile(self.volume, m)
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

    def repeat(self, m):
        """
        """
        self *= m
        return self
    
    def __mul__(self, m):
        self.repeat(m)
        return self

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
        # geometry nodes
        obj.modifiers[0].node_group = self.obj.modifiers[0].node_group.copy()
        bpy.data.collections['Collection'].objects.link(obj)
        # species instancer
        self.species.copy(name)
        # cell
        batoms = self.__class__(name)
        batoms._cell = self.cell.copy(name)
        batoms.translate([2, 2, 2])
        return batoms
    
    def extend(self, other):
        """
        Extend batom object by appending batoms from *other*.
        
        >>> from batoms.batoms import Batom
        >>> h2o = Batoms('h2o', species = ['O', 'H', 'H'], positions= [[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]])
        >>> co = Batoms('co', species = ['C', 'O'], positions= [[1, 0, 0], [2.3, 0, 0]])
        >>> h2o = h2o + co
        """
        # could also use self.add_vertices(other.positions)
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        other.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.join()
        # update species and species_indextsa
        self._species.extend(other._species)
        self.build_geometry_node()
    
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
        text.append('species=%s, '%(list(self.species.species)))
        text.append('cell={0}, '.format(self.cell))
        text.append('pbc={0}'.format(self.pbc))
        # text.append('positions={0}'.format(self.positions))
        text = "".join(text)
        text = "Batoms(%s)"%text
        return text
    
    def replace(self, indices, species):
        """Replace species.
        Parameters:

        indices: list
            indices of atoms will be replaced.
        species: str
            atoms will be changed to this species.
        
        >>> from ase.build import molecule, fcc111
        >>> from batoms.batoms import Batoms
        >>> pt111 = fcc111('Pt', (5, 5, 4), vacuum = 5.0)
        >>> pt111 = Batoms('pt111', from_ase = pt111)
        >>> pt111.replace([0], 'Au')
        >>> pt111.replace(range(20), 'Au')
        """
        # if kind exists, merger, otherwise build a new kind and add.
        object_mode()
        if isinstance(species, str):
            species = [species, {'elements': {species.split('_')[0]: 1.0}}]
        if species[0] not in self.species:
            self.species[species[0]] = species[1]
            # add geometry node, select
            select_array = self.attributes['select'][indices]
            selnames = np.unique(select_array)
            for selname in selnames:
                self.add_geometry_node(species[0], number2String(selname))
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
    
    def get_distances(self, i, indices, mic=False):
        """
        Return distances of atom No.i with a list of atoms.

        Use mic=True to use the Minimum Image Convention.

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
        return self.atoms.get_center_of_mass(scaled = scaled)

    def get_center_of_geometry(self, colls = None):
        """
        """
        vertices = self.get_all_vertices(colls = colls, cell = self.show_unit_cell)
        canvas = np.zeros([2, 3])
        canvas[0] = np.min(vertices, axis = 0)
        canvas[1] = np.max(vertices, axis = 0)
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
    
    def set_show(self, show, only_atoms = True):
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
        return self.attributes['wrap']
    
    def set_wrap(self, wrap):
        #
        if isinstance(wrap, bool):
            wrap = [wrap]*3
        wrap = np.array(wrap)
        if not wrap.any():
            # switch off
            n = len(self.gnodes.node_group.nodes['%s_CombineXYZ_'%self.label].outputs['Vector'].links)
            if n > 0:
                l = self.gnodes.node_group.nodes['%s_CombineXYZ_'%self.label].outputs['Vector'].links[0]
                self.gnodes.node_group.links.remove(l)
        else:
            self.gnodes.node_group.links.new(self.gnodes.node_group.nodes['%s_CombineXYZ_'%self.label].outputs[0], 
                    self.gnodes.node_group.nodes['%s_SetPosition'%self.label].inputs['Position'])
        self.gnodes.node_group.update_tag()
        #todo: support selection
        
    def get_spacegroup_number(self, symprec = 1e-5):
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
        import spglib
        atoms = self.atoms
        lattice = atoms.cell
        points = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, points, numbers)
        lattice, points,  numbers= spglib.find_primitive(cell)
        atoms = Atoms(numbers=numbers, scaled_positions=points, cell = lattice)
        return atoms
    
    def get_all_vertices(self, colls = None, cell = True):
        """
        Get position of all vertices from all mesh in batoms.
        Used for plane boundary and calc_camera_data
        """
        positions = self.positions
        # isosurface, plane
        if colls is None:
            colls = subcollections
        for coll in colls:
            if not cell and coll == 'cell': continue
            if 'atom' == coll: continue
            if 'instancer' in coll: continue
            if 'render' in coll: continue
            for obj in self.coll.children['%s_%s'%(self.label, coll)].all_objects:
                if obj.type != 'MESH': continue
                if 'volume' in obj.name: continue
                if 'instancer' in obj.name: continue
                n = len(obj.data.vertices)
                vertices = np.empty(n*3, dtype=np.float64)
                obj.data.vertices.foreach_get('co', vertices)  
                vertices = vertices.reshape((n, 3))
                vertices = np.append(vertices, np.ones((n, 1)), axis = 1)
                mat= np.array(obj.matrix_world)
                vertices = mat.dot(vertices.T).T
                # (natom, 4) back to (natom, 3)
                vertices = vertices[:, :3]
                positions = np.concatenate((positions, vertices), axis = 0)
        return positions
    
    def get_canvas_box(self, direction = [0, 0, 1], padding = None, colls = None):
        """
        Calculate the canvas box from [0, 0, 1] and other direction.

        """
        from batoms.tools import get_canvas
        vertices = self.get_all_vertices(colls = colls, cell = self.show_unit_cell)
        canvas = get_canvas(vertices, direction = direction, padding = padding)
        width = canvas[1][0] - canvas[0][0]
        height = canvas[1][1] - canvas[0][1]
        depth = canvas[1][2] - canvas[0][2]
        return width, height, depth
    
    def lock_to_camera(self, obj):
        from batoms.butils import lock_to
        for sel, instancers in self.species.instancers.items():
            for sp, instancer in instancers.items():
                lock_to(instancer, obj, location = False, rotation = True)
    
    @property
    def bonds(self):
        """bonds object."""
        if self._bonds is not None:
            return self._bonds
        bonds = Bonds(self.label, bond_datas = default_bond_datas, 
                    batoms = self)
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
        polyhedras = Polyhedras(self.label, polyhedra_datas = default_polyhedra_datas, 
                    batoms = self)
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
        boundary = Boundary(self.label, boundary_datas = default_boundary_datas, 
                    batoms = self)
        self.boundary = boundary
        return boundary

    @boundary.setter
    def boundary(self, boundary):
        self._boundary = boundary

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

    def get_image(self, viewport = None, engine = None, 
                    frame = 1, 
                    animation = False, 
                    output = None,
                    center = None,
                    padding = None, 
                    canvas = None, 
                    gpu = False):
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
            output = '%s.png'%self.label
        if self.render is None:
            raise RuntimeError('Batoms object has no render.')
        self.render.run_render = True
        if viewport is not None:
            self.render.viewport = viewport
        if engine is not None:
            self.render.engine = engine
        self.render.gpu = gpu
        image = self.render.run(self, frame = frame, 
                                    animation = animation, 
                                    output = output,
                                    center = center,
                                    padding = padding,
                                    canvas = canvas, 
                                    )
        return image
    
    def make_real(self):
        """
        """
        self.select = True
        bpy.ops.object.duplicates_make_real()

    def draw(self, model_style = None, draw_isosurface = True):
        """
        Draw atoms, bonds, polyhedra, .

        Parameters:

        model_style: str
        draw_isosurface: bool
        """
        from batoms.butils import clean_coll_objects
        # clean_coll_objects(self.coll.children['%s_bond'%self.label], 'bond')
        clean_coll_objects(self.coll.children['%s_polyhedra'%self.label], 'polyhedra')
        # self.draw_cell()
        self.draw_space_filling()
        self.draw_ball_and_stick()
        if 2 in model_style:
            self.draw_polyhedra()
        self.draw_wireframe()
        """
        if model_style == 2:

        elif model_style == 3:
            self.scale = 0.01
            self.update_bonds()
        """
    
    def draw_space_filling(self, scale = 1.0):
        mask = np.where(self.model_style == 0, True, False)
        self.set_attribute_with_indices('scale', mask, scale)

    def draw_ball_and_stick(self, scale = 0.4):
        mask = np.where(self.model_style == 1, True, False)
        self.set_attribute_with_indices('scale', mask, scale)
        self.bonds.update_bonds()
    
    def draw_polyhedra(self, scale = 0.4):
        mask = np.where(self.model_style == 2, True, False)
        self.polyhedras.update_polyhedras()
        self.set_attribute_with_indices('show', mask, True)
        if self.polyhedra_style == 0:
            self.set_attribute_with_indices('scale', mask, scale)
            # self.bonds.update_bonds()
        if self.polyhedra_style == 1:
            self.set_attribute_with_indices('scale', mask, scale)
        elif self.polyhedra_style == 2:
            for b in self.bondsetting:
                if b.polyhedra:
                    mask1 = np.where(self.attributes['species'] == b.species1, True, False)
                    self.set_attribute_with_indices('scale', mask1, scale)
                    mask[mask1] = False
            self.set_attribute_with_indices('show', mask, False)
        elif self.polyhedra_style == 3:
            self.set_attribute_with_indices('show', mask, False)

    def draw_wireframe(self):
        mask = np.where(self.model_style == 3, True, False)
        self.set_attribute_with_indices('show', mask, 0)
        self.set_attribute_with_indices('scale', mask, 0.0001)
        # self.update_bonds(mask)


