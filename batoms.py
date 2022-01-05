"""Definition of the Batom class.

This module defines the Batom object in the batoms package.

"""

from numpy.core.records import array
import bpy
import bmesh
from batoms.cell import Bcell
from batoms.render.render import Render
from batoms.bselect import Selects
from batoms.tools import get_default_species_data
from batoms.base import BaseObject
from batoms.bondsetting import BondSetting, build_bondlists, calc_bond_data
from batoms.polyhedrasetting import PolyhedraSetting, build_polyhedralists
from batoms.isosurfacesetting import IsosurfaceSetting
from batoms.planesetting import PlaneSetting
from batoms.mssetting import MSsetting
from batoms.ribbon import Ribbon
from batoms.boundary import search_boundary, search_bond
from batoms.bdraw import draw_cylinder, draw_surface_from_vertices
from batoms.butils import object_mode, show_index
import numpy as np
from time import time


shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]

subcollections = ['atom', 'cell', 'bond', 'polyhedra', 'instancer', 
'instancer_atom', 'volume', 'surface', 'ghost', 'boundary', 'search', 'text', 'plane', 'ribbon']


class Batoms(BaseObject):
    """Batom Class
    
    The Batoms object is a interface to a batoms object in Blender.

    Parameters:

    label: str
        Name for the object in Blender.
    species: list of str
        ['O', 'H', 'H']
        ['Fe_up', 'Fe_down', 'O']
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
                elements = {},
                info = {},
                pbc = False, cell = None,
                location = np.array([0, 0, 0]),
                boundary = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                show_unit_cell = True,
                volume = None,
                scale = 1.0, 
                segments = [32, 16],
                shape = 0,
                subdivisions = 2,
                props = {},
                color = None,
                model_style = 0, 
                polyhedra_style = 0, 
                radius_style = 'covalent',
                color_style = 'JMOL',
                material_style = 'default',
                materials = {},
                node_inputs = None,
                instancer_data = {},
                atoms = None, 
                metaball = False,
                movie = False,
                 ):
        #
        obj_name = label
        bobj_name = 'batom'
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        if atoms is not None:
            if isinstance(atoms, list):
                frames = atoms
                atoms = frames[0]
            else:
                frames = [atoms]
            if 'ase' in str(type(atoms)):
                species, positions, attributes, cell, pbc, info = self.from_ase(atoms)
            elif 'pymatgen' in str(type(atoms)):
                species, positions, attributes, cell, pbc, info = self.from_pymatgen(atoms)
            self._frames = frames
        natom = len(positions)
        if natom > 0:
            positions = np.array(positions)
            if len(positions.shape) == 2:
                self._frames = np.array([positions])
            elif len(positions.shape) == 3:
                self._frames = positions
                positions = self._frames[0]
            else:
                raise Exception('Shape of positions is wrong!')
            if not elements:
                elements = {sp: {sp.split('_')[0]:1.0} for sp in species}
            if isinstance(scale, (int, float)):
                scale = np.ones(natom)*scale
            # elements = self.check_elements(elements)
            species_props = get_default_species_data(elements,
                                radius_style = radius_style, 
                                color_style = color_style,
                                props = props)
            self.set_collection(label)
            self.build_object(
                        positions = positions,
                        elements = elements, 
                        location = location,
                        label = label,
                        metaball = metaball)
            elements = self.elements
            species_index = [elements[sp]['index'] for sp in species]
            attributes.update({'species': species, 
                               'species_index': species_index, 
                               'scale': scale})
            self.set_attributes(attributes)
            self.build_materials(species_props, label = label, 
                            node_inputs = node_inputs, 
                            material_style = material_style, 
                            materials = materials)
            self.build_instancers(species_props, 
                            scale = scale,
                            segments = segments, 
                            subdivisions = subdivisions, shape = shapes[shape],
                            instancer_data = instancer_data)
            # self.set_frames(self._frames, only_basis = True)
            self.build_geometry_node()
            self._cell = Bcell(self.label, cell)
            # self.coll.children['%s_cell'%self.label].objects.link(self._cell.obj)
            self.set_pbc(pbc)
            self.selects = Selects(self.label, self)
        else:
            self.from_batoms(label)
        self.from_batoms(label)
        self.bondsetting = BondSetting(self.label, batoms = self)
        self.polyhedrasetting = PolyhedraSetting(self.label, batoms = self)
        self.isosurfacesetting = IsosurfaceSetting(self.label, volume = volume)
        self.planesetting = PlaneSetting(self.label, batoms = self)
        self.mssetting = MSsetting(self.label, probe = 1.4, batoms = self)
        self.ribbon = Ribbon(self.label, batoms = self, datas = info)
        show_index()
    
    def set_collection(self, label):
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
    
    def build_materials(self, species_props, label = 'batoms', node_inputs = None, 
                material_style = 'default', materials = {}):
        """
        """
        from batoms.material import create_material
        for sp, data in species_props.items():
            for ele, color in data['color'].items():
                name = '%s_%s_%s'%(self.label, sp, ele)
                if name not in bpy.data.materials:
                    create_material(name,
                                color = color,
                                node_inputs = node_inputs,
                                material_style = material_style,
                                backface_culling = True)
    
    def build_object(self, positions, elements, location, label = 'batoms', metaball = False):
        """
        build child object and add it to main objects.
        """
        if self.obj_name not in bpy.data.objects:
            mesh = bpy.data.meshes.new(self.obj_name)
            obj = bpy.data.objects.new(self.obj_name, mesh)
            obj.data.from_pydata(positions, [], [])
            obj.location = location
            obj.batoms.batom.flag = True
            obj.batoms.batom.label = label
            for sp, data in elements.items():
                spdata = obj.batoms.batom.species.add()
                spdata.name = sp
                for ele, occupancy in data.items():
                    eledata = spdata.elements.add()
                    eledata.name = ele
                    eledata.occupancy = occupancy
            self.coll.children['%s_atom'%label].objects.link(obj)
            # Assigning attributes to vertexes
        elif bpy.data.objects[self.obj_name].batoms.batom.flag:
            obj = bpy.data.objects[self.obj_name]
        else:
            raise Exception("Failed, the name %s already in use and is not Batom object!"%self.obj_name)
    
    def build_instancers(self, species_props, scale = [1, 1, 1], segments = [32, 16], subdivisions = 2, 
                        shape = 'UV_SPHERE', shade_smooth = True, instancer_data = {}):
        object_mode()
        for sp, data in species_props.items():
            name = '%s_%s_instancer'%(self.label, sp)
            radius = data['radius']
            if name in bpy.data.objects:
                obj = bpy.data.objects.get(name)
                bpy.data.objects.remove(obj, do_unlink = True)
            if sp in instancer_data:
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
                # if isinstance(scale, float):
                    # scale = [scale]*3
                # obj.scale = scale[sp]
                obj.batoms.batom.radius = radius
                obj.name = name
                obj.data.name = name
            #
            self.coll.children['%s_instancer'%self.label].objects.link(obj)
            # obj.data.materials.append(self.material)
            if shade_smooth:
                bpy.ops.object.shade_smooth()
            if shape.upper() != 'METABALL':
                obj.hide_set(True)
            # obj.parent = self.obj
            # self.obj.instance_type = 'VERTS'
            #
            self.assign_materials(sp)
            bpy.context.view_layer.update()
        return obj
    
    @staticmethod
    def from_ase(atoms):
        """
        Import structure from ASE atoms.
        """
        from batoms.tools import npbool2bool
        if 'species' not in atoms.arrays:
            atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype = 'U20'))
        if 'scale' not in atoms.arrays:
            atoms.new_array('scale', np.ones(len(atoms)))
        arrays = atoms.arrays
        positions = arrays.pop('positions')
        species = arrays.pop('species')
        info = atoms.info
        return species, positions, arrays, atoms.cell, npbool2bool(atoms.pbc), info
    
    @staticmethod
    def from_pymatgen(structure):
        """
        Import structure from Pymatgen structure.
        """
        attributes = {}
        symbols = np.array([str(site.specie.symbol) for site in structure])
        natom = len(symbols)
        attributes['species'] = np.array(symbols, dtype='U20')
        if hasattr(structure, "lattice"):
            cell = structure.lattice.matrix
            pbc = True
        else:
            cell = None
            pbc = False
        positions = [structure[i].coords for i in range(natom)]
        info = {}
        return symbols, positions, attributes, cell, pbc, info

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

    def build_geometry_node(self):
        """
        """
        from batoms.butils import get_nodes_by_name
        obj = self.obj
        name = 'GeometryNodes_%s'%self.label
        gn = obj.modifiers.get(name)
        if gn is None:
            gn = obj.modifiers.new(name = name, type = 'NODES')
        # print(gn.name)
        gn['Input_2_attribute_name'] = 'species_index'
        gn['Input_2_use_attribute'] = 1
        gn
        GroupInput = gn.node_group.nodes.get('Group Input')
        GroupInput.outputs[1].type = 'BOOLEAN'
        GroupInput.outputs[1].name = 'Selection'
        # print(GroupInput.outputs[:])
        GroupOutput = gn.node_group.nodes.get('Group Output')
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        'JoinGeometry_%s'%self.label, 
                        'GeometryNodeJoinGeometry')
        attributes = self.attributes
        for sp, data in self.elements.items():
            mask = np.where(attributes == sp)[0]
            InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                        'InstanceOnPoint_%s'%sp, 
                        'GeometryNodeInstanceOnPoints')
            ObjectInfo = get_nodes_by_name(gn.node_group.nodes, 
                        'ObjectInfo_%s_%s'%(self.label, sp),
                        'GeometryNodeObjectInfo')
            ObjectInfo.inputs['Object'].default_value = self.instancers[sp]
            #
            CompareFloats = get_nodes_by_name(gn.node_group.nodes, 
                        'CompareFloats_%s_%s'%(self.label, sp),
                        'FunctionNodeCompareFloats')
            CompareFloats.operation = 'EQUAL'
            CompareFloats.inputs[1].default_value = data['index']
            BooleanMath = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_%s'%(self.label, sp),
                        'FunctionNodeBooleanMath')
            BooleanMath.inputs[1].default_value = True
            gn.node_group.links.new(GroupInput.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
            gn.node_group.links.new(GroupInput.outputs[1], CompareFloats.inputs[0])
            gn.node_group.links.new(CompareFloats.outputs[0], BooleanMath.inputs[0])
            gn.node_group.links.new(BooleanMath.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
            gn.node_group.links.new(ObjectInfo.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
            gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])
            gn.node_group.links.new(JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])

    def assign_materials(self, species):
        # sort element by occu
        mesh = self.instancers[species].data
        mesh.materials.clear()
        sorted_ele = sorted(self.elements[species]['elements'].items(), key=lambda x: -x[1])
        materials = self.materials[species]
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

    def from_batoms(self, label):
        if label not in bpy.data.objects:
            raise Exception("%s is not a object!"%label)
        elif not bpy.data.objects[label].batoms.batom.flag:
            raise Exception("%s is not Batom object!"%label)
        obj = bpy.data.objects[label]
        self._cell = Bcell(label = label)
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        return bpy.data.collections.get(self.label)
    
    @property
    def instancers(self):
        return self.get_instancers()
    
    def get_instancers(self):
        instancers = {}
        for sp in self.elements:
            name = '%s_%s_instancer'%(self.label, sp)
            instancers[sp] = bpy.data.objects.get(name)
        return instancers
    
    @property
    def materials(self):
        return self.get_materials()
    
    def get_materials(self):
        materials = {}
        for sp, data in self.elements.items():
            materials[sp] = {}
            for ele in data['elements']:
                name = '%s_%s_%s'%(self.label, sp, ele)
                mat = bpy.data.materials.get(name)
                materials[sp][ele] = mat
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
        scale = {}
        for sp in self.elements:
            name = 'InstanceOnPoint_%s'%sp
            node = self.gnodes.get(name)
            scale[sp] = node.inputs['Scale'].default_value
        return scale
    
    def set_scale(self, scale):
        if isinstance(scale, float) or isinstance(scale, int):
            scale_dict = {}
            for sp in self.elements:
                scale_dict[sp] = [scale]*3
        for sp in self.elements:
            name = 'InstanceOnPoint_%s'%sp
            node = self.gnodes.get(name)
            node.inputs['Scale'].default_value = scale[sp]
    
    @property
    def species_props(self):
        return self.get_species_props()
    
    def get_species_props(self):
        species_props = {}
        radius = self.radius
        color = self.color
        for sp in self.elements:
            species_props[sp] = {'radius': radius[sp], 'color': color[sp]}
        return species_props
    
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
        species = {}
        collection = self.obj.batoms.batom.species
        i = 0
        for data in collection:
            species[data.name] = {'index': i, 'elements': {}}
            elecollection = data.elements
            for eledata in elecollection:
                species[data.name]['elements'][eledata.name] = round(eledata.occupancy, 3)
            i += 1
        return species
    
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
        species_props = get_default_species_data(elements)
        self.build_materials(species_props)
        self.assign_materials()

    @property
    def main_elements(self):
        main_elements = {}
        for sp in self.elements:
            sorted_ele = sorted(self.elements[sp]['elements'].items(), key=lambda x: -x[1])
            if sorted_ele[0][0] == 'X':
                main_elements[sp] = sorted_ele[1][0]
            else:
                main_elements[sp] = sorted_ele[0][0]
        return main_elements
    
    @property
    def model_style(self):
        return self.get_model_style()
    
    @model_style.setter
    def model_style(self, model_style):
        self.set_model_style(model_style)
    
    def get_model_style(self):
        return int(self.obj.batoms.model_style)
    
    def set_model_style(self, model_style):
        self.obj.batoms.model_style = str(model_style)
        self.draw(draw_isosurface = False)
    
    @property
    def polyhedra_style(self):
        return int(self.obj.batoms.polyhedra_style)
    
    @polyhedra_style.setter
    def polyhedra_style(self, polyhedra_style):
        self.obj.batoms.polyhedra_style = str(polyhedra_style)
        self.draw()
    
    @property
    def show_unit_cell(self):
        return self.obj.batoms.show_unit_cell
    
    @show_unit_cell.setter
    def show_unit_cell(self, show_unit_cell):
        self.obj.batoms.show_unit_cell = show_unit_cell
        self.draw_cell()

    @property
    def radius(self):
        return self.get_radius()
    
    def get_radius(self):
        radius = {}
        for sp in self.elements:
            radius[sp] = self.instancers[sp].batoms.batom.radius
        return radius
    
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
        species_props = get_default_species_data(self.elements,
                                radius_style = radius_style)
        # print(species_props)
        instancer = self.build_instancer(radius = species_props['radius'], 
                                scale = scale)
        instancer.parent = self.obj

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
        return radii

    @property
    def size(self):
        return self.get_size()
    
    @size.setter
    def size(self, size):
        self.set_size(size)
    
    def get_size(self):
        size = {}
        radius = self.radius
        scale = self.scale
        for sp in radius:
            size[sp] = radius[sp]*scale[sp][0]
        return size
    
    def set_size(self, size):
        scale = {}
        radius = self.radius
        for sp in self.elements:
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
        me = self.obj.data
        for key, data in attributes.items():
            # print(key)
            att = me.attributes.get(key)
            dtype = type(attributes[key][0])
            if np.issubdtype(dtype, int):
                dtype = 'INT'
            elif np.issubdtype(dtype, float):
                dtype = 'FLOAT'
            elif np.issubdtype(dtype, str):
                dtype = 'STRING'
            else:
                raise KeyError('%s is not supported.'%dtype)
            if att is None:
                att = me.attributes.new(name = key, type = dtype, domain = 'POINT')
            if dtype == 'STRING':
                nvert = len(me.vertices)
                for i in range(nvert):
                    att.data[i].value = data[i]
            else:
                att.data.foreach_set("value", data)
        me.update()
    
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
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        # main elements
        main_elements = self.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        arrays.update({'elements': np.array(elements, dtype='U20')})
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
            for ba in self.batoms.values():
                ba.positions = np.dot(ba.positions(), M)
    
    @property
    def pbc(self):
        return self.get_pbc()
    
    @pbc.setter
    def pbc(self, pbc):
        self.set_pbc(pbc)
    
    def get_pbc(self):
        return list(self.obj.batoms.pbc)
    
    def set_pbc(self, pbc):
        if isinstance(pbc, bool):
            pbc = [pbc]*3
        self.obj.batoms.pbc = pbc
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
        color = {}
        for sp in self.elements:
            Viewpoint_color = self.materials[sp][self.main_elements[sp]].diffuse_color
            for node in self.materials[sp][self.main_elements[sp]].node_tree.nodes:
                if 'Base Color' in node.inputs:
                    node_color = node.inputs['Base Color'].default_value[:]
                if 'Alpha' in node.inputs:
                    Alpha = node.inputs['Alpha'].default_value
            color[sp] = [node_color[0], node_color[1], node_color[2], Alpha]
        return color
    
    def set_color(self, color):
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.materials[self.main_elements].diffuse_color = color
        for node in self.materials[self.main_elements].node_tree.nodes:
            if 'Base Color' in node.inputs:
                node.inputs['Base Color'].default_value = color
            if 'Alpha' in node.inputs:
                node.inputs['Alpha'].default_value = color[3]
    
    @property
    def gnodes(self):
        return self.get_gnodes()
    
    @gnodes.setter
    def gnodes(self, gnodes):
        self.set_gnodes(gnodes)
    
    def get_gnodes(self):
        print(self.obj.modifiers[0])
        gnodes = self.obj.modifiers.get("GeometryNodes_%s"%self.label)
        return gnodes.node_group.nodes
    
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
        if self.isosurfacesetting.volume is not None:
            self.isosurfacesetting.volume = np.tile(self.isosurfacesetting.volume, m)
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
        #
        for sp, instancer in self.instancers.items():
            instancer = instancer.copy()
            bpy.data.collections['Collection'].objects.link(instancer)
            instancer.hide_set(True)
            instancer.name = '%s_%s_instancer'%(name, sp)
            batom = self.__class__(name)
        batom.translate([2, 2, 2])
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
    
    def __repr__(self) -> str:
        text = []
        text.append('label={0}, '.format(self.label))
        text.append('species=%s, '%(list(self.elements)))
        text.append('cell={0}, '.format(self.cell))
        text.append('pbc={0}'.format(self.pbc))
        # text.append('positions={0}'.format(self.positions))
        text = "".join(text)
        text = "Batoms(%s)"%text
        return text
    
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
    def hide(self):
        return self.get_hide()
    
    @hide.setter
    def hide(self, state):
        self.set_hide(state)
    
    def get_hide(self):
        return 'Unknown'
    
    def set_hide(self, state, only_atoms = False):
        names = self.coll.all_objects.keys()
        for name in names:
            obj = bpy.data.objects.get(name)
            if only_atoms and not obj.batoms.batom.flag: continue
            obj.hide_render = state
            obj.hide_set(state)
    
    def get_spacegroup_number(self, symprec = 1e-5):
        """
        """
        try:
            import spglib
        except ImportError:
            return 1
        atoms = self.atoms
        sg = spglib.get_spacegroup((atoms.get_cell(), atoms.get_scaled_positions(),
                                    atoms.get_atomic_numbers()),
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
        positions = self.atoms.positions
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
        for sp, batom in self.batoms.items():
            lock_to(batom.instancer, obj, location = False, rotation = True)
    
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