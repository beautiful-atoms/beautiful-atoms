from numpy.lib import select
import bpy
import numpy as np
from batoms.base import Setting
from batoms.tools import get_default_species_data, string2Number, number2String


class species():
    def __init__(self, name, species, elements, 
                radius_style = 'colvent', 
                color_style = 'ASE',
                
                ) -> None:
        self.species = species
        self.elements = elements
        self.name = name
    
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
        self._species.build_instancers(segments = segments)
        self.build_geometry_node()
        

class Bspecies(Setting):
    """Bspecies Class
    """
    def __init__(self, label, coll_name, species = {}, 
                batoms = None,
                material_style = 'default',
                segments = [16, 16],
                subdivisions = 2,
                ) -> None:
        Setting.__init__(self, label, coll_name=coll_name)
        self.label = label
        self.name = 'bspecies'
        self.batoms = batoms
        self.segments = segments
        self.subdivisions = subdivisions
        for sp, data in species.items():
            self[sp] = data
    
    def __getitem__(self, index):
        return self.species[index]

    def __setitem__(self, name, data):
        """
        Set species
        """
        self.add(name, data)
        
    def add(self, name, data):
        """
        """
        if 'color_style' not in data: 
            data['color_style'] = '0'
        if 'radius_style' not in data: 
            data['radius_style'] = '0'
        sp = self.find(name)
        if sp is None:
            sp = self.collection.add()
        sp.name = name
        sp.label = self.label
        # add default props
        data['elements'] = self.check_elements(data['elements'])
        props = get_default_species_data(data['elements'],
                                radius_style = data['radius_style'], 
                                color_style = data['color_style'])
        props.update(data)
        #
        for key, value in props.items():
            if key == 'elements': continue
            setattr(sp, key, value)
        # add elements
        sp.elements.clear()
        for ele, occupancy in props['elements'].items():
            eledata = sp.elements.add()
            eledata.name = ele
            eledata.occupancy = occupancy
            eledata.color = props['color'][ele]
        for sel in self.batoms.selects:
            self.build_instancer(sp, select = sel.name)
            self.batoms.add_geometry_node(sp.name, sel.name)

    def keys(self):
        return self.species.keys()

    def items(self):
        return self.species.items()

    def __repr__(self):
        s = str(self.species)
        # s = "Bspecies(species = '%s', index = %s, elements = %s, \
                    # radius = %s" % (self.species,  self.index, \
                    # str(self.elements), self.radius)
        return s
    
    @property
    def instancers(self):
        return self.get_instancers()
    
    def get_instancers(self):
        instancers = {}
        for sel in self.batoms.selects:
            instancers[sel.name] = {}
            for sp in self.batoms.species:
                name = '%s_%s_%s'%(self.label, sp.name, sel.name)
                instancers[sel.name][sp.name] = bpy.data.objects.get(name)
        return instancers
    
    @property
    def materials(self):
        return self.get_materials()
    
    def get_materials(self):
        materials = {}
        for sel in self.batoms.selects:
            materials[sel.name] = {}
            for sp, data in self.batoms.species.items():
                materials[sel.name][sp] = {}
                for ele in data['elements']:
                    name = '%s_%s_%s_%s'%(self.label, sp, ele, sel.name)
                    mat = bpy.data.materials.get(name)
                    materials[sel.name][sp][ele] = mat
        return materials

    @property
    def species(self):
        return self.get_species()
    
    def get_species(self):
        species = {}
        for sp in self.collection:
            species[sp.name] = {'radius': sp.radius, 'elements': {}}
            elecollection = sp.elements
            for eledata in elecollection:
                species[sp.name]['elements'][eledata.name] = round(eledata.occupancy, 3)
        return species
    
    @property
    def index(self):
        return self.batoms.obj.batoms.batom.species.find(self.species)
    
    @property
    def elements(self):
        return self.get_elements()
    
    def get_elements(self):
        elements = {}
        data = self.collection.elements
        for eledata in data:
            elements[eledata.name] = round(eledata.occupancy, 3)
        return elements
    
    @property
    def main_elements(self):
        main_elements = {}
        for sp in self.species:
            sorted_ele = sorted(self.species[sp]['elements'].items(), key=lambda x: -x[1])
            if sorted_ele[0][0] == 'X':
                main_elements[sp] = sorted_ele[1][0]
            else:
                main_elements[sp] = sorted_ele[0][0]
        return main_elements
    
    @property
    def species_props(self):
        return self.get_species_props()
    
    def get_species_props(self):
        species_props = {}
        instancers = self.instancers
        for sel in self.batoms.selects:
            species_props[sel.name] = {}
            for sp in self.species:
                radius = instancers[sel.name][sp].batoms.batom.radius
                matname = instancers[sel.name][sp].material_slots[0].name
                color = bpy.data.materials[matname].node_tree.nodes[0].inputs[0].default_value[:]
                species_props[sel.name][sp] = {'radius': radius, 'color': color}
        return species_props
    
    def __add__(self, other):
        self += other
        return self
    
    def __iadd__(self, other):
        self.extend(other)
        return self
    
    def copy(self, name):
        for sel in self.batoms.selects:
            for sp, instancer in self.instancers[sel.name].items():
                instancer = instancer.copy()
                bpy.data.collections['Collection'].objects.link(instancer)
                instancer.hide_set(True)
                instancer.name = '%s_%s_instancer'%(name, sp)
        bspecies = self.__class__(name, name)
        return bspecies
        
    def extend(self, other):
        """

        """
        species1 = self.species
        species2 = other.species
        species3 = {}
        for key, data in species2.items():
            if key not in species1:
                instancer = other.instancers[key]
                name = '%s_%s_instancer'%(self.label, key)
                instancer.name = name
                # materials
                # for mat in instancer.materials_slots:
                    # mat.name.replace()
                species3[key] = [data, instancer]
            elif data != species1[key]:
                instancer = other.instancers[key]
                spname = '%s_%s'%(key, other.label)
                name = '%s_%s_instancer'%(self.label, spname)
                instancer.name = name
                species3[spname] = [data, instancer]
        #
        for key, data in species3.items():
            self.add(key, data[0], data[1])
        return self
    
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    
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
    
    def build_materials(self, sp, select = 'all', node_inputs = None, 
                material_style = 'default'):
        """
        """
        from batoms.material import create_material
        for element in sp.elements:
            name = '%s_%s_%s_%s'%(self.label, sp.name, element.name, select)
            if name not in bpy.data.materials:
                create_material(name,
                            color = element.color,
                            node_inputs = node_inputs,
                            material_style = material_style,
                            backface_culling = True)
    
    def build_instancer(self, sp, select = 'all', 
                        shape = 'UV_SPHERE', 
                        shade_smooth = True):
        name = '%s_%s_%s'%(self.label, sp.name, select)
        radius = sp.radius
        segments = sp.segments
        natom = len(self.batoms)
        if natom >= 1e3:
            segments = [16, 16]
        if natom >= 1e4:
            segments = [10, 10]
        if natom >= 1e5:
            segments = [8, 8]
        if natom >= 1e6:
            segments = [6, 6]
        if self.segments:
            segments = self.segments
        self.delete_obj(name)
        if shape.upper() == 'UV_SPHERE':
            bpy.ops.mesh.primitive_uv_sphere_add(segments = segments[0], 
                                ring_count = segments[1], 
                                radius = radius)
        if shape.upper() == 'ICO_SPHERE':
            shade_smooth = False
            bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = self.subdivisions, 
                        radius = radius)
        if shape.upper() == 'CUBE':
            bpy.ops.mesh.primitive_cube_add(size = radius)
            shade_smooth = False
        if shape.upper() == 'METABALL':
            bpy.ops.object.metaball_add(type = 'BALL', location = [0, 0, 0])
        obj = bpy.context.view_layer.objects.active
        obj.name = name
        obj.data.name = name
        obj.batoms.batom.radius = radius
        #
        obj.users_collection[0].objects.unlink(obj)
        self.coll.children['%s_instancer'%self.coll_name].objects.link(obj)
        if shade_smooth:
            bpy.ops.object.shade_smooth()
        if shape.upper() != 'METABALL':
            obj.hide_set(True)
            obj.hide_render = True
        #
        self.build_materials(sp, select)
        self.assign_materials(sp, select)
        bpy.context.view_layer.update()
        return obj

    def assign_materials(self, sp, select):
        # sort element by occu
        mesh = self.instancers[select][sp.name].data
        mesh.materials.clear()
        sorted_ele = sorted(self.batoms.species[sp.name]['elements'].items(), key=lambda x: -x[1])
        materials = self.materials[select][sp.name]
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

    def build_instancers(self, sel = None, segments = None):
        if segments:
            self.segments = segments
        if select is not None:
            for sp in self:
                self.build_instancer(sp, select = sel.name)
        else:
            for sel in self.batoms.selects:
                for sp in self:
                    self.build_instancer(sp, select = sel.name)
            
    @property
    def radius_style(self):
        return self.get_radius_style()
    
    @radius_style.setter
    def radius_style(self, radius_style):
        self.set_radius_style(radius_style)
    
    def get_radius_style(self):
        return self.coll.batoms.radius_style
    
    def set_radius_style(self, radius_style):
        # print(species_props)
        for sel in self.batoms.selects:
            for sp in self:
                sp.radius_style = str(radius_style)
                self.build_instancer(sp, select = sel.name)
        self.batoms.build_geometry_node()