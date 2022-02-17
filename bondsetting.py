"""
"""
from re import I
import bpy
from batoms.butils import object_mode
import numpy as np
from time import time
from batoms.tools import string2Number
from batoms.base import Setting, tuple2string
from pprint import pprint


class BondSetting():

    def __init__(self, label, name, bonds = None) -> None:
        self.label = label
        self.coll_name = label
        self.name = name
        self.bonds = bonds
        self.type = 'bbond'
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            coll = bpy.data.collections.new(self.coll_name)
            bpy.data.scenes['Scene'].collection.children.link(coll)
        return coll

    @property
    def collection(self):
        return self.get_collection()
    
    def get_collection(self):
        collection = getattr(self.coll.batoms, self.type)
        return collection
    
    @property
    def species1(self):
        return self.collection[self.name].species1
    
    @property
    def species2(self):
        return self.collection[self.name].species2

    @property
    def min(self):
        return self.collection[self.name].min
    
    @property
    def max(self):
        return self.collection[self.name].max
    
    @property
    def search(self):
        return self.collection[self.name].search
    
    @property
    def polyhedra(self):
        return self.collection[self.name].polyhedra
    
    @search.setter
    def search(self, search):
        self.collection[self.name].search = search
        # if search not exist, add one

    @property
    def order(self):
        return self.collection[self.name].order
    
    @order.setter
    def order(self, order):
        order0 = self.bonds.attributes['order']
        order0[self.indices] = order
        self.bonds.set_attributes({'order': order0})
        # if order not exist, add one
        self.collection[self.name].order = order
        sp = self.as_dict()
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)
    
    @property
    def style(self):
        return self.collection[self.name].style
    
    @style.setter
    def style(self, style):
        self.bonds.obj.data.attributes['style'].data[self.indices].value = style
        # if stule not exist, add one
        sp = self.species
        sp['style'] = style
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)
    
    @property
    def color1(self):
        return self.collection[self.name].color1
    
    @color1.setter
    def color1(self, color1):
        self.collection[self.name].color1 = color1
        sp = self.as_dict()
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)
    
    @property
    def color2(self):
        return self.collection[self.name].color2
    
    @color2.setter
    def color2(self, color2):
        self.collection[self.name].color2 = color2
        sp = self.as_dict()
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)
    
    @property
    def style(self):
        return self.collection[self.name].style
    
    @style.setter
    def style(self, style):
        self.bonds.obj.data.attributes['style'].data[self.indices].value = style
        # if stule not exist, add one
        sp = self.species
        sp['style'] = style
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)

    @property
    def indices(self):
        return self.get_indices()
    
    def get_indices(self):
        sp1 = self.bonds.arrays['species_index1']
        sp2 = self.bonds.arrays['species_index2']
        indices = np.where((sp1 == string2Number(self.species1)) & 
                           (sp2 == string2Number(self.species2)))[0]
        return indices

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair   min     max   Search_bond    Polyhedra  Order Style\n'
        s += '{:10s} {:4.3f}   {:4.3f}    {:10s}   {:10s}  {}   {}\n'.format(\
                self.name, self.min, self.max, str(self.search), str(self.polyhedra),
                self.order, self.style)
        s += '-'*60 + '\n'
        return s

    def as_dict(self) -> dict:
        data = self.collection[self.name].as_dict()
        return data


class BondSettings(Setting):
    """
    BondSetting object

    The BondSetting object store the bondpair information.

    Parameters:

    label: str
        The label define the batoms object that a Bondsetting belong to.
    cutoff: float
        Cutoff used to calculate the maxmium bondlength for bond pairs.
    """
    
    def __init__(self, label, batoms = None, 
                bonds = None,
                bondsetting = None, 
                cutoff = 1.3) -> None:
        Setting.__init__(self, label, coll_name='%s'%label)
        self.label = label
        self.name = 'bbond'
        self.cutoff = cutoff
        self.batoms = batoms
        self.bonds = bonds
        if len(self) == 0:
            self.set_default(self.batoms.species.species_props, cutoff)
        if bondsetting is not None:
            for key, data in bondsetting.items():
                self[key] = data
        #

    def build_materials(self, sp, node_inputs = None, 
                material_style = 'default'):
        """
        """
        from batoms.material import create_material
        colors = [sp['color1'], sp['color2']]
        for i in range(2):
            name = '%s_%s_%s_%s_%s'%(self.label, sp['name'], sp['order'], sp['style'], i)
            if name in bpy.data.materials:
                mat = bpy.data.materials.get(name)
                bpy.data.materials.remove(mat, do_unlink = True)
            create_material(name,
                        color = colors[i],
                        node_inputs = node_inputs,
                        material_style = material_style,
                        backface_culling = True)
    
    def build_instancer(self, sp, vertices = 32, shade_smooth = True):
        # only build the needed one
        order = sp['order']
        style = int(sp['style'])
        name = 'bond_%s_%s_%s_%s'%(self.label, sp['name'], order, style)
        radius = sp['width']
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        self.cylinder(order = order, style = style, vertices = vertices, depth = 1, radius = radius)
        obj = bpy.context.view_layer.objects.active
        obj.name = name
        obj.data.name = name
        obj.batoms.batom.radius = radius
        #
        obj.users_collection[0].objects.unlink(obj)
        self.batoms.coll.children['%s_instancer'%self.label].objects.link(obj)
        # bpy.data.scenes['Scene'].objects.unlink(bb.obj)
        if shade_smooth:
            bpy.ops.object.shade_smooth()
        obj.hide_set(True)
        obj.hide_render = True
        #
        self.build_materials(sp)
        self.assign_materials(sp, order, style)
        # obj.data.materials.append(materials[data[0]])
        bpy.context.view_layer.update()
        return obj
    
    def cylinder(self, order = 1, style = 1, vertices = 32, depth = 1.0, radius = 1.0):
        """
        create cylinder and subdivde to 2 parts
        todo: subdivde by the a ratio = radius1/radius2
        """
        style = int(style)
        radius = radius/order
        if style == 2:
            depth = depth/20
            radius = radius
        bpy.ops.mesh.primitive_cylinder_add(vertices = vertices, depth = depth, radius = radius)
        obj = bpy.context.view_layer.objects.active
        me = obj.data
        # select edges for subdivde
        n = len(me.vertices)
        vertices = np.zeros(n*3, dtype=np.float64)
        me.vertices.foreach_get('co', vertices)  
        vertices = vertices.reshape((n, 3))
        #
        me.update()
        m = len(me.edges)
        selects = np.zeros(m, dtype=bool)
        for i in range(m):
            v0 = me.edges[i].vertices[0]
            v1 = me.edges[i].vertices[1]
            center = (vertices[v0] + vertices[v1])/2
            if np.isclose(center[2], 0):
                selects[i] = True
        #
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
        me.edges.foreach_set('select', selects)
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.subdivide(number_cuts = 1, smoothness = 0, fractal_along_normal = 0)
        bpy.ops.object.mode_set(mode = 'OBJECT')
        # order 
        mod = obj.modifiers.new('arrayOrder', 'ARRAY')
        mod.relative_offset_displace = (2, 0, 0)
        mod.count = order
        nv = len(obj.data.vertices)
        for i in range(nv):
            obj.data.vertices[i].co[0] -= 2*(order - 1)*radius
        bpy.ops.object.modifier_apply(modifier='Array')
        #
        if style == 2:
            mod = obj.modifiers.new('arrayStyle', 'ARRAY')
            mod.relative_offset_displace = (0, 0, 2)
            mod.count = 10
            nv = len(obj.data.vertices)
            for i in range(nv):
                obj.data.vertices[i].co[2] -= (10 - 1)*depth
        bpy.ops.object.modifier_apply(modifier='Array')
        return obj

    def assign_materials(self, sp, order, style = 0):
        # sort element by occu
        style = int(style)
        me = self.instancers[sp["name"]]['%s_%s'%(order, style)].data
        me.materials.clear()
        if style in [0, 2]:
            for i in range(2):
                me.materials.append(self.materials[sp["name"]]['%s_%s_%s'%(order, style, 0)])
        elif style == 1:
            for i in range(2):
                me.materials.append(self.materials[sp["name"]]['%s_%s_%s'%(order, style, i)])
        # find the face index for sp1 and sp2
        #
        npoly = len(me.polygons)
        centers = np.zeros(npoly*3, dtype=np.float64)
        me.polygons.foreach_get('center', centers)  
        centers = centers.reshape((npoly, 3))
        material_indexs = np.zeros(npoly, dtype='int')
        index = np.where(centers[:, 2] < 0)[0]
        material_indexs[index] = 1
        me.polygons.foreach_set('material_index', material_indexs)

    @property
    def instancers(self):
        return self.get_instancers()
    
    def get_instancers(self):
        instancers = {}
        for sp in self:
            data = sp.as_dict()
            instancers[data['name']] = {}
            for order in [1, 2, 3]:
                for style in [0, 1, 2]:
                    name = 'bond_%s_%s_%s_%s'%(self.label, data['name'], order, style)
                    instancers[data['name']]['%s_%s'%(order, style)] = bpy.data.objects.get(name)
            data = sp.as_dict(reversed = True)
            if sp.search == 2:
                instancers[data['name']] = {}
                for order in [1, 2, 3]:
                    for style in [0, 1, 2]:
                        name = 'bond_%s_%s_%s_%s'%(self.label, data['name'], order, style)
                        instancers[data['name']]['%s_%s'%(order, style)] = bpy.data.objects.get(name)
        return instancers
    
    @property
    def materials(self):
        return self.get_materials()
    
    def get_materials(self):
        materials = {}
        for sp in self:
            data = sp.as_dict()
            materials[data['name']] = {}
            for order in [1, 2, 3]:
                for style in [0, 1, 2]:
                    for i in range(2):
                        name = '%s_%s_%s_%s_%s'%(self.label, data['name'], order, style, i)
                        mat = bpy.data.materials.get(name)
                        materials[data['name']]['%s_%s_%s'%(order, style, i)] = mat
            data = sp.as_dict(reversed = True)
            if sp.search == 2:
                materials[data['name']] = {}
                for order in [1, 2, 3]:
                    for style in [0, 1, 2]:
                        for i in range(2):
                            name = '%s_%s_%s_%s_%s'%(self.label, data['name'], order, style, i)
                            mat = bpy.data.materials.get(name)
                            materials[data['name']]['%s_%s_%s'%(order, style, i)] = mat
        return materials

    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        if isinstance(index, str):
            raise Exception("Bond index should be a tuple or list, \
    e.g. h2o.bondseeting[('O', 'H')]")
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        subset.species1 = index[0]
        subset.species2 = index[1]
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_instancer(subset.as_dict())
        if subset.search == 2 and subset.species1 != subset.species2:
            self.build_instancer(subset.as_dict(reversed = True))

    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting'%(name, self.name))
        item = BondSetting(label = self.label, name = name, bonds = self.bonds)
        return item
        
    def set_default(self, species_props, cutoff = 1.3, self_interaction = True):
        """
        """
        for sel, data in species_props.items():
            bondtable = get_bondtable(self.label, sel, data, cutoff=cutoff, 
                        self_interaction = self_interaction)
            for key, value in bondtable.items():
                self[key] = value
    
    def extend(self, other):
        for b in other:
            self[(b.species1, b.species2)] = b.as_dict()
        # new
        species1 = set(self.batoms.species.species_props)
        species2 = set(other.batoms.species.species_props)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                species = {sp1: self.batoms.species.species_props[sp1], sp2: other.batoms.species.species_props[sp2]}
                self.set_default(species, self.cutoff, self_interaction = False)
    
    def add(self, bondpairs):
        if isinstance(bondpairs, tuple):
            bondpairs = [bondpairs]
        if isinstance(bondpairs[0], (str, int, float)):
            bondpairs = [bondpairs]
        for bondpair in bondpairs:
            species = {sp: self.batoms.species.species_props['sel0'][sp] for sp in bondpair}
            maxlength = species[bondpair[0]]['radius'] + \
                        species[bondpair[1]]['radius']
            self[bondpair] = {'max': maxlength*self.cutoff, 
                            'color1':species[bondpair[0]]['color'],
                            'color2':species[bondpair[1]]['color'],
                            }
            self.set_default(species)
    
    def copy(self, label):
        object_mode()
        bondsetting = self.__class__(label)
        for b in self:
            bondsetting[(b.species1, b.species2)] = b.as_dict()
        return bondsetting
    
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair    min     max   Search_bond    Polyhedra style\n'
        for b in self.collection:
            s += '{:10s}  {:4.3f}   {:4.3f}      {:10s}   {:10s}  {:4s} \n'.format(\
                b.name, b.min, b.max, str(b.search), str(b.polyhedra), b.style)
        s += '-'*60 + '\n'
        return s
    
    def as_dict(self) -> dict:
        bondsetting = {}
        for b in self.collection:
            bondsetting[(b.species1, b.species2)] = b.as_dict()
        return bondsetting

    @property    
    def cutoff_dict(self):
        cutoff_dict = {}
        for b in self.collection:
            cutoff_dict[(b.species1, b.species2)] = [b.min, b.max]
            # if b.search == 2:
                # cutoff_dict[(b.species2, b.species1)] = [b.min, b.max]
        return cutoff_dict
    
    @property    
    def search_dict(self):
        search_dict = {}
        for b in self.collection:
            search_dict[(b.species1, b.species2)] = b.search
            # if b.search == 2:
                # search_dict[(b.species2, b.species1)] = b.search
        return search_dict
    
    @property    
    def polyhedra_dict(self):
        polyhedra_dict = {}
        for b in self.collection:
            polyhedra_dict[(b.species1, b.species2)] = b.polyhedra
            # if b.polyhedra == 2:
                # polyhedra_dict[(b.species2, b.species1)] = b.polyhedra
        return polyhedra_dict
    
    @property    
    def maxcutoff(self):
        maxcutoff = 2
        for bl in self.cutoff_dict.values():
            if bl[1] > maxcutoff:
                maxcutoff = bl[1]
        return maxcutoff
    
    def search_bond_list(self, atoms, bondlists0, offsets_skin0):
        bondlist1 = []
        offsets_skin1 = []
        bondlist2 = []
        offsets_skin2 = []
        speciesarray = np.array(atoms.arrays['species'])
        for b in self:
            if b.search == 1:
                temp = bondlists0[(speciesarray[bondlists0[:, 0]] == b.species1) 
                              & (speciesarray[bondlists0[:, 1]] == b.species2)]
                bondlist1.extend(temp)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.species1)]
                offsets_skin1.extend(temp)
            elif b.search == 2:
                #  recursively if either sp1 or sp2
                temp = bondlists0[((speciesarray[bondlists0[:, 0]] == b.species1) 
                                  & (speciesarray[bondlists0[:, 1]] == b.species2))]
                bondlist2.extend(temp)
                temp1 = temp.copy()
                temp1[:, 0] = temp[:, 1]
                temp1[:, 1] = temp[:, 0]
                temp1[:, 2:] = -temp[:, 2:]
                bondlist2.extend(temp1)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.species1)
                                   | (speciesarray[offsets_skin0[:, 0]] == b.species2)]
                offsets_skin2.extend(temp)
        return np.array(offsets_skin1), np.array(bondlist1), np.array(offsets_skin2), np.array(bondlist2)

def get_bondtable(label, select, speciesdict, cutoff = 1.3, self_interaction = True):
    """
    """
    from batoms.data import default_bonds
    bondtable = {}
    for species1 in speciesdict:
        element1 = species1.split('_')[0]
        color1 = speciesdict[species1]['color']
        radius1 = cutoff * speciesdict[species1]['radius']
        for species2 in speciesdict:
            if not self_interaction and species1 == species2: continue
            element2 = species2.split('_')[0]
            element_pair = (element1, element2)
            if element_pair not in default_bonds: continue
            pair12 = (species1, species2)
            pair21 = (species2, species1)
            # only add bond once, except for hedrogen bond
            if pair21 in bondtable and ((element1 != 'H') and (element2 !='H')): continue
            color2 = speciesdict[species2]['color']
            radius2 = cutoff * speciesdict[species2]['radius']
            bondmax = radius1 + radius2
            bondtable[pair12] = {
                    'flag': True,
                    'label': label,
                    'select': select,
                    'species1': species1, 
                    'species2': species2, 
                    'min': 0.0, 
                    'max': bondmax,
                    'search': default_bonds[element_pair][0], 
                    'polyhedra': default_bonds[element_pair][1],
                    'color1': color1,
                    'color2': color2,
                    'width': 0.10,
                    'order': 1,
                    'order_offset': 0.15,
                    'style': '1',
                    }
    # special for hydrogen bond
    hbs = [('H', 'O'), ('H', 'N')]
    for pair in bondtable:
        key = (pair[0].split('_')[0], pair[1].split('_')[0])
        if key in hbs:
            bondtable[pair]['min'] = 1.2
            bondtable[pair]['max'] = 2.1
            bondtable[pair]['search'] = 0
            bondtable[pair]['color1'] = [0.1, 0.1, 0.1, 1.0]
            bondtable[pair]['width'] = 0.03
            bondtable[pair]['style'] = '2'
    return bondtable

if __name__ == "__main__":
    pass