"""
"""
import bpy
import numpy as np
from time import time
from batoms.utils import string2Number
from batoms.base.collection import Setting, tuple2string
# from pprint import pprint
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

class BondSetting():

    def __init__(self, label, name, bonds=None) -> None:
        self.label = label
        self.coll_name = label
        self.name = name
        self.bonds = bonds
        self.type = 'Bbond'

    @property
    def coll(self):
        return self.get_coll()

    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            coll = bpy.data.collections.new(self.coll_name)
            bpy.context.scene.collection.children.link(coll)
        return coll

    @property
    def bpy_setting(self):
        return self.get_bpy_setting()

    def get_bpy_setting(self):
        bpy_setting = getattr(self.coll, self.type).settings
        return bpy_setting

    @property
    def species1(self):
        return self.bpy_setting[self.name].species1

    @property
    def species2(self):
        return self.bpy_setting[self.name].species2

    @property
    def min(self):
        return self.bpy_setting[self.name].min

    @min.setter
    def min(self, min):
        self.bpy_setting[self.name].min = min

    @property
    def max(self):
        return self.bpy_setting[self.name].max

    @max.setter
    def max(self, max):
        self.bpy_setting[self.name].max = max

    @property
    def search(self):
        return self.bpy_setting[self.name].search

    @search.setter
    def search(self, search):
        self.bpy_setting[self.name].search = search

    @property
    def polyhedra(self):
        return self.bpy_setting[self.name].polyhedra

    @polyhedra.setter
    def polyhedra(self, polyhedra):
        self.bpy_setting[self.name].polyhedra = polyhedra

    @property
    def width(self):
        return self.bpy_setting[self.name].width

    @width.setter
    def width(self, width):
        self.bpy_setting[self.name].width = width

    @property
    def order(self):
        return self.bpy_setting[self.name].order

    @order.setter
    def order(self, order):
        self.bpy_setting[self.name].order = order
        self.bonds.set_attribute_with_indices('order', self.indices, order)
        # if instancer with this order not exist, add one
        sp = self.as_dict()
        self.bonds.settings.build_instancer(sp)
        self.bonds.add_geometry_node(sp)

    @property
    def style(self):
        return self.bpy_setting[self.name].style

    @style.setter
    def style(self, style):
        self.bpy_setting[self.name].style = str(style)
        self.bonds.set_attribute_with_indices('style', self.indices, style)
        # if instancer with this style not exist, add one
        sp = self.as_dict()
        self.bonds.settings.build_instancer(sp)
        self.bonds.add_geometry_node(sp)

    @property
    def color1(self):
        return self.bpy_setting[self.name].color1

    @color1.setter
    def color1(self, color1):
        self.bpy_setting[self.name].color1 = color1
        sp = self.as_dict()
        self.bonds.settings.build_instancer(sp)
        self.bonds.add_geometry_node(sp)

    @property
    def color2(self):
        return self.bpy_setting[self.name].color2

    @color2.setter
    def color2(self, color2):
        self.bpy_setting[self.name].color2 = color2
        sp = self.as_dict()
        self.bonds.settings.build_instancer(sp)
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
        s = 'Bondpair   min     max   width    Search_bond    Polyhedra  Order Style\n'
        s += '{:10s} {:4.3f}   {:4.3f}  {:4.3f}   {:10s}   {:10s}  {}   {}\n'.format(
            self.name, self.min, self.max, self.width, str(
                self.search), str(self.polyhedra),
            self.order, self.style)
        s += '-'*60 + '\n'
        return s

    def as_dict(self) -> dict:
        data = self.bpy_setting[self.name].as_dict()
        return data


class BondSettings(Setting):
    def __init__(self, label, batoms=None,
                 bonds=None,
                 bondsettings=None,
                 dcutoff=0.5) -> None:
        """
        BondSettings object
        The BondSettings object store the bondpair information.

        Parameters:
        label: str
            The label define the batoms object that a Bondsetting belong to.
        dcutoff: float
            extra cutoff used to calculate the maxmium bondlength for bond pairs.
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Bbond'
        self.dcutoff = dcutoff
        self.batoms = batoms
        self.bonds = bonds
        if len(self) == 0:
            self.add_species_list(self.batoms.species.keys(), dcutoff)
        if bondsettings is not None:
            for key, data in bondsettings.items():
                self[key] = data

    def get_collection(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            collection = getattr(coll, self.name)
        elif self.obj_name:
            obj = bpy.data.objects.get(self.obj_name)
            collection = getattr(obj, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return collection.settings


    def build_materials(self, sp, order=None, style = None, node_inputs=None,
                        material_style='default'):
        """
        """
        from batoms.material import create_material
        if not order:
            order = sp["order"]
        if not style:
            style = int(sp["style"])
        colors = [sp['color1'], sp['color2']]
        for i in range(2):
            name = 'bond_%s_%s_%s_%s_%s' % (
                self.label, sp['name'], order, style, i)
            if name in bpy.data.materials:
                mat = bpy.data.materials.get(name)
                bpy.data.materials.remove(mat, do_unlink=True)
            create_material(name,
                            color=colors[i],
                            node_inputs=node_inputs,
                            material_style=material_style,
                            backface_culling=True)

    def build_instancer(self, sp, order=None, style=None, vertices=32, shade_smooth=True):
        """_summary_

        Args:
            sp (_type_): _description_
            vertices (int, optional): _description_. Defaults to 32.
            shade_smooth (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        from batoms.utils.butils import get_nodes_by_name

        # only build the needed one
        if not order:
            order = sp["order"]
        if not style:
            style = int(sp["style"])
        name = 'bond_%s_%s_%s_%s' % (self.label, sp['name'], order, style)
        radius = sp['width']
        self.delete_obj(name)
        if style == 3:
            obj = self.build_spring_bond(name, radius=radius)
        else:
            obj = self.cylinder(name, order=order, style=style,
                                vertices=vertices, depth=1, radius=radius)
        obj.Bbond.width = sp['width']
        obj.batoms.type = 'INSTANCER'
        #
        for coll in obj.users_collection:
            coll.objects.unlink(obj)
        self.batoms.coll.children['%s_instancer' %
                                  self.label].objects.link(obj)
        # bpy.context.scene.objects.unlink(bb.obj)
        if shade_smooth:
            bpy.ops.object.shade_smooth()
        obj.hide_set(True)
        obj.hide_render = True
        obj.scale = [0.001, 0.001, 0.001]

        #
        self.build_materials(sp, order, style, material_style=sp['material_style'])
        self.assign_materials(sp, order, style)
        # update geometry nodes
        ObjectInstancer = get_nodes_by_name(self.bonds.gnodes.node_group.nodes,
                                            'ObjectInfo_%s' % name,
                                            'GeometryNodeObjectInfo')
        if ObjectInstancer is not None:
            ObjectInstancer.inputs['Object'].default_value = obj
        # obj.data.materials.append(materials[data[0]])
        bpy.context.view_layer.update()
        return obj

    def cylinder(self, name, order=1, style=1, vertices=32, depth=1.0, radius=1.0):
        """
        create cylinder and subdivde to 2 parts
        todo: subdivde by the a ratio = radius1/radius2
        """
        style = int(style)
        radius = radius/order
        if style == 2:
            depth = depth/20
            radius = radius
        bpy.ops.mesh.primitive_cylinder_add(
            vertices=vertices, depth=depth, radius=radius)
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
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        me.edges.foreach_set('select', selects)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.subdivide(
            number_cuts=1, smoothness=0, fractal_along_normal=0)
        bpy.ops.object.mode_set(mode='OBJECT')
        # order
        mod = obj.modifiers.new('arrayOrder', 'ARRAY')
        mod.relative_offset_displace = (2, 0, 0)
        mod.count = order
        nv = len(obj.data.vertices)
        for i in range(nv):
            obj.data.vertices[i].co[0] -= 2*(order - 1)*radius
        # bpy.ops.object.modifier_apply(modifier='Array')
        #
        if style == 2:
            mod = obj.modifiers.new('arrayStyle', 'ARRAY')
            mod.relative_offset_displace = (0, 0, 2)
            mod.count = 10
            nv = len(obj.data.vertices)
            for i in range(nv):
                obj.data.vertices[i].co[2] -= (10 - 1)*depth
        # bpy.ops.object.modifier_apply(modifier='Array')
        obj.name = name
        obj.data.name = name
        return obj

    def build_spring_bond(self, name, radius=1.0, resolution=5, n=10,
                          bevel_radius=0.02):
        """Build a spring bond instancer

        Args:
            name (str): name of bond
            radius (float, optional):
                radius of the spiral (bond width). Defaults to 1.0.
            resolution (int, optional):
                number of point in one rotation of the spiral.
                Defaults to 5.
            n (int, optional):
                number of rotations. Defaults to 10.
            bevel_radius (float, optional):
                radius of the bevel object. Defaults to 0.02.

        Returns:
            bpy.data.object: the instancer.
        """

        # Prepare arrays x, y, z
        theta = np.linspace(-n * np.pi, n * np.pi, resolution*n)
        z = np.linspace(-0.5, 0.5, resolution*n)
        r = radius*1.5
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        vertices = np.concatenate((x.reshape(-1, 1),
                                   y.reshape(-1, 1),
                                   z.reshape(-1, 1)), axis=1)
        #
        crv = bpy.data.curves.new('spring', 'CURVE')
        crv.dimensions = '3D'
        crv.resolution_u = 3
        crv.fill_mode = 'FULL'
        crv.use_fill_caps = True
        crv.twist_mode = 'Z_UP'  # ''
        spline = crv.splines.new(type='NURBS')
        nvert = len(vertices)
        spline.points.add(nvert-1)
        vertices = np.append(vertices, np.ones((nvert, 1)), axis=1)
        vertices = vertices.reshape(-1, 1)
        spline.points.foreach_set('co', vertices)
        # bevel
        bpy.ops.curve.primitive_bezier_circle_add(
            radius=bevel_radius, enter_editmode=False)
        bevel_control = bpy.context.active_object
        bevel_control.data.resolution_u = 2
        bevel_control.hide_set(True)
        bevel_control.hide_render = True
        bevel_control.data.name = bevel_control.name = '%s_bevel' % \
            name
        crv.bevel_object = bevel_control
        crv.bevel_mode = 'OBJECT'
        obj = bpy.data.objects.new('spring', crv)
        bpy.data.collections['Collection'].objects.link(obj)
        depsgraph = bpy.context.evaluated_depsgraph_get()
        object_eval = obj.evaluated_get(depsgraph)
        mesh = bpy.data.meshes.new_from_object(object_eval)
        obj_m = bpy.data.objects.new(name, mesh)
        self.delete_obj('spring')
        self.delete_obj('%s_bevel' % name)
        return obj_m

    def assign_materials(self, sp, order, style=0):
        # sort element by occu
        style = int(style)
        me = self.instancers[sp["name"]]['%s_%s' % (order, style)].data
        me.materials.clear()
        if style in [0, 2, 3]:
            for i in range(2):
                me.materials.append(
                    self.materials[sp["name"]]['%s_%s_%s' % (order, style, 0)])
        elif style == 1:
            for i in range(2):
                me.materials.append(
                    self.materials[sp["name"]]['%s_%s_%s' % (order, style, i)])
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
                for style in [0, 1, 2, 3]:
                    name = 'bond_%s_%s_%s_%s' % (
                        self.label, data['name'], order, style)
                    instancers[data['name']]['%s_%s' %
                                             (order, style)] = \
                        bpy.data.objects.get(name)
            data = sp.as_dict(reversed=True)
            if sp.search == 2:
                instancers[data['name']] = {}
                for order in [1, 2, 3]:
                    for style in [0, 1, 2, 3]:
                        name = 'bond_%s_%s_%s_%s' % (
                            self.label, data['name'], order, style)
                        instancers[data['name']]['%s_%s' %
                                                 (order, style)] = \
                            bpy.data.objects.get(name)
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
                for style in [0, 1, 2, 3]:
                    for i in range(2):
                        name = 'bond_%s_%s_%s_%s_%s' % (
                            self.label, data['name'], order, style, i)
                        mat = bpy.data.materials.get(name)
                        materials[data['name']]['%s_%s_%s' %
                                                (order, style, i)] = mat
            data = sp.as_dict(reversed=True)
            if sp.search == 2:
                materials[data['name']] = {}
                for order in [1, 2, 3]:
                    for style in [0, 1, 2, 3]:
                        for i in range(2):
                            name = 'bond_%s_%s_%s_%s_%s' % (
                                self.label, data['name'], order, style, i)
                            mat = bpy.data.materials.get(name)
                            materials[data['name']]['%s_%s_%s' %
                                                    (order, style, i)] = mat
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
            subset = self.bpy_setting.add()
            self.ui_list_index = len(self) - 1
        subset.name = name
        subset.species1 = index[0]
        subset.species2 = index[1]
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_instancer(subset.as_dict())
        if subset.search == 2 and subset.species1 != subset.species2:
            self.build_instancer(subset.as_dict(reversed=True))

    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting' % (name, self.name))
        item = BondSetting(label=self.label, name=name, bonds=self.bonds)
        return item

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
                pair = (sp1, sp2)
                props = {sp1: self.batoms.species.species_props[sp1],
                         sp2: other.batoms.species.species_props[sp2]}
                self.add(pair, props)

    def add_species_list(self, speciesList,
                         only_default=True):
        for sp in speciesList:
            self.add_species(sp, only_default)

    def add_species(self, name,
                         only_default=True):
        """_summary_

        Args:
            speciesList (_type_): _description_
            self_interaction (bool, optional): _description_. Defaults to True.
            only_default (bool, optional): _description_. Defaults to False.
        """
        from batoms.data import default_bonds
        species_props = self.batoms.species.species_props
        pairs = []
        ele = species_props[name]['element']
        for name1, props in species_props.items():
            ele1 = props['element']
            if only_default:
                pair = (ele, ele1)
                if pair in default_bonds:
                    pairs.append((name, name1))
                pair = (ele1, ele)
                if pair in default_bonds:
                    pairs.append((name1, name))
            else:
                pairs.append((name, name1))
                pairs.append((name1, name))
        for pair in pairs:
            self.add(pair)

    def add(self, key, value = None):
        """_summary_

        Args:
            key (str): _description_
        """
        species_props0 = self.batoms.species.species_props
        species_props = {}
        for sp in key:
            if sp in self.batoms.species.keys():
                species_props[sp] = species_props0[sp]
        if value:
            species_props.update(value)
        bond = self.get_bondtable(key, species_props, dcutoff=0.5)
        self[key] = bond

    def copy(self, label):
        from batoms.utils.butils import object_mode
        object_mode()
        bondsetting = self.__class__(label)
        for b in self:
            bondsetting[(b.species1, b.species2)] = b.as_dict()
        return bondsetting

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Bondpair    min     max   Search_bond    Polyhedra style\n"
        for b in self.bpy_setting:
            s += "{:10s}  {:4.3f}   {:4.3f}      ".format(b.name, b.min, b.max)
            s += "{:10s}   {:10s}  {:4s} \n".format(str(b.search),
                                                    str(b.polyhedra), b.style)
        s += "-"*60 + "\n"
        return s

    def as_dict(self) -> dict:
        bondsetting = {}
        for b in self.bpy_setting:
            bondsetting[(b.species1, b.species2)] = b.as_dict()
        return bondsetting

    @property
    def cutoff_dict(self):
        cutoff_dict = {}
        for b in self.bpy_setting:
            cutoff_dict[(b.species1, b.species2)] = [b.min, b.max]
            # if b.search == 2:
            # cutoff_dict[(b.species2, b.species1)] = [b.min, b.max]
        return cutoff_dict

    @property
    def search_dict(self):
        search_dict = {}
        for b in self.bpy_setting:
            search_dict[(b.species1, b.species2)] = b.search
            # if b.search == 2:
            # search_dict[(b.species2, b.species1)] = b.search
        return search_dict

    @property
    def polyhedra_dict(self):
        polyhedra_dict = {}
        for b in self.bpy_setting:
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
                temp = bondlists0[
                    (speciesarray[bondlists0[:, 0]] == b.species1)
                    & (speciesarray[bondlists0[:, 1]] == b.species2)]
                bondlist1.extend(temp)
                temp = offsets_skin0[(
                    speciesarray[offsets_skin0[:, 0]] == b.species1)]
                offsets_skin1.extend(temp)
            elif b.search == 2:
                #  recursively if either sp1 or sp2
                temp = bondlists0[(
                    (speciesarray[bondlists0[:, 0]] == b.species1)
                    & (speciesarray[bondlists0[:, 1]] == b.species2))]
                bondlist2.extend(temp)
                temp1 = temp.copy()
                temp1[:, 0] = temp[:, 1]
                temp1[:, 1] = temp[:, 0]
                temp1[:, 2:] = -temp[:, 2:]
                bondlist2.extend(temp1)
                temp = offsets_skin0[(
                    speciesarray[offsets_skin0[:, 0]] == b.species1)
                    | (speciesarray[offsets_skin0[:, 0]] == b.species2)]
                offsets_skin2.extend(temp)
        return np.array(offsets_skin1), np.array(bondlist1), \
            np.array(offsets_skin2), np.array(bondlist2)

    def get_bondtable(self, pair, props, dcutoff=0.5):
        """
        """
        from batoms.data import default_bonds

        sp1 = pair[0]
        sp2 = pair[1]
        bondmax = props[sp1]['radius'] + \
            props[sp2]['radius'] + dcutoff
        pair_ele = (props[sp1]['element'], props[sp2]['element'])
        if pair_ele in default_bonds:
            search = default_bonds[pair_ele][0]
            polyhedra = default_bonds[pair_ele][1]
            bondtype = default_bonds[pair_ele][2]
        else:
            search = 1
            polyhedra = 0
            bondtype = 0
        bond = {
            'label': self.label,
            'species1': sp1,
            'species2': sp2,
            'min': 0.0,
            'max': bondmax,
            'search': search,
            'polyhedra': polyhedra,
            'color1': props[sp1]['color'],
            'color2': props[sp2]['color'],
            'order': 1,
            'order_offset': 0.15,
            'style': '1',
            'type': bondtype,
        }
        # special for hydrogen bond
        if bondtype == 1:
            bond['min'] = 1.5
            bond['max'] = 2.5
            bond['search'] = 0
            bond['color1'] = [0.1, 0.1, 0.1, 1.0]
            bond['style'] = '2'
        return bond


if __name__ == "__main__":
    pass
