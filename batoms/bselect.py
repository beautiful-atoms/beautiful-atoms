"""
#TODO add panel for select
#TODO active radius_style
#TODO merge select and SliceBatoms
"""
import bpy
from batoms.base.collection import Setting
import numpy as np
from batoms.utils import string2Number
# from time import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

class Select():
    def __init__(self, label, parent, name='all',
                 batoms=None, indices=None,
                 ) -> None:
        """_summary_

        Args:
            label (_type_): _description_
            name (str, optional): _description_. Defaults to 'all'.
            batoms (_type_, optional): _description_. Defaults to None.
            indices (_type_, optional): _description_. Defaults to None.

        Examples:

        >>>sel = batoms.selects[key]
        >>>sel.model_style = 1

        """
        self.label = label
        self.parent = parent
        self.name = name
        self.batoms = batoms
        self.coll_name = '%s_%s' % (label, name)
        if indices is not None:
            self.indices = indices

    def get_bpy_setting(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            data = getattr(coll.batoms, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.setting

    @property
    def positions(self):
        return self.batoms.positions[self.indices]

    @property
    def radii_vdw(self):
        return self.batoms.radii_vdw[self.indices]

    @property
    def model_style(self):
        return self.get_model_style()

    @model_style.setter
    def model_style(self, model_style):
        self.set_model_style(model_style)

    def get_model_style(self):
        model_style = self.batoms.arrays['model_style']
        return model_style[self.indices]

    def set_model_style(self, model_style):
        if int(model_style) == 0:
            self.scale = 1
        elif int(model_style) == 1:
            self.scale = 0.4
        elif int(model_style) == 2:
            self.scale = 0.4
        elif int(model_style) == 3:
            self.scale = 0.0001
        model_style0 = self.batoms.attributes['model_style']
        model_style0[self.indices] = model_style
        model_style = {'model_style': model_style0}
        self.batoms.set_attributes(model_style)
        self.batoms.draw()

    @property
    def radius_style(self):
        return self.get_radius_style()

    @radius_style.setter
    def radius_style(self, radius_style):
        self.set_radius_style(radius_style)

    def get_radius_style(self):
        return self.parent.bpy_setting[self.name].radius_style

    def set_radius_style(self, radius_style):
        self.parent.bpy_setting[self.name].radius_style = radius_style

    @property
    def show(self):
        return self.get_show()

    @show.setter
    def show(self, state):
        self.set_show(state)

    def get_show(self):
        return self.batoms.attributes['show'][self.indices]

    def set_show(self, show):
        show0 = self.batoms.show
        show0[self.indices] = show
        self.batoms.set_attributes({'show': show0})

    @property
    def mask(self):
        return self.get_mask()

    def get_mask(self):
        mask = self.batoms.arrays['select']
        mask = np.where(mask == string2Number(self.name))[0]
        return mask

    @property
    def vg(self):
        vg = self.batoms.obj.vertex_groups.get(self.name)
        if vg is None:
            vg = self.batoms.obj.vertex_groups.new(name=self.name)
        return vg

    @property
    def indices(self):
        return self.get_indices()

    @indices.setter
    def indices(self, indices):
        self.set_indices(indices)

    def get_indices(self):
        vg_idx = self.vg.index
        obj = self.batoms.obj
        indices = [v.index for v in obj.data.vertices if vg_idx in [
            vg.group for vg in v.groups]]
        return indices

    def set_indices(self, indices):
        if isinstance(indices, np.ndarray):
            indices = indices.tolist()
        self.vg.add(indices, 1.0, 'ADD')
        select = self.batoms.attributes['select']
        select[self.indices] = string2Number('all')
        select[indices] = string2Number(self.name)
        select = {'select': select}
        self.batoms.set_attributes(select)

    @property
    def scale(self):
        return self.get_scale()

    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)

    def get_scale(self):
        scale = self.batoms.attributes['scale'][self.indices]
        return scale

    def set_scale(self, scale):
        """
        """
        scale0 = self.batoms.attributes['scale']
        scale0[self.indices] = scale
        # for species
        if isinstance(scale, dict):
            species = self.batoms.attributes['species']
            for key, value in scale.items():
                scale0[np.where((species == key) & self.mask)] = value
            scale = scale0
        self.batoms.set_attributes({'scale': scale0})

    @property
    def color(self):
        return self.get_color()

    @color.setter
    def color(self, color):
        self.set_color(color)

    def get_color(self):
        """
        """
        color = {}
        materials = self.materials
        main_elements = self.main_elements
        for sp in self.batoms.species:
            for node in materials[sp][main_elements[sp]].node_tree.nodes:
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
    def radius(self):
        """
        radius depends on radius style, thus could be with select
        """
        radius = {}
        for sp in self.bpy_setting:
            radius[sp.name] = sp.radius
        return radius

    @property
    def species_props(self):
        return self.get_species_props()

    def get_species_props(self):
        species_props = {}
        radius = self.radius
        color = self.color
        for sp in self.species:
            species_props[sp] = {'radius': radius[sp], 'color': color[sp]}
        return species_props

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   model_style   radius_style     show \n'
        s += '{0:10s} {1:10s}   {2:10s}   {3:10s} \n'.format(
            self.name, str(self.model_style[0]),
            self.radius_style[0], str(self.show[0]))
        s += '-'*60 + '\n'
        return s


class Selects(Setting):
    """
    sel = batoms.select.add(name, index)
    sel = batoms.select.from_dict({name: index})
    sel.model_style = 1
    """

    def __init__(self, label, batoms=None) -> None:
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'batoms'
        self.batoms = batoms

    def get_ui_list_index(self):
        return self.bpy_data.ui_list_index_select

    def set_ui_list_index(self, value):
        self.bpy_data.ui_list_index_select = value

    def get_bpy_setting(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            data = getattr(coll, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.settings_select

    @property
    def selects(self):
        return self.get_selects()

    def get_selects(self):
        selects = {}
        for b in self.bpy_setting:
            selects[b.name] = Select(self.label, parent=self, name=b.name,
                                     batoms=self.batoms)
        for vg in self.batoms.obj.vertex_groups:
            if vg.name not in selects:
                self.add(vg.name)
            selects[vg.name] = Select(self.label, parent=self, name=vg.name,
                                      batoms=self.batoms)
        return selects

    def __getitem__(self, key):
        """Return a subset of the Batom.

        key -- str, describing which batom to return.
        """
        if isinstance(key, str):
            if key not in self.selects:
                raise KeyError('%s is not in this structure' % key)
            return self.selects[key]
        elif isinstance(key, list):
            raise KeyError('dict not supported yet!')

    def __len__(self):
        return len(self.selects)

    def __repr__(self) -> str:
        s = '     Name        model_style    show  \n'
        i = 0
        for name, sel in self.selects.items():
            s += '{:3d}   {:10s}  \n'.format(i, name)
            i += 1
        return s

    def add(self, name, expre=None):
        if expre is not None:
            indices = elect_expression(expre, self.batoms)
        else:
            indices = None
        if len(name) > 4:
            print("Warning: the name of the select: {}, longer than four characters. Has been renamed to {}".format(name, name[:4]))
            logger.warning("The name of the select: {}, longer than four characters. Has been renamed to {}".format(name, name[:4]))
            name = name[:4]
        # select = self.batoms.attributes['select']
        # select[indices] = string2Number(name)
        # select = {'select': select}
        # self.batoms.set_attributes(select)
        sel = None
        subset = self.find(name)
        if subset is None:
            subset = self.bpy_setting.add()
            self.ui_list_index = len(self.bpy_setting) - 1
            subset.name = name
        subset.label = self.label
        subset.flag = True
        # species_props = self.batoms.
        sel = Select(self.label, parent=self, name=name,
                     batoms=self.batoms, indices=indices)
        return sel


def elect_expression(expre, batoms):
    """
    select expression
all
none

Logical
not Sele
sele1 and sele2

Properties
element O

Coordinates
z<5

Chemcial classses

solvent: water
hydrogens:
metals

"""
    indices = None
    if isinstance(expre, (list, np.ndarray)):
        return expre
    if isinstance(expre, str):
        if 'chain' in expre:
            chainid = expre.split()[1]
            logger.debug('chain %s' % chainid)
            indices = np.where(batoms.attributes['chainids'] == chainid)[0]
        if 'sheet' in expre:
            sheet = expre.split()[1]
            logger.debug('sheet %s' % sheet)
            indices = batoms.ribbon.protein.sheets[sheet].indices
        if 'helix' in expre:
            helix = expre.split()[1]
            logger.debug('helix %s' % helix)
            indices = batoms.ribbon.protein.helixs[helix].indices
        if 'hetatm' in expre:
            logger.debug('hetatm')
            indices = np.where(batoms.attributes['types'] == 'HETATM')[0]
        if 'species' in expre:
            species = expre.split()[1]
            logger.debug('species %s' % species)
            indices = np.where(batoms.attributes['species'] == species)[0]
    return indices
