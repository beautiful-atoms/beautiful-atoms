"""

"""

from math import exp
import bpy
from batoms.base import Setting, BaseCollection
from batoms.bspecies import Bspecies
from batoms.bondsetting import BondSetting, build_bondlists, calc_bond_data
from batoms.polyhedrasetting import PolyhedraSetting, build_polyhedralists
from batoms.isosurfacesetting import IsosurfaceSetting
from batoms.planesetting import PlaneSetting
from batoms.mssetting import MSsetting
from batoms.ribbon import Ribbon
import numpy as np
from batoms.butils import object_mode, show_index
from time import time

from batoms.tools import string2Number

subcollections = []

class Select():
    """

    sel = batoms.selects[key]
    sel.model_style = 1

    """
    def __init__(self, label, name = 'sel0', 
                batoms = None,
                ) -> None:
        self.batoms = batoms
        self.label = label
        self.name = name
        self.coll_name = '%s_%s'%(label, name)
        # BaseCollection.__init__(self, coll_name = '%s_%s'%(label, name))
        
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
        return self.batoms.selects.collection[self.name].radius_style
    
    def set_radius_style(self, radius_style):
        self.batoms.selects.collection[self.name].radius_style = radius_style
        self.draw(draw_isosurface = False)
    
    @property
    def show(self):
        return self.get_show()
    
    @show.setter
    def show(self, state):
        self.set_show(state)
    
    def get_show(self):
        return self.batoms.attributes['show'][self.indices]
    
    def set_show(self, show, only_atoms = True):
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
    def indices(self):
        return self.get_indices()
    
    @indices.setter
    def indices(self, indices):
        self.set_indices(indices)
    
    def get_indices(self):
        indices = self.batoms.arrays['select']
        indices = np.where(indices == string2Number(self.name))[0]
        return indices
    
    def set_indices(self, indices):
        select = self.batoms.attributes['select']
        select[self.indices] = string2Number('sel0')
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
        for sp in self.batoms.species:
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
    def radius(self):
        """
        radius depends on radius style, thus could be with select
        """
        radius = {}
        for sp in self.collection:
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
        s += '{0:10s} {1:10s}   {2:10s}   {3:10s} \n'.format(\
                self.name, str(self.model_style[0]), self.radius_style[0], str(self.show[0]))
        s += '-'*60 + '\n'
        return s

class Selects(Setting):
    """
    sel = batoms.select.add(name, index)
    sel = batoms.select.from_dict({name: index})
    sel.model_style = 1
    """
    def __init__(self, label, batoms = None) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bselect'
        self.batoms = batoms

    @property
    def selects(self):
        return self.get_selects()
    
    def get_selects(self):
        selects = {}
        for b in self.collection:
            selects[b.name] = Select(self.label, name = b.name, 
                    batoms = self.batoms)
        return selects

    def __getitem__(self, key):
        """Return a subset of the Batom.

        key -- str, describing which batom to return.
        """
        if isinstance(key, str):
            if key not in self.selects:
                raise KeyError('%s is not in this structure'%key)
            return self.selects[key]
        elif isinstance(key, list):
            raise KeyError('dict not supported yet!')
    
    def __len__(self):
        n = 0
        for ba in self.selects.values():
            n += len(ba)
        return n
    
    def __repr__(self) -> str:
        s = '     Name        model_style    show  \n'
        i = 0
        for name, sel in self.selects.items():
            s += '{:3d}   {:10s}  {:5s}   {} \n'.format(i, 
                    name, str(sel.model_style[0]), str(sel.show[0]))
            i += 1
        return s
    
    def add(self, name, expre):
        indices = elect_expression(expre, self.batoms)
        # print(name, string2Number(name))
        select = self.batoms.attributes['select']
        select[indices] = string2Number(name)
        select = {'select': select}
        self.batoms.set_attributes(select)
        sel = None
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        subset.label = self.label
        subset.flag = True
        # species_props = self.batoms.
        self.batoms.species.build_instancers()
        self.batoms.build_geometry_node()
        sel = Select(self.label, name = name, batoms = self.batoms)
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
    if isinstance(expre, (list, np.ndarray)):
        return expre
    if isinstance(expre, str):
        if 'chain' in expre:
            chainid = expre.split()[1]
            print('chain %s'%chainid)
            indices = np.where(batoms.attributes['chainids'] == chainid)[0]
        if 'species' in expre:
            species = expre.split()[1]
            print('species %s'%species)
            indices = np.where(batoms.attributes['species'] == species)[0]
        
    return indices
        