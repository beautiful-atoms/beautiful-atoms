"""
"""
import bpy
from batoms.base.collection import Setting, tuple2string
import numpy as np
from batoms.data import default_polyhedras
# from time import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class PolyhedraSettings(Setting):
    """
    PolyhedraSetting object

    The PolyhedraSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """

    def __init__(self, label, batoms=None,
                 polyhedras=None,
                 polyhedrasetting=None,
                 ) -> None:
        Setting.__init__(self, label, coll_name='%s' % label)
        self.name = 'Bpolyhedra'
        self.batoms = batoms
        self.polyhedras = polyhedras
        if len(self) == 0:
            self.add_species_list(self.batoms.species.keys())
        if polyhedrasetting is not None:
            for key, data in polyhedrasetting.items():
                self[key] = data

    def build_materials(self, sp, node_inputs=None):
        """
        """
        from batoms.material import create_material
        name = '%s_%s' % (self.label, sp['species'])
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                        color=sp['color'],
                        node_inputs=node_inputs,
                        material_style=sp["material_style"],
                        backface_culling=False)
        # logger.debug('build polyhedra materials: {}'.format(sp['color']))
        return mat

    def __setitem__(self, index, setdict):
        """
        Set properties
        # TODO can not set color by setting['C'].color = [1, 0, 0, 1]
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.bpy_setting.add()
        subset.species = index
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_materials(self[name].as_dict())

    def set_default(self, species_props):
        """
        """
        for sp, props in species_props.items():
            if sp not in default_polyhedras:
                continue
            props['color'] = (props['color'][0], props['color'][1], props['color'][2], 0.8)
            self[sp] = props


    def add(self, key, value = None):
        props = self.batoms.species.species_props[key]
        if value:
            props.update(value)
        props['color'] = (props['color'][0], props['color'][1], props['color'][2], 0.8)
        props['species'] = key
        self[key] = props

    def add_species_list(self, species_list, only_default=True):
        for sp in species_list:
            self.add_species(sp, only_default)

    def add_species(self, name, only_default=True):
        species_props = self.batoms.species.species_props[name]
        if only_default and species_props['element'] not in default_polyhedras:
            return
        self.add(name, species_props)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center                color         width \n"
        for p in self.bpy_setting:
            s += "{:10s}   ".format(p.species)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}]   ".format(p.color[0],
                                                                  p.color[1],
                                                                  p.color[2],
                                                                  p.color[3])
            s += "{:1.3f} \n".format(p.width)
        s += "-"*60 + "\n"
        return s

    @property
    def materials(self):
        return self.get_materials()

    def get_materials(self):
        materials = {}
        index = 0
        for sp in self:
            name = '%s_%s' % (self.label, sp.species)
            mat = bpy.data.materials.get(name)
            materials[sp.species] = [mat, index]
            index += 1
        return materials
