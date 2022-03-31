"""
"""
import bpy
from batoms.base.collection import Setting, tuple2string
import numpy as np
from batoms.data import default_polyhedras
# from time import time


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
        self.name = 'bpolyhedra'
        self.batoms = batoms
        self.polyhedras = polyhedras
        if len(self) == 0:
            self.set_default(self.batoms.species.species_props)
        if polyhedrasetting is not None:
            for key, data in polyhedrasetting.items():
                self[key] = data

    def build_materials(self, sp, node_inputs=None,
                        material_style='default'):
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
                        material_style=material_style,
                        backface_culling=False)
        return mat

    def __setitem__(self, index, setdict):
        """
        Set properties
        # TODO can not set color by setting['C'].color = [1, 0, 0, 1]
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
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
        for sp, data in species_props.items():
            if sp not in default_polyhedras:
                continue
            self[sp] = {
                'flag': True,
                'label': self.label,
                'species': sp,
                'color': np.append(data['color'][:3], 0.8),
                'width': 0.005,
            }

    def add(self, polyhedras):
        if isinstance(polyhedras, str):
            polyhedras = [polyhedras]
        species = {sp: self.batoms.species.species_props[sp] for sp in polyhedras}
        self.set_default(species)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center                color         width \n"
        for p in self.collection:
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
