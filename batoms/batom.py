import bpy
from batoms.base.object import childObjectGN
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class Batom(childObjectGN):
    def __init__(self, label, index=0, batoms=None) -> None:
        """Batom class

        Access properties of one atoms by reading attribute 
        from obj directly. This will be faster than Batoms[index]

        Args:
            label (str):
                Name of the Batoms object
            index (int, optional):
                Index of the atom. Defaults to 0.
            batoms (Batoms, optional):
                Batoms object. Defaults to None.
        """
        childObjectGN.__init__(self, label, index, parent=batoms)

    def __repr__(self):
        bpy.context.view_layer.objects.active = self.obj
        mode = self.obj.mode
        bpy.ops.object.mode_set(mode='OBJECT')
        s = "Batom(species = '%s', " % self.species
        # s += "elements = %s, " % str(self.elements)
        s += "positions = %s)" % np.round(self.position, 2)
        bpy.ops.object.mode_set(mode=mode)
        return s

    @property
    def species(self):
        return self.attributes['species'].data[self.index].value

    @species.setter
    def species(self, species):
        self.set_species(species)

    def set_species(self, species):
        # self.attributes['species'].data[self.index].value = species
        self.parent.replace([self.index], species)

    @property
    def elements(self):
        elements = self.parent.species[self.species].elements
        occupancies = {}
        for name, eledata in elements.items():
            occupancies[name] = round(eledata["occupancy"], 3)
        return occupancies

    @property
    def polyhedra(self):
        polyhedra = self.attributes['model_style'].data[self.index].value == 2
        return polyhedra

    @polyhedra.setter
    def polyhedra(self, polyhedra):
        if polyhedra:
            self.attributes['model_style'].data[self.index].value = 2
        else:
            self.attributes['model_style'].data[self.index].value = 1
        self.parent.draw_polyhedra()

    @property
    def bond(self):
        bond = self.attributes['model_style'].data[self.index].value > 0
        return bond

    @bond.setter
    def bond(self, bond):
        if bond:
            self.attributes['model_style'].data[self.index].value = 1
        else:
            self.attributes['model_style'].data[self.index].value = 0
        self.parent.draw_ball_and_stick()
