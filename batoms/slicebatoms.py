import bpy
from batoms.base.object import childObjectGN
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class SliceBatoms(childObjectGN):
    def __init__(self, label, indices=None, batoms=None) -> None:
        """SliceBatom class

        Two cases:
        1) len(indices) > 1: return Batoms[indices]
        1) one atoms: access properties of one atoms by reading attribute
        from obj directly. This will be faster than Batoms[indices]

        Args:
            label (str):
                Name of the Batoms object
            indices (int, optional):
                indices of the atom. Defaults to 0.
            batoms (Batoms, optional):
                Batoms object. Defaults to None.
        """
        childObjectGN.__init__(self, label, indices, parent=batoms)

    def __repr__(self):
        # bpy.context.view_layer.objects.active = self.obj
        # mode = self.obj.mode
        # bpy.ops.object.mode_set(mode='OBJECT')
        s = "SliceBatoms(species = %s, " % np.unique(self.species)
        # s += "elements = %s, " % str(self.elements)
        s += "positions = %s)" % np.round(self.position, 2)
        # bpy.ops.object.mode_set(mode=mode)
        return s

    @property
    def species(self):
        return self.get_attribute('species')

    # @species.setter
    # def species(self, species):
    #     self.set_species(species)

    # def set_species(self, species):
    #     # self.attributes['species'].data[self.indices].value = species
    #     self.parent.replace([self.indices], species)

    # @property
    # def elements(self):
    #     elements = self.parent.species[self.species].elements
    #     occupancies = {}
    #     for name, eledata in elements.items():
    #         occupancies[name] = round(eledata["occupancy"], 3)
    #     return occupancies

    @property
    def polyhedra(self):
        model_style = self.get_attribute('model_style')
        return model_style  == 2

    @polyhedra.setter
    def polyhedra(self, polyhedra):
        model_style = np.array(polyhedra).astype(int)
        if len(self.indices) == 1:
            self.attributes['model_style'].data[self.indices[0]].value = 2 if polyhedra else 1
        else:
            self.parent.set_attribute_with_indices('model_style',
                self.indices, [2 if polyhedra else 1]*len(self.indices))
        self.parent.draw_polyhedra()

    @property
    def bond(self):
        model_style = self.get_attribute('model_style')
        return model_style > 0

    @bond.setter
    def bond(self, bond):
        model_style = np.array(bond).astype(int)
        self.set_attribute('model_style', model_style)
        self.parent.draw_ball_and_stick()
