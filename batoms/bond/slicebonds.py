import bpy
import bmesh
import numpy as np
from batoms.utils import number2String
from batoms.base.object import childObjectGN
import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

class SliceBonds(childObjectGN):
    """Batom Class
    """

    def __init__(self, label, indices=0, bonds=None) -> None:
        childObjectGN.__init__(self, label, indices,
                    obj_name="{}_bond".format(label),
                    parent=bonds)

    def __repr__(self):
        # bpy.context.view_layer.objects.active = self.parent.obj
        # mode = self.parent.obj.mode
        # bpy.ops.object.mode_set(mode='OBJECT')
        s = "SliceBonds(species = %s, " % self.species["name"]
        s += "order = %s, " % self.order
        s += "style = %s)" % self.style
        # bpy.ops.object.mode_set(mode=mode)
        return s

    @property
    def species(self):
        #TODO: support slice
        sp1 = self.get_attribute('species_index1')
        sp2 = self.get_attribute('species_index2')
        if len(self.indices) == 1:
            name = np.array(['%s-%s' % (number2String(sp1), number2String(sp2))])
        else:
            n = len(sp1)
            name = np.zeros(n, dtype="U20")
            for i in range(n):
                name[i] = '%s-%s' % (number2String(sp1[i]), number2String(sp2[i]))
        order = self.get_attribute('order')
        style = self.get_attribute('style')
        # name = '%s-%s' % (sp1, sp2)
        sp = self.parent.settings[name[0]].as_dict()
        sp.update({'order': order[0],
                   'style': style[0]})
        return sp

    @species.setter
    def species(self, species):
        self.set_species(species)

    def set_species(self, species):
        # self.attributes['species'].data[self.indices].value = species
        self.parent.replace([self.indices], species)

    @property
    def order(self):
        order = self.get_attribute('order')
        return order

    @order.setter
    def order(self, value):
        self.set_attribute('order', value)
        # find plane for high order
        a3, a4 = self.secondBond()
        self.set_attribute('atoms_index3', a3)
        self.set_attribute('atoms_index4', a4)
        # if order not exist, add one
        sp = self.species
        sp['order'] = value
        self.parent.settings.build_instancer(sp)
        self.parent.add_geometry_node(sp)

    @property
    def style(self):
        style = self.get_attribute('style')
        return style

    @style.setter
    def style(self, value):
        self.set_attribute('style', value)
        # if stule not exist, add one
        sp = self.species
        sp['style'] = value
        self.parent.settings.build_instancer(sp)
        self.parent.add_geometry_node(sp)

    def secondBond(self):
        """
        determine the plane of high order bond
        """
        arrays = self.parent.arrays
        ai = arrays['atoms_index1'][self.indices[0]]
        aj = arrays['atoms_index2'][self.indices[0]]
        indi = np.where((arrays['atoms_index1'] == ai) |
                        (arrays['atoms_index2'] == ai) |
                        (arrays['atoms_index1'] == aj) |
                        (arrays['atoms_index2'] == aj))[0]
        if len(indi) == 1:
            second_bond = [0, 1]
        else:
            for i in indi:
                if arrays['atoms_index1'][i] == ai and \
                        arrays['atoms_index2'][i] == aj:
                    continue
                second_bond = [arrays['atoms_index1']
                               [i], arrays['atoms_index2'][i]]
                break
        return second_bond

    def high_order_bond_plane(self):
        """
        determine the plane of high order bond
        """
        arrays = self.parent.arrays
        ai = arrays['atoms_index1'][self.indices]
        aj = arrays['atoms_index2'][self.indices]
        indi = np.where(arrays['atoms_index1'] == ai)[0]
        if len(indi) == 1:
            second_bond = np.array([0.0, 0.0, 1])
        else:
            for i in indi:
                if arrays['atoms_index2'][i] == aj:
                    continue
                second_bond = arrays['vectors'][i]
                break
        vec = np.cross(second_bond,
                       arrays['vectors'][self.indices]) + np.array([1e-8, 0, 0])
        offset = np.cross(vec, arrays['vectors']
                          [self.indices]) + np.array([1e-8, 0, 0])
        return offset
