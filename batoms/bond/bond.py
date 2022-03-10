import bpy
import numpy as np
from batoms.utils.butils import update_object
from batoms.utils import number2String
from batoms.base.object import childObjectGN


class Bond(childObjectGN):
    """Batom Class
    """

    def __init__(self, label, index=0, bonds=None) -> None:
        childObjectGN.__init__(self, label, index, parent=bonds)

    def __repr__(self):
        bpy.context.view_layer.objects.active = self.parent.obj
        mode = self.parent.obj.mode
        bpy.ops.object.mode_set(mode='OBJECT')
        s = "Bond(species = '%s', " % self.species['name']
        s += "order = %s, style = %s)" % (str(self.species['order']),
                                          str(self.species['style']))
        bpy.ops.object.mode_set(mode=mode)
        return s

    @property
    def species(self):
        sp1 = self.attributes['species_index1'].data[self.index].value
        sp2 = self.attributes['species_index2'].data[self.index].value
        sp1 = number2String(sp1)
        sp2 = number2String(sp2)
        order = self.attributes['order'].data[self.index].value
        style = self.attributes['style'].data[self.index].value
        name = '%s-%s' % (sp1, sp2)
        sp = self.parent.setting[name].as_dict()
        sp.update({'order': order,
                   'style': style})
        return sp

    @species.setter
    def species(self, species):
        self.set_species(species)

    def set_species(self, species):
        # self.attributes['species'].data[self.index].value = species
        self.parent.replace([self.index], species)

    @property
    def order(self):
        return self.attributes['order'].data[self.index].value

    @order.setter
    def order(self, order):
        self.attributes['order'].data[self.index].value = order
        # find plane for high order
        a3, a4 = self.secondBond()
        self.attributes['atoms_index3'].data[self.index].value = a3
        self.attributes['atoms_index4'].data[self.index].value = a4
        # if order not exist, add one
        sp = self.species
        sp['order'] = order
        self.parent.setting.build_instancer(sp)
        self.parent.add_geometry_node(sp)

    @property
    def style(self):
        return self.attributes['style'].data[self.index].value

    @style.setter
    def style(self, style):
        self.attributes['style'].data[self.index].value = int(style)
        # if stule not exist, add one
        sp = self.species
        sp['style'] = style
        self.parent.setting.build_instancer(sp)
        self.parent.add_geometry_node(sp)

    def secondBond(self):
        """
        determine the plane of high order bond
        """
        arrays = self.parent.arrays
        ai = arrays['atoms_index1'][self.index]
        aj = arrays['atoms_index2'][self.index]
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
        ai = arrays['atoms_index1'][self.index]
        aj = arrays['atoms_index2'][self.index]
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
                       arrays['vectors'][self.index]) + np.array([1e-8, 0, 0])
        offset = np.cross(vec, arrays['vectors']
                          [self.index]) + np.array([1e-8, 0, 0])
        return offset
