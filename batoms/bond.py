import bpy
import numpy as np
from batoms.butils import update_object
from batoms.tools import number2String

class Bond():
    """Batom Class
    """
    def __init__(self, label, index = 0, bonds = None) -> None:
        self.label = label
        self.bonds = bonds
        self.index = index

    def __repr__(self):
        bpy.context.view_layer.objects.active = self.bonds.obj
        mode = self.bonds.obj.mode
        bpy.ops.object.mode_set(mode='OBJECT')
        s = "Bond(species = '%s', order = %s, style = %s" % (self.species['name'],  \
                    str(self.species['order']), str(self.species['style']))
        bpy.ops.object.mode_set(mode=mode)
        return s
    
    @property
    def vertice(self):
        return self.bonds.obj.data.vertices[self.index]
    
    @property
    def attribute(self):
        return self.bonds.obj.data.attributes[self.index]
    
    @property
    def species(self):
        sp1 = self.bonds.obj.data.attributes['species_index1'].data[self.index].value
        sp2 = self.bonds.obj.data.attributes['species_index2'].data[self.index].value
        sp1 = number2String(sp1)
        sp2 = number2String(sp2)
        order = self.bonds.obj.data.attributes['order'].data[self.index].value
        style = self.bonds.obj.data.attributes['style'].data[self.index].value
        name = '%s-%s'%(sp1, sp2)
        sp = self.bonds.setting[name].as_dict()
        sp.update({'order': order,
                   'style': style})
        return sp
    
    @species.setter
    def species(self, species):
        self.set_species(species)
    
    def set_species(self, species):
        # self.bonds.obj.data.attributes['species'].data[self.index].value = species
        self.bonds.replace([self.index], species)
    
    @property
    def local_positions(self):
        return self.vertice.co

    @property
    def positions(self):
        return self.get_positions()
    
    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)
    
    def get_positions(self):
        """
        Get global positions.
        """
        from batoms.tools import local2global
        positions = local2global(np.array([self.local_positions]), 
                np.array(self.bonds.obj.matrix_world))
        return positions[0]
    
    def set_positions(self, positions):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        natom = len(self)
        if len(positions) != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natom))
        positions = local2global(positions, 
                np.array(self.obj.matrix_world), reversed = True)
        # rashpe to (natoms*3, 1) and use forseach_set
        positions = positions.reshape((natom*3, 1))
        # I don't know why 'Basis' shape keys is not updated when editing mesh,
        # so we edit the 'Basis' shape keys directly.
        # self.obj.data.vertices.foreach_set('co', positions)
        self.obj.data.shape_keys.key_blocks[0].data.foreach_set('co', positions)
        self.obj.data.update()
        # bpy.context.view_layer.update()
        # I don't why this is need to update the mesh positions
        update_object(self.obj)

    @property
    def order(self):
        return self.bonds.obj.data.attributes['order'].data[self.index].value
    
    @order.setter
    def order(self, order):
        self.bonds.obj.data.attributes['order'].data[self.index].value = order
        # find plane for high order
        a3, a4 = self.secondBond()
        self.bonds.obj.data.attributes['atoms_index3'].data[self.index].value = a3
        self.bonds.obj.data.attributes['atoms_index4'].data[self.index].value = a4
        # if order not exist, add one
        sp = self.species
        sp['order'] = order
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)
    
    @property
    def style(self):
        return self.bonds.obj.data.attributes['style'].data[self.index].value
    
    @style.setter
    def style(self, style):
        self.bonds.obj.data.attributes['style'].data[self.index].value = style
        # if stule not exist, add one
        sp = self.species
        sp['style'] = style
        self.bonds.setting.build_instancer(sp)
        self.bonds.add_geometry_node(sp)

    @property
    def scale(self):
        return self.bonds.obj.data.attributes['scale'].data[self.index].value
    
    @scale.setter
    def scale(self, scale):
        self.bonds.obj.data.attributes['scale'].data[self.index].value = scale
    
    @property
    def show(self):
        obj = self.bonds.obj
        mode = obj.mode
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode = 'OBJECT')
        show = obj.data.attributes['show'].data[self.index].value
        bpy.ops.object.mode_set(mode = mode)
        return show
    
    @show.setter
    def show(self, show):
        obj = self.bonds.obj
        mode = obj.mode
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode = 'OBJECT')
        self.bonds.obj.data.attributes['show'].data[self.index].value = show
        bpy.ops.object.mode_set(mode = mode)
        update_object(self.bonds.obj)
    

    def secondBond(self):
        """
        determine the plane of high order bond
        """ 
        arrays = self.bonds.arrays
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
                # print(arrays['atoms_index1'][i] == ai)
                # print(arrays['atoms_index2'][i] == aj)
                if arrays['atoms_index1'][i] == ai and arrays['atoms_index2'][i] == aj:
                    continue
                second_bond = [arrays['atoms_index1'][i], arrays['atoms_index2'][i]]
                break
        return second_bond

    def high_order_bond_plane(self):
        """
        determine the plane of high order bond
        """ 
        arrays = self.bonds.arrays
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
        vec = np.cross(second_bond, arrays['vectors'][self.index]) + np.array([1e-8, 0, 0])
        offset = np.cross(vec, arrays['vectors'][self.index]) + np.array([1e-8, 0, 0])
        return offset
        