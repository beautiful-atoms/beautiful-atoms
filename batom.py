
import numpy as np

class Batom():
    """Batom Class
    """
    def __init__(self, label, index = 0, batoms = None) -> None:
        self.label = label
        self.batoms = batoms
        self.index = index

    def __repr__(self):
        s = "Batom(species = '%s', elements = %s, \
                    positions = %s" % (self.species,  \
                    str(self.elements), self.positions)
        return s
    
    @property
    def vertice(self):
        return self.batoms.obj.data.vertices[self.index]
    
    @property
    def attribute(self):
        return self.batoms.obj.data.attributes[self.index]
    
    @property
    def species(self):
        return self.batoms.obj.data.attributes['species'].data[0].value
    
    @species.setter
    def species(self, species):
        self.set_species(species)
    
    def set_species(self, species):
        self.batoms.obj.data.attributes['species'].data[0].value = species
    
    @property
    def elements(self):
        return self.batoms.species[self.species]['elements']
    
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
                np.array(self.batoms.obj.matrix_world))
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
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')

    @property
    def scale(self):
        return self.batoms.obj.data.attributes['scale'].data[self.index].value
    
    @scale.setter
    def scale(self, scale):
        self.batoms.obj.data.attributes['scale'].data[self.index].value = scale
    
    @property
    def show(self):
        return self.batoms.obj.data.attributes['show'].data[self.index].value
    
    @show.setter
    def show(self, show):
        self.batoms.obj.data.attributes['show'].data[self.index].value = show