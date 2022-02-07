import bpy
import numpy as np
from batoms.butils import update_object, object_mode

class Batom():
    """Batom Class
    """
    def __init__(self, label, index = 0, batoms = None) -> None:
        self.label = label
        self.batoms = batoms
        self.index = index

    def __repr__(self):
        bpy.context.view_layer.objects.active = self.batoms.obj
        mode = self.batoms.obj.mode
        bpy.ops.object.mode_set(mode='OBJECT')
        s = "Batom(species = '%s', elements = %s, positions = %s" % (self.species,  \
                    str(self.elements), self.position)
        bpy.ops.object.mode_set(mode=mode)
        return s
    
    @property
    def vertice(self):
        # vertice = self.batoms.obj.data.vertices[self.index]
        vertice = self.batoms.obj.data.shape_keys.key_blocks[0].data[self.index]
        return vertice
    
    @property
    def attribute(self):
        return self.batoms.obj.data.attributes[self.index]
    
    @property
    def species(self):
        return self.batoms.obj.data.attributes['species'].data[self.index].value
    
    @species.setter
    def species(self, species):
        self.set_species(species)
    
    def set_species(self, species):
        # self.batoms.obj.data.attributes['species'].data[self.index].value = species
        self.batoms.replace([self.index], species)
    
    @property
    def elements(self):
        return self.batoms.species[self.species]['elements']
    
    @property
    def local_position(self):
        return self.vertice.co

    @property
    def position(self):
        return self.get_position()
    
    @position.setter
    def position(self, position):
        self.set_position(position)
    
    def get_position(self):
        """
        Get global position.
        """
        from batoms.tools import local2global
        position = local2global(np.array([self.local_position]), 
                np.array(self.batoms.obj.matrix_world))
        return position[0]
    
    def set_position(self, position):
        """
        Set global position to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        position = np.array([position])
        position = local2global(position, 
                np.array(self.batoms.obj.matrix_world), reversed = True)
        # rashpe to (natoms*3, 1) and use forseach_set
        self.vertice.co = position[0]
        self.batoms.obj.data.update()
        update_object(self.batoms.obj)

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
        update_object(self.batoms.obj)

