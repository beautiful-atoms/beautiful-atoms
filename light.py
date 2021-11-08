"""
Light setting

"""
import bpy
from mathutils import Vector
import os
import numpy as np


class Light():
    """Light Class
    
    A Light object is linked to a Light object in Blender. 

    Parameters:

    label: str
        Name of the Batoms.
    name: str
        name of the light.
    light_type: str
        * POINT, Omnidirectional point light source.
        * SUN, Constant direction parallel ray light source.
        * SPOT, Directional cone light source.
        * AREA, Directional area light source.
    energy: float
    location: array
        location
    direction: array
        The light direction related to camera.
    lock_to_camera: bool
        lock_to_camera
    """
    def __init__(self, label, name, light_type = 'SUN', energy = 5, 
                location = (0, 0, 100),
                direction = [0, 0, 1],
                lock_to_camera = False) -> None:
        self.label = label
        self.name = name
        self.obj_name = "%s_light_%s"%(label, name)
        self.create_light(light_type, location, energy)
        self.direction = direction
        self.lock_to_camera = lock_to_camera
    @property
    def light(self):
        return self.get_light()
    def get_light(self):
        return bpy.data.objects[self.obj_name]
    @property
    def type(self):
        return self.get_type()
    @type.setter
    def type(self, type):
        self.set_type(type)
    def get_type(self):
        return self.light.data.type
    def set_type(self, type):
        self.light.data.type = type.upper()
    @property
    def energy(self):
        return self.get_energy()
    @energy.setter
    def energy(self, energy):
        self.set_energy(energy)
    def get_energy(self):
        return self.light.data.energy
    def set_energy(self, energy):
        self.light.data.energy = energy
    @property
    def location(self):
        return self.get_location()
    @location.setter
    def location(self, location):
        self.set_location(location)
    def get_location(self):
        return self.light.location
    def set_location(self, location):
        self.light.location = location
    @property
    def direction(self):
        return self.get_direction()
    @direction.setter
    def direction(self, direction):
        self.set_direction(direction)
    def get_direction(self):
        return self.light.blight.direction[:]
    def set_direction(self, direction):
        self.lock_to_camera = False
        self.light.blight.direction = direction
    @property
    def lock_to_camera(self):
        return self.get_lock_to_camera()
    @lock_to_camera.setter
    def lock_to_camera(self, lock_to_camera):
        self.set_lock_to_camera(lock_to_camera)
    def get_lock_to_camera(self):
        return self.light.blight.lock_to_camera
    def set_lock_to_camera(self, lock_to_camera):
        self.light.blight.lock_to_camera = lock_to_camera
    def create_light(self, light_type, location, energy):
        '''
        light_type: str
        energy: float
        '''
        # check light exist or not
        if self.obj_name in bpy.data.objects:
            if not bpy.data.objects[self.obj_name].blight.flag:
                raise Exception("%s is not a Light for Batoms!"%self.obj_name)
            light = bpy.data.objects[self.obj_name]
        else:
            if self.obj_name in bpy.data.lights:
                light_data = bpy.data.lights[self.obj_name]
            else:
                light_data = bpy.data.lights.new(self.obj_name, type = light_type)
            light = bpy.data.objects.new(self.obj_name, light_data)
            light.data.energy = energy
            light.data.use_nodes = True
            light.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
            light.location = Vector(location)
            light.blight.flag = True
            light.blight.label = self.label
            light.blight.name = self.name

    def lock_to(self, obj = None):
        """
        track to obj
        """
        if obj is not None:
            self.light.constraints.new(type = 'COPY_LOCATION')
            self.light.constraints["Copy Location"].target = obj
            self.light.constraints.new(type = 'COPY_ROTATION')
            self.light.constraints["Copy Rotation"].target = obj
        else:
            for c in self.light.constraints:
                self.light.constraints.remove(c)
    def __repr__(self) -> str:
        s = "Light('%s', energy = %s, direction = %s, lock_to_camera = %s)" \
                    % (self.type, self.energy, list(self.direction), self.lock_to_camera)
        return s

class Lights():
    """
    A collection of Light
    """
    def __init__(self, label) -> None:
        self.label = label
        self.name = "%s_light"%label
        self.set_collection()
    def set_collection(self):
        """
        build main collection and its child collections.
        """
        if self.name not in bpy.data.collections:
            bpy.data.collections.new(self.name)
    @property
    def coll(self):
        return self.get_coll()
    def get_coll(self):
        return bpy.data.collections[self.name]
    @property
    def lights(self):
        return self.get_lights()
    def get_lights(self):
        lights = {}
        for l in self.coll.objects:
            lights[l.blight.name] = Light(l.blight.label, 
                                        l.blight.name,
                                        direction=l.blight.direction,
                                        lock_to_camera=l.blight.lock_to_camera)
        return lights
    def __getitem__(self, name):
        """Return a subset of the Batom.

        species -- str, describing which batom to return.
        """
        if isinstance(name, str):
            if name not in self.lights:
                raise SystemExit('%s is not in this light'%name)
            return self.lights[name]
        elif isinstance(name, list):
            raise SystemExit('dict not supported yet!')
    def __setitem__(self, name, setdict):
        """
        Set properties
        """
        light = self.lights(name)
        for key, value in setdict.items():
            setattr(light, key, value)
    def __repr__(self) -> str:
        s = '       Name        Type    Direction      Lock_to_camera\n'
        i = 0
        for name, l in self.lights.items():
            s += '{:3d}   {:10s}  {:5s}   {}    {}\n'.format(i, 
                    name, l.type, l.direction, l.lock_to_camera)
            i += 1
        return s
    def __len__(self):
        return len(self.lights)
    def add(self, name, **parameters):
        light = Light(self.label, name, **parameters)
        if light.light.name not in self.coll.objects:
            self.coll.objects.link(light.light)
    def remove(self, name):
        light = self.lights[name]
        bpy.data.objects.remove(light.light, do_unlink = True)

