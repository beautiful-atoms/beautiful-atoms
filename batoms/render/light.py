"""
Light setting

"""
import bpy
from mathutils import Vector
import numpy as np
from batoms.base.collection import BaseCollection
from batoms.base.object import BaseObject
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class Light(BaseObject):
    def __init__(self, label, name, type='SUN', energy=5,
                 location=(0, 0, 30),
                 direction=[0, 0, 1],
                 lock_to_camera=False) -> None:
        """Light Class

        A Light object is linked to a Light object in Blender.

        Parameters:

        label: str
            Name of the Batoms.
        name: str
            name of the light.
        type: str
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
        self.label = label
        self.name = name
        obj_name = "%s_light_%s" % (label, name)
        btype = "Blight"
        self.coll_name = "%s_render" % (label)
        BaseObject.__init__(self, obj_name=obj_name, btype=btype)
        self.create_light(type, location, energy)
        self.direction = direction
        self.lock_to_camera = lock_to_camera

    @property
    def coll(self):
        return self.get_coll()

    def get_coll(self):
        return bpy.data.collections.get(self.coll_name)

    def get_bobj(self):
        bobj = getattr(self.obj, self.btype)
        return bobj

    @property
    def energy(self):
        return self.get_energy()

    @energy.setter
    def energy(self, energy):
        self.set_energy(energy)

    def get_energy(self):
        return self.obj.data.energy

    def set_energy(self, energy):
        self.obj.data.energy = energy

    @property
    def direction(self):
        return self.get_direction()

    @direction.setter
    def direction(self, direction):
        self.set_direction(direction)

    def get_direction(self):
        return self.obj.Blight.direction[:]

    def set_direction(self, direction):
        from batoms.utils import rotate_frame
        self.lock_to_camera = False
        self.obj.Blight.direction = direction
        #
        new_frame = rotate_frame(self.coll.Brender.viewport)
        dirction = np.dot(direction, new_frame)
        dirction = dirction/np.linalg.norm(dirction)
        location = self.look_at + dirction*self.coll.Brender.distance
        self.location = location
        self.look_at = self.coll.Brender.center

    @property
    def lock_to_camera(self):
        return self.get_lock_to_camera()

    @lock_to_camera.setter
    def lock_to_camera(self, lock_to_camera):
        self.set_lock_to_camera(lock_to_camera)

    def get_lock_to_camera(self):
        return self.obj.Blight.lock_to_camera

    def set_lock_to_camera(self, lock_to_camera):
        self.obj.Blight.lock_to_camera = lock_to_camera

    def create_light(self, type, location, energy):
        '''
        type: str
        energy: float
        '''
        # check light exist or not
        if self.obj_name in bpy.data.objects:
            if not bpy.data.objects[self.obj_name].Blight.flag:
                raise Exception("%s is not a Light for Batoms!" %
                                self.obj_name)
            light = bpy.data.objects[self.obj_name]
        else:
            if self.obj_name in bpy.data.lights:
                light_data = bpy.data.lights[self.obj_name]
            else:
                light_data = bpy.data.lights.new(self.obj_name, type=type)
            light = bpy.data.objects.new(self.obj_name, light_data)
            light.data.energy = energy
            light.data.use_nodes = True
            light.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
            light.location = Vector(location)
            light.Blight.flag = True
            light.Blight.label = self.label
            light.Blight.name = self.name

    def __repr__(self) -> str:
        s = "Light('%s', energy = %s, direction = %s, lock_to_camera = %s)" \
            % (self.type, self.energy, list(self.direction),
               self.lock_to_camera)
        return s


class Lights(BaseCollection):
    """
    A collection of Light
    """

    def __init__(self, label) -> None:
        self.label = label
        self.name = "%s_light" % label
        BaseCollection.__init__(self, coll_name=self.name)

    @property
    def lights(self):
        return self.get_lights()

    def get_lights(self):
        lights = {}
        for light in self.coll.objects:
            lights[light.Blight.name] = Light(light.Blight.label,
                                                    light.Blight.name,
                                                    direction=light.Blight.direction,
                                                    lock_to_camera=light.Blight.lock_to_camera)
        return lights

    def __getitem__(self, name):
        """Return a subset of the Batom.

        species -- str, describing which batom to return.
        """
        if isinstance(name, str):
            if name not in self.lights:
                raise SystemExit('%s is not in this light' % name)
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
        if light.obj.name not in self.coll.objects:
            self.coll.objects.link(light.obj)

    def remove(self, name):
        light = self.lights[name]
        bpy.data.objects.remove(light.obj, do_unlink=True)
