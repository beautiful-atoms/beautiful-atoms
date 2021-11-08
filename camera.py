"""
Camera setting

"""
import bpy
from mathutils import Vector
import os
import numpy as np

def lock_camera_to_view(switch):
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.lock_camera = switch

class Camera():
    def __init__(self, label, name = 'Default', 
                        coll = None,
                        camera_type = 'ORTHO', 
                        lens = 50, 
                        location = (0, 0, 100), 
                        ) -> None:
        self.label = label
        self.name = "%s_camera_%s"%(label, name)
        self.create_camera(camera_type, location, lens, coll = coll)
    @property
    def camera(self):
        return self.get_camera()
    def get_camera(self):
        return bpy.data.objects[self.name]
    @property
    def type(self):
        return self.get_type()
    @type.setter
    def type(self, type):
        self.set_type(type)
    def get_type(self):
        return self.camera.data.type
    def set_type(self, type):
        self.camera.data.type = type.upper()
    @property
    def lens(self):
        return self.get_lens()
    @lens.setter
    def lens(self, lens):
        self.set_lens(lens)
    def get_lens(self):
        return self.camera.data.lens
    def set_lens(self, lens):
        self.camera.data.lens = lens
    @property
    def ortho_scale(self):
        return self.get_ortho_scale()
    @ortho_scale.setter
    def ortho_scale(self, ortho_scale):
        self.set_ortho_scale(ortho_scale)
    def get_ortho_scale(self):
        return self.camera.data.ortho_scale
    def set_ortho_scale(self, ortho_scale):
        self.camera.data.ortho_scale = ortho_scale
    @property
    def location(self):
        return self.get_location()
    @location.setter
    def location(self, location):
        self.set_location(location)
    def get_location(self):
        return self.camera.location
    def set_location(self, location):
        self.camera.location = location
    def create_camera(self, camera_type = None, 
                    location = None, lens = 50, coll = None):
        '''
        camera_type: str
            * PERSP
            * ORTHO
        lens: float
        '''
        if not camera_type:
            camera_type = self.camera_type
        if not lens:
            lens = self.lens
        if self.name in bpy.data.objects:
            camera = bpy.data.objects[self.name]
        else:
            if self.name in bpy.data.cameras:
                camera_data = bpy.data.cameras[self.name]
            else:
                camera_data = bpy.data.cameras.new(self.name)
            camera = bpy.data.objects.new(self.name, camera_data)
            camera.data.lens = lens
            camera.data.type = camera_type
            camera.location = Vector(location)
            if coll is not None:
                coll.objects.link(camera)
    def lock_camera_to_camera(self, camera = None):
        """
        track to camera
        """
        if camera is not None:
            self.camera.constraints.new(type = 'COPY_LOCATION')
            self.camera.constraints["Copy Location"].target = camera
            self.camera.constraints.new(type = 'COPY_ROTATION')
            self.camera.constraints["Copy Rotation"].target = camera
        else:
            for c in self.camera.constraints:
                self.camera.constraints.remove(c)
