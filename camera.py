"""
Camera setting

"""
import bpy
from mathutils import Vector
import os
import numpy as np
from batoms.base import BaseObject

def set_lock_camera_to_view(switch):
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.lock_camera = switch

class Camera(BaseObject):
    def __init__(self, label, name = 'Default', 
                        coll = None,
                        camera_type = 'ORTHO', 
                        lens = 50, 
                        location = (0, 0, 100), 
                        camera_target = [0, 0, 0],
                        from_coll = False,
                        lock_camera_to_view = True,
                        ) -> None:
        self.label = label
        self.name = name
        obj_name = "%s_camera_%s"%(label, name)
        BaseObject.__init__(self, obj_name = obj_name)
        if not from_coll:
            self.create_camera(camera_type, location, lens, coll = coll)
        set_lock_camera_to_view(lock_camera_to_view)
    @property
    def lens(self):
        return self.get_lens()
    @lens.setter
    def lens(self, lens):
        self.set_lens(lens)
    def get_lens(self):
        return self.obj.data.lens
    def set_lens(self, lens):
        self.obj.data.lens = lens
    @property
    def ortho_scale(self):
        return self.get_ortho_scale()
    @ortho_scale.setter
    def ortho_scale(self, ortho_scale):
        self.set_ortho_scale(ortho_scale)
    def get_ortho_scale(self):
        return self.obj.data.ortho_scale
    def set_ortho_scale(self, ortho_scale):
        self.obj.data.ortho_scale = ortho_scale
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
        if self.obj_name in bpy.data.objects:
            camera = bpy.data.objects[self.obj_name]
        else:
            if self.obj_name in bpy.data.cameras:
                camera_data = bpy.data.cameras[self.obj_name]
            else:
                camera_data = bpy.data.cameras.new(self.obj_name)
            camera = bpy.data.objects.new(self.obj_name, camera_data)
            camera.data.lens = lens
            camera.data.type = camera_type
            camera.location = Vector(location)
            if coll is not None:
                coll.objects.link(camera)