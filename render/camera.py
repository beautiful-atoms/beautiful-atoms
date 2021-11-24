"""
Camera setting

"""
import bpy
from batoms.butils import set_look_at
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
    """Camera Class

    Parameters:

    label: str
    name: str
    coll: str
    type: str
        'ORTHO' or 'PERSP'
    lens: float
    location: array
    look_at : array
    from_coll: bool
    lock_camera_to_view: bool

    """
    
    def __init__(self, label, name = 'Default', 
                        coll = None,
                        type = 'ORTHO', 
                        lens = 39.6, 
                        location = (0, 0, 30), 
                        look_at = None,
                        lock_camera_to_view = True,
                        ) -> None:
        self.label = label
        self.name = name
        obj_name = "%s_camera_%s"%(label, name)
        bobj_name = "blight"
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        self.create_camera(type, location, 
                        lens = lens, coll = coll, look_at = look_at,
                        lock_camera_to_view = lock_camera_to_view)
        self.obj.data.lens_unit = 'FOV'
        
    
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
    
    def create_camera(self, type = None, 
                    location = None, lens = 50, 
                    coll = None, look_at = None,
                    lock_camera_to_view = False):
        '''
        type: str
            * PERSP
            * ORTHO
        lens: float
        '''
        if not self.obj:
            if self.obj_name in bpy.data.cameras:
                camera_data = bpy.data.cameras[self.obj_name]
            else:
                camera_data = bpy.data.cameras.new(self.obj_name)
            camera = bpy.data.objects.new(self.obj_name, camera_data)
            camera.data.lens = lens
            camera.data.type = type
            camera.location = Vector(location)
            if coll is not None:
                coll.objects.link(camera)
            set_lock_camera_to_view(lock_camera_to_view)
            if look_at is not None:
                self.look_at = look_at
    def __repr__(self) -> str:
        s = 'Camera:   type = %s, location = %s, ortho_scale = %s \n' \
                %(self.type, self.location, self.ortho_scale)
        return s
