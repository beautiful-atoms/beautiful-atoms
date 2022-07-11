"""
Camera setting

"""
import bpy
from mathutils import Vector
from batoms.base.object import BaseObject
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def set_lock_camera_to_view(switch):
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.lock_camera = switch


class Camera(BaseObject):
    def __init__(self, label, name='Default',
                 coll=None,
                 type='ORTHO',
                 lens=39.6,
                 location=(0, 0, 30),
                 look_at=None,
                 lock_camera_to_view=True,
                 ) -> None:
        """Camera Class

        Args:
            label (_type_): _description_
            name (str, optional): _description_. Defaults to 'Default'.
            coll (_type_, optional): _description_. Defaults to None.
            type (str, optional): _description_. Defaults to 'ORTHO'. 'ORTHO' or 'PERSP'
            lens (float, optional): _description_. Defaults to 39.6.
            location (tuple, optional): _description_. Defaults to (0, 0, 30).
            look_at (_type_, optional): _description_. Defaults to None.
            lock_camera_to_view (bool, optional): _description_. Defaults to True.
        """
        self.label = label
        self.name = name
        obj_name = "%s_camera_%s" % (label, name)
        btype = "Bcamera"
        BaseObject.__init__(self, obj_name=obj_name, btype=btype)
        self.create_camera(type, location,
                           lens=lens, coll=coll, look_at=look_at,
                           lock_camera_to_view=lock_camera_to_view)
        self.obj.data.lens_unit = 'FOV'

    def get_bpy_data(self):
        bpy_data = getattr(self.obj, self.btype)
        return bpy_data

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

    def create_camera(self, type=None,
                      location=None, lens=50,
                      coll=None, look_at=None,
                      lock_camera_to_view=False):
        '''
        type: str
            * PERSP
            * ORTHO
        lens: float
        '''
        if self.obj_name in bpy.data.cameras:
            camera_data = bpy.data.cameras[self.obj_name]
        else:
            camera_data = bpy.data.cameras.new(self.obj_name)
        camera = bpy.data.objects.get(self.obj_name)
        # create camera if not exist
        if not camera:
            # bpy.data.objects.remove(obj, do_unlink=True)
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
            % (self.type, self.location, self.ortho_scale)
        return s
