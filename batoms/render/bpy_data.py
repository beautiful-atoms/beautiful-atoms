import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       CollectionProperty,
                       )

from batoms.internal_data import Base


class LightSetting(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    name: StringProperty(name="name", default='X')
    type: StringProperty(name="type", default='SUN')
    strength: FloatProperty(name="strength", soft_min=0.01, soft_max=1, default=0.01)
    lock_to_camera: BoolProperty(name="lock_to_camera", default=False)
    direction: FloatVectorProperty(name="direction", default=[0, 0, 1], size=3)
    look_at: FloatVectorProperty(name="look_at", default=[0, 0, 0], size=3)


class CameraSetting(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    name: StringProperty(name="name", default='X')
    type: StringProperty(name="type", default='ORTHO')
    direction: FloatVectorProperty(name="direction", default=[0, 0, 1], size=3)
    look_at: FloatVectorProperty(name="look_at", default=[0, 0, 0], size=3)


class Render(bpy.types.PropertyGroup):
    """This module defines the Render properties to extend
    Blender’s internal data.

    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    engine: StringProperty(name="engine", default='BLENDER_EEVEE')
    compute_device_type: StringProperty(
        name="compute_device_type", default='CUDA')
    animation: BoolProperty(name="animation", default=False)
    run_render: BoolProperty(name="run_render", default=True)
    gpu: BoolProperty(name="gpu", default=False)
    viewport: FloatVectorProperty(name="viewport", default=[0, 0, 1], size=3)
    center: FloatVectorProperty(name="center", default=[0, 0, 1], size=3)
    distance: FloatProperty(name="distance",
                            description="Distance from camera",
                            default=-1)
    padding: FloatVectorProperty(name="padding", default=[1, 1, 1, 1], size=4)


    light_ui_list_index: IntProperty(name="light_ui_list_index",
                              default=0)
    camera_ui_list_index: IntProperty(name="camera_ui_list_index",
                              default=0)
    # colleciton properties
    light_settings: CollectionProperty(name='lightsetting',
                                type=LightSetting)
    camera_settings: CollectionProperty(name='camerasetting',
                                type=CameraSetting)
