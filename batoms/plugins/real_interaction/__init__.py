import bpy
from bpy.props import PointerProperty

from . import (
            modal_rigid_body,
            modal_force_field,
            )


classes = [
    modal_rigid_body.Rigid_Body_Operator,
    modal_rigid_body.Rigid_Body_Modal_Panel,
    modal_rigid_body.RigidBodyProperties,
    modal_force_field.Force_Field_Operator,
    modal_force_field.Force_Field_Modal_Panel,
    modal_force_field.ForceFieldProperties,

]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    scene = bpy.types.Scene
    scene.rbpanel = PointerProperty(type=modal_rigid_body.RigidBodyProperties)
    scene.ffpanel = PointerProperty(type=modal_force_field.ForceFieldProperties)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    scene = bpy.types.Scene
    del scene.rbpanel
    del scene.ffpanel
