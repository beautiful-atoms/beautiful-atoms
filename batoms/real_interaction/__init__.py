import bpy
from bpy.props import PointerProperty

from . import (
            rigid_body,
            force_field,
            )


classes = [
    rigid_body.Rigid_Body_Operator,
    rigid_body.Rigid_Body_Modal_Panel,
    rigid_body.RigidBodyProperties,
    force_field.Force_Field_Operator,
    force_field.Force_Field_Modal_Panel,
    force_field.ForceFieldProperties,
    
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    scene = bpy.types.Scene
    scene.rbpanel = PointerProperty(type=rigid_body.RigidBodyProperties)
    scene.ffpanel = PointerProperty(type=force_field.ForceFieldProperties)
    
    
def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    scene = bpy.types.Scene
    del scene.rbpanel
    del scene.ffpanel
    