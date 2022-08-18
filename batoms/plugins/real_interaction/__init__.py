import bpy
from bpy.props import PointerProperty

from . import (
            modal_force_field_ase,
            modal_force_field_openbabel,
            modal_rigid_body,
            )


classes = [
    modal_rigid_body.Rigid_Body_Operator,
    modal_rigid_body.Rigid_Body_Modal_Panel,
    modal_rigid_body.RigidBodyProperties,
    modal_force_field_ase.ASE_Force_Field_Operator,
    modal_force_field_ase.BATOMS_PT_ASE_Force_Field_Modal,
    modal_force_field_ase.ASEForceFieldProperties,
    modal_force_field_openbabel.OB_Force_Field_Operator,
    modal_force_field_openbabel.BATOMS_PT_OB_Force_Field_Modal,
    modal_force_field_openbabel.OBForceFieldProperties,
    
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    scene = bpy.types.Scene
    scene.rbpanel = PointerProperty(type=modal_rigid_body.RigidBodyProperties)
    scene.aseffpanel = PointerProperty(type=modal_force_field_ase.ASEForceFieldProperties)
    scene.obffpanel = PointerProperty(type=modal_force_field_openbabel.OBForceFieldProperties)
    
    
def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    scene = bpy.types.Scene
    del scene.rbpanel
    del scene.aseffpanel
    del scene.obffpanel
    