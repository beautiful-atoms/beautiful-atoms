from .crystal_shape import CrystalShape

from . import (
    ops,
    ui_list_crystal_shape,
)


classes = [
    ops.CrystalShapeAdd,
    ops.CrystalShapeRemove,
    ops.CrystalShapeDraw,
    ops.CrystalShapeModify,
    ui_list_crystal_shape.BATOMS_MT_crystal_shape_context_menu,
    ui_list_crystal_shape.BATOMS_UL_crystal_shape,
    ui_list_crystal_shape.BATOMS_PT_crystal_shape,
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
