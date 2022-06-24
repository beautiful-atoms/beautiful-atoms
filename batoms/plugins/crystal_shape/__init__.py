from .crystal_shape import CrystalShape

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.CrystalShapeSetting,
    bpy_data.CrystalShape,
    ops.CrystalShapeAdd,
    ops.CrystalShapeRemove,
    ops.CrystalShapeDraw,
    ops.CrystalShapeModify,
    ui_list.BATOMS_MT_crystal_shape_context_menu,
    ui_list.BATOMS_UL_crystal_shape,
    ui_list.BATOMS_PT_crystal_shape,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.bCrystalShape = PointerProperty(name='bCrystalShape',
                                        type=bpy_data.CrystalShape)
    Object.bCrystalShape = PointerProperty(name='bCrystalShape',
                                    type=bpy_data.CrystalShape)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.bCrystalShape
    del Object.bCrystalShape
