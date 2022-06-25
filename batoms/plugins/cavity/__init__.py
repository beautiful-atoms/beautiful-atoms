from .cavity import Cavity

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.CavitySetting,
    bpy_data.Cavity,
    ops.CavityAdd,
    ops.CavityRemove,
    ops.CavityDraw,
    ui_list.BATOMS_MT_cavity_context_menu,
    ui_list.BATOMS_UL_cavity,
    ui_list.BATOMS_PT_cavity,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bcavity = PointerProperty(name='Bcavity',
                                        type=bpy_data.Cavity)
    Object.Bcavity = PointerProperty(name='Bcavity',
                                    type=bpy_data.Cavity)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bcavity
    del Object.Bcavity
