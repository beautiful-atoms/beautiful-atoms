from batoms.cavity.cavity import Cavity

from . import (
    ops,
    ui_list_cavity,
)


classes = [
    ops.CavityAdd,
    ops.CavityRemove,
    ops.CavityDraw,
    ui_list_cavity.BATOMS_MT_cavity_context_menu,
    ui_list_cavity.BATOMS_UL_cavity,
    ui_list_cavity.BATOMS_PT_cavity,
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
