from batoms.magres.magres import Magres

from . import (
    ops,
    ui_list_magres,
)


classes = [
    ops.MagresAdd,
    ops.MagresRemove,
    ops.MagresDraw,
    ui_list_magres.BATOMS_MT_magres_context_menu,
    ui_list_magres.BATOMS_UL_magres,
    ui_list_magres.BATOMS_PT_magres,
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
