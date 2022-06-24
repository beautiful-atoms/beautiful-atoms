from .magres import Magres

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    bpy_data.MagresSetting,
    bpy_data.Magres,
    ops.MagresAdd,
    ops.MagresRemove,
    ops.MagresDraw,
    ui_list.BATOMS_MT_magres_context_menu,
    ui_list.BATOMS_UL_magres,
    ui_list.BATOMS_PT_magres,
]

def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bmagres = PointerProperty(name='Bmagres',
                                        type=bpy_data.Magres)
    Object.Bmagres = PointerProperty(name='Bmagres',
                                    type=bpy_data.MagresSetting)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bmagres
    del Object.Bmagres
