from .isosurface import Isosurface

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.IsosurfaceSetting,
    bpy_data.Isosurface,
    ops.IsosurfaceAdd,
    ops.IsosurfaceRemove,
    ops.IsosurfaceDraw,
    ui_list.BATOMS_MT_isosurface_context_menu,
    ui_list.BATOMS_UL_isosurface,
    ui_list.BATOMS_PT_isosurface,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bisosurface = PointerProperty(name='Bisosurface',
                                        type=bpy_data.Isosurface)
    Object.Bisosurface = PointerProperty(name='Bisosurface',
                                    type=bpy_data.Isosurface)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bisosurface
    del Object.Bisosurface
