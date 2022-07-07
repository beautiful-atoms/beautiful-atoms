from .polyhedra import Polyhedra

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.PolyhedraSetting,
    bpy_data.Polyhedra,
    ops.PolyhedraAdd,
    ops.PolyhedraRemove,
    ops.PolyhedraDraw,
    ops.PolyhedraModify,
    ui_list.BATOMS_MT_polyhedra_context_menu,
    ui_list.BATOMS_UL_polyhedra,
    ui_list.BATOMS_PT_polyhedra,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bpolyhedra = PointerProperty(name='Bpolyhedra',
                                        type=bpy_data.Polyhedra)
    Object.Bpolyhedra = PointerProperty(name='Bpolyhedra',
                                    type=bpy_data.PolyhedraSetting)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bpolyhedra
    del Object.Bpolyhedra
