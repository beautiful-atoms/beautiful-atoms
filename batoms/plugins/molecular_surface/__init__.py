from .molecular_surface import MolecularSurface

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.MolecularSurfaceSetting,
    bpy_data.MolecularSurface,
    ops.MolecularSurfaceAdd,
    ops.MolecularSurfaceRemove,
    ops.MolecularSurfaceDraw,
    ui_list.BATOMS_MT_molecular_surface_context_menu,
    ui_list.BATOMS_UL_molecular_surface,
    ui_list.BATOMS_PT_molecular_surface,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bmolecularsurface = PointerProperty(name='Bmolecularsurface',
                                        type=bpy_data.MolecularSurface)
    Object.Bmolecularsurface = PointerProperty(name='Bmolecularsurface',
                                    type=bpy_data.MolecularSurface)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bmolecularsurface
    del Object.Bmolecularsurface
