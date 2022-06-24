from .lattice_plane import LatticePlane

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.LatticePlaneSetting,
    bpy_data.LatticePlane,
    ops.LatticePlaneAdd,
    ops.LatticePlaneRemove,
    ops.LatticePlaneDraw,
    ops.LatticePlaneModify,
    ui_list.BATOMS_MT_lattice_plane_context_menu,
    ui_list.BATOMS_UL_lattice_plane,
    ui_list.BATOMS_PT_lattice_plane,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.bLatticePlane = PointerProperty(name='bLatticePlane',
                                        type=bpy_data.LatticePlane)
    Object.bLatticePlane = PointerProperty(name='bLatticePlane',
                                    type=bpy_data.LatticePlane)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.bLatticePlane
    del Object.bLatticePlane
