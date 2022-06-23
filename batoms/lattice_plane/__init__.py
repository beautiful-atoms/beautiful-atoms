from batoms.lattice_plane.lattice_plane import LatticePlane

from . import (
    ops,
    ui_list_lattice_plane,
)


classes = [
    ops.LatticePlaneAdd,
    ops.LatticePlaneRemove,
    ops.LatticePlaneDraw,
    ops.LatticePlaneModify,
    ui_list_lattice_plane.BATOMS_MT_lattice_plane_context_menu,
    ui_list_lattice_plane.BATOMS_UL_lattice_plane,
    ui_list_lattice_plane.BATOMS_PT_lattice_plane,
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
