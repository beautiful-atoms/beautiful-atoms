from . import gui
from .lattice_plane import LatticePlane

from . import (
    bpy_data,
    ops,
    gui,
)


classes_bpy_data = [
    # internal data first
    bpy_data.LatticePlaneSetting,
    bpy_data.LatticePlane,
    gui.LatticePlaneProperties,
]
classes = [
    ops.LatticePlaneAdd,
    ops.LatticePlaneRemove,
    ops.LatticePlaneDraw,
    ops.LatticePlaneModify,
    gui.VIEW3D_PT_Batoms_lattice_plane,
    gui.BATOMS_MT_lattice_plane_context_menu,
    gui.BATOMS_UL_lattice_plane,
    gui.BATOMS_PT_lattice_plane,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    #
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Blatticeplane = PointerProperty(name='Blatticeplane',
                                        type=bpy_data.LatticePlane)
    Object.Blatticeplane = PointerProperty(name='Blatticeplane',
                                    type=bpy_data.LatticePlane)
    Scene.Blatticeplane = PointerProperty(type=gui.LatticePlaneProperties)
    #
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class


    del Collection.Blatticeplane
    del Object.Blatticeplane
    del Scene.Blatticeplane
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
