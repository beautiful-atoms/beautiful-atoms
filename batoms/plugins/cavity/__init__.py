from .cavity import Cavity

from . import (
    bpy_data,
    ops,
    gui,
)


classes_bpy_data = [
    # internal data first
    bpy_data.CavitySetting,
    bpy_data.Cavity,
    gui.CavityProperties,
]

classes = [
    ops.CavityAdd,
    ops.CavityRemove,
    ops.CavityDraw,
    gui.VIEW3D_PT_Batoms_cavity,
    gui.BATOMS_MT_cavity_context_menu,
    gui.BATOMS_UL_cavity,
    gui.BATOMS_PT_cavity,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Bcavity = PointerProperty(name='Bcavity',
                                        type=bpy_data.Cavity)
    Object.Bcavity = PointerProperty(name='Bcavity',
                                    type=bpy_data.Cavity)
    Scene.Bcavity = PointerProperty(type=gui.CavityProperties)
#
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class

    del Collection.Bcavity
    del Object.Bcavity
    del Scene.Bcavity
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
