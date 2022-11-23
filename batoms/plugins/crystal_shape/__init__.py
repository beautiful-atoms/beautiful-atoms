from .crystal_shape import CrystalShape

from . import (
    bpy_data,
    ops,
    gui,
)


classes_bpy_data = [
    # internal data first
    bpy_data.CrystalShapeSetting,
    bpy_data.CrystalShape,
    gui.CrystalShapeProperties,
]
classes = [
    ops.CrystalShapeAdd,
    ops.CrystalShapeRemove,
    ops.CrystalShapeDraw,
    ops.CrystalShapeModify,
    gui.VIEW3D_PT_Batoms_crystal_shape,
    gui.BATOMS_MT_crystal_shape_context_menu,
    gui.BATOMS_UL_crystal_shape,
    gui.BATOMS_PT_crystal_shape,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    #
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Bcrystalshape = PointerProperty(name='Bcrystalshape',
                                        type=bpy_data.CrystalShape)
    Object.Bcrystalshape = PointerProperty(name='Bcrystalshape',
                                    type=bpy_data.CrystalShape)
    Scene.Bcrystalshape = PointerProperty(type=gui.CrystalShapeProperties)
    #
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class

    del Collection.Bcrystalshape
    del Object.Bcrystalshape
    del Scene.Bcrystalshape
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
