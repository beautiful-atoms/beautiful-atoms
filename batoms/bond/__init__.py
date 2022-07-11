from .bond import Bond

from . import (
    bpy_data,
    ops,
    gui,
    gui_slicebonds,
)


classes_bpy_data = [
    # internal data first
    bpy_data.BondSetting,
    bpy_data.Bond,
    gui.BondProperties,
    gui_slicebonds.BondProperties,
]

classes = [
    ops.BondPairAdd,
    ops.BondPairRemove,
    ops.BondModify,
    ops.BondDraw,
    ops.BondOrderAutoSet,
    ops.BondShowHydrogenBond,
    ops.BondShowSearch,
    gui.VIEW3D_PT_Batoms_bond,
    gui.BATOMS_MT_bond_pair_context_menu,
    gui.BATOMS_UL_bond_pair,
    gui.BATOMS_PT_bond_pair,
    gui_slicebonds.Bond_PT_prepare,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Bbond = PointerProperty(name='Bbond',
                                        type=bpy_data.Bond)
    Object.Bbond = PointerProperty(name='Bbond',
                                    type=bpy_data.BondSetting)
    Scene.Bbond = PointerProperty(type=gui_slicebonds.BondProperties)
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class

    del Collection.Bbond
    del Object.Bbond
    del Scene.Bbond
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
