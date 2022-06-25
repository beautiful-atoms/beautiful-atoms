from .bond import Bond

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.BondSetting,
    bpy_data.Bond,
    ops.BondPairAdd,
    ops.BondPairRemove,
    ops.BondModify,
    ops.BondDraw,
    ops.BondOrderAutoSet,
    ops.BondShowHydrogenBond,
    ops.BondShowSearch,
    ui_list.BATOMS_MT_bond_pair_context_menu,
    ui_list.BATOMS_UL_bond_pair,
    ui_list.BATOMS_PT_bond_pair,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bbond = PointerProperty(name='Bbond',
                                        type=bpy_data.Bond)
    Object.Bbond = PointerProperty(name='Bbond',
                                    type=bpy_data.BondSetting)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bbond
    del Object.Bbond
