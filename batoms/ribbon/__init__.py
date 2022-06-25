from .ribbon import Ribbon

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.SheetSetting,
    bpy_data.HelixSetting,
    bpy_data.TurnSetting,
    bpy_data.Protein,
    # ops.ProteinAdd,
    # ops.ProteinRemove,
    # ops.ProteinDraw,
    # ui_list.BATOMS_MT_sheet_context_menu,
    # ui_list.BATOMS_UL_sheet,
    # ui_list.BATOMS_PT_sheet,
    # ui_list.BATOMS_MT_helix_context_menu,
    # ui_list.BATOMS_UL_helix,
    # ui_list.BATOMS_PT_helix,
    # ui_list.BATOMS_MT_turn_context_menu,
    # ui_list.BATOMS_UL_turn,
    # ui_list.BATOMS_PT_turn,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Bprotein = PointerProperty(name='Bprotein',
                                        type=bpy_data.Protein)
    Object.Bprotein = PointerProperty(name='Bprotein',
                                    type=bpy_data.Protein)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Bprotein
    del Object.Bprotein
