from batoms.render.render import Render
from batoms.render.light import Light, Lights
from batoms.render.camera import Camera


from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    # internal data first
    bpy_data.LightSetting,
    bpy_data.CameraSetting,
    bpy_data.Render,
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
    Collection.Brender = PointerProperty(name='Brender',
                                        type=bpy_data.Render)
    Object.Blight = PointerProperty(name='Blight',
                                    type=bpy_data.LightSetting)
    Object.Bcamera = PointerProperty(name='Bcamera',
                                    type=bpy_data.CameraSetting)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Brender
    del Object.Blight
    del Object.Bcamera
