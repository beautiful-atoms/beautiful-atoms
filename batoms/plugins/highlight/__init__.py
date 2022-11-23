from .highlight import Highlight

from . import (
    bpy_data,
    ops,
    gui,
)


classes_bpy_data = [
    # internal data first
    bpy_data.HighlightSetting,
    bpy_data.Highlight,
    gui.HighlightProperties,
]

classes = [
    ops.HighlightAdd,
    ops.HighlightRemove,
    ops.HighlightDraw,
    gui.VIEW3D_PT_Batoms_highlight,
    gui.BATOMS_MT_highlight_context_menu,
    gui.BATOMS_UL_highlight,
    gui.BATOMS_PT_highlight,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Bhighlight = PointerProperty(name='Bhighlight',
                                        type=bpy_data.Highlight)
    Object.Bhighlight = PointerProperty(name='Bhighlight',
                                    type=bpy_data.Highlight)
    Scene.Bhighlight = PointerProperty(type=gui.HighlightProperties)
#
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class

    del Collection.Bhighlight
    del Object.Bhighlight
    del Scene.Bhighlight
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
