from .template import Template

from . import (
    bpy_data,
    ops,
    gui,
)


classes_bpy_data = [
    bpy_data.TemplateSetting,
    bpy_data.Template,
    gui.TemplateProperties,
]
classes = [
    ops.TemplateAdd,
    ops.TemplateRemove,
    ops.TemplateDraw,
    gui.VIEW3D_PT_Batoms_template,
    gui.BATOMS_MT_template_context_menu,
    gui.BATOMS_UL_template,
    gui.BATOMS_PT_template,
]


def register_class():
    from bpy.types import Collection, Object, Scene
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes_bpy_data:
        register_class(cls)
    # attach to blender internal data
    Collection.Btemplate = PointerProperty(name='Btemplate',
                                        type=bpy_data.Template)
    Object.Btemplate = PointerProperty(name='Btemplate',
                                    type=bpy_data.Template)
    Scene.Btemplate = PointerProperty(type=gui.TemplateProperties)
    #
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.types import Collection, Object, Scene
    from bpy.utils import unregister_class
    #
    del Collection.Btemplate
    del Object.Btemplate
    del Scene.Btemplate
    for cls in reversed(classes):
        unregister_class(cls)
    for cls in reversed(classes_bpy_data):
        unregister_class(cls)
