from .template import Template

from . import (
    bpy_data,
    ops,
    ui_list_template,
)


classes = [
    bpy_data.TemplateSetting,
    bpy_data.Template,
    ops.TemplateAdd,
    ops.TemplateRemove,
    ops.TemplateDraw,
    ui_list_template.BATOMS_MT_template_context_menu,
    ui_list_template.BATOMS_UL_template,
    ui_list_template.BATOMS_PT_template,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.btemplate = PointerProperty(name='Btemplate',
                                        type=bpy_data.Template)
    Object.btemplate = PointerProperty(name='Btemplate',
                                    type=bpy_data.Template)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.btemplate
    del Object.btemplate
