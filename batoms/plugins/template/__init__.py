from .template import Template

from . import (
    bpy_data,
    ops,
    ui_list,
)


classes = [
    bpy_data.TemplateSetting,
    bpy_data.Template,
    ops.TemplateAdd,
    ops.TemplateRemove,
    ops.TemplateDraw,
    ui_list.BATOMS_MT_template_context_menu,
    ui_list.BATOMS_UL_template,
    ui_list.BATOMS_PT_template,
]


def register_class():
    from bpy.types import Collection, Object
    from bpy.props import PointerProperty
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # attach to blender internal data
    Collection.Btemplate = PointerProperty(name='Btemplate',
                                        type=bpy_data.Template)
    Object.Btemplate = PointerProperty(name='Btemplate',
                                    type=bpy_data.Template)


def unregister_class():
    from bpy.types import Collection, Object
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.Btemplate
    del Object.Btemplate
