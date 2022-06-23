from .template import Template

from . import (
    ops,
    ui_list_template,
)


classes = [
    ops.TemplateAdd,
    ops.TemplateRemove,
    ops.TemplateDraw,
    ui_list_template.BATOMS_MT_template_context_menu,
    ui_list_template.BATOMS_UL_template,
    ui_list_template.BATOMS_PT_template,
]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
