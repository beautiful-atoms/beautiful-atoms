import bpy
from bpy.props import (BoolProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class TemplateAdd(OperatorBatoms):
    bl_idname = "surface.template_add"
    bl_label = "Add template"
    bl_description = ("Add template to a Batoms")

    name: StringProperty(
        name="name", default='2',
        description="Name of template to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.template.setting.add(self.name)
        batoms.coll.batoms.template_index = len(batoms.template.setting) - 1
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class TemplateRemove(OperatorBatoms):
    bl_idname = "surface.template_remove"
    bl_label = "Remove template"
    bl_description = ("Remove template to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of template to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all template")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.batoms.template_index
        batoms.template.setting.remove((self.name))
        batoms.coll.batoms.template_index = min(max(0, index - 1),
                                          len(batoms.template.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class TemplateDraw(OperatorBatoms):
    bl_idname = "surface.template_draw"
    bl_label = "Draw template"
    bl_description = ("Draw template to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of template surface to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.template.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

