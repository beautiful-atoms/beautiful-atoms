import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms

class HighlightAdd(OperatorBatoms):
    bl_idname = "batoms.highlight_add"
    bl_label = "Add highlight"
    bl_description = ("Add highlight to a Batoms")

    name: StringProperty(
        name="name", default='0',
        description="Name of highlight to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.highlight.settings.add(self.name)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class HighlightRemove(OperatorBatoms):
    bl_idname = "batoms.highlight_remove"
    bl_label = "Remove highlight"
    bl_description = ("Remove highlight to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of highlight to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all highlights")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.highlight.settings.remove((self.name))
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class HighlightDraw(OperatorBatoms):
    bl_idname = "batoms.highlight_draw"
    bl_label = "Draw highlight"
    bl_description = ("Draw highlight to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of highlight to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.highlight.draw()
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}
