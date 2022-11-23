import bpy
from bpy.props import (BoolProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class MagresAdd(OperatorBatoms):
    bl_idname = "surface.magres_add"
    bl_label = "Add Magres"
    bl_description = ("Add Magres to a Batoms")

    name: StringProperty(
        name="name", default='2',
        description="Name of Magres to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.magres.settings.add(self.name)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class MagresRemove(OperatorBatoms):
    bl_idname = "surface.magres_remove"
    bl_label = "Remove Magres"
    bl_description = ("Remove Magres to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Magres to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all Magres")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.magres.settings.remove((self.name))
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class MagresDraw(OperatorBatoms):
    bl_idname = "surface.magres_draw"
    bl_label = "Draw Magres"
    bl_description = ("Draw Magres to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of Molecular Surface to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.magres.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}
