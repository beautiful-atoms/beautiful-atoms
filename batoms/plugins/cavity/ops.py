import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms

class CavityAdd(OperatorBatoms):
    bl_idname = "surface.cavity_add"
    bl_label = "Add Cavity"
    bl_description = ("Add cavity to a Batoms")

    name: StringProperty(
        name="name", default='0',
        description="Name of cavity to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.cavity.settings.add(self.name)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class CavityRemove(OperatorBatoms):
    bl_idname = "surface.cavity_remove"
    bl_label = "Remove Cavity"
    bl_description = ("Remove cavity to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of cavity to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all cavitys")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.cavity.settings.remove((self.name))
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class CavityDraw(OperatorBatoms):
    bl_idname = "surface.cavity_draw"
    bl_label = "Draw Cavity"
    bl_description = ("Draw cavity to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of cavity to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.cavity.draw()
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}
