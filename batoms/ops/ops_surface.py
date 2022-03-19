import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class MSAdd(OperatorBatoms):
    bl_idname = "surface.ms_add"
    bl_label = "Add Molecular Surface"
    bl_description = ("Add Molecular Surface to a Batoms")

    name: StringProperty(
        name="name", default='2',
        description="Name of Molecular Surface to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.ms.setting.add(self.name)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class MSRemove(OperatorBatoms):
    bl_idname = "surface.ms_remove"
    bl_label = "Remove Molecular Surface"
    bl_description = ("Remove Molecular Surface to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Molecular Surface to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all Molecular Surfaces")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.batoms.ms_index
        batoms.ms.setting.remove((self.name))
        batoms.coll.batoms.ms_index = min(max(0, index - 1),
                                          len(batoms.ms.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class MSDraw(OperatorBatoms):
    bl_idname = "surface.ms_draw"
    bl_label = "Draw Molecular Surface"
    bl_description = ("Draw Molecular Surface to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of Molecular Surface to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.ms.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class MSModify(OperatorBatoms):
    bl_idname = "surface.ms_modify"
    bl_label = "Modify ms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify ms")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this species")

    slice: BoolProperty(name="slice", default=False,
                        )
    boundary: BoolProperty(name="boundary", default=False,
                           )
    distance: FloatProperty(name="distance",
                            description="Distance from origin",
                            default=1)

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'MS' and obj.mode == 'EDIT'
        else:
            return False

    def execute(self, context):
        obj = context.object
        data = obj.data
        bm = bmesh.from_edit_mesh(data)
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        for i in v:
            setattr(batoms.bonds[i], self.key, getattr(self, self.key))
        # batoms.draw()
        return {'FINISHED'}


class IsosurfaceAdd(OperatorBatoms):
    bl_idname = "surface.isosurface_add"
    bl_label = "Add Isosurface"
    bl_description = ("Add Isosurface to a Batoms")

    name: StringProperty(
        name="name", default='2',
        description="Name of Isosurface to be added")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.isosurfaces.setting.add(self.name)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class IsosurfaceRemove(OperatorBatoms):
    bl_idname = "surface.isosurface_remove"
    bl_label = "Remove Isosurface"
    bl_description = ("Remove Isosurface to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Isosurface to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all Isosurfaces")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.batoms.isosurface_index
        batoms.isosurfaces.setting.remove((self.name))
        batoms.coll.batoms.isosurface_index = min(max(0, index - 1),
                                                  len(batoms.isosurfaces.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class IsosurfaceDraw(OperatorBatoms):
    bl_idname = "surface.isosurface_draw"
    bl_label = "Draw Isosurface"
    bl_description = ("Draw Isosurface to a Batoms")

    name: StringProperty(
        name="name", default='ALL',
        description="Name of Isosurface to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.isosurfaces.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class IsosurfaceModify(Operator):
    bl_idname = "surface.isosurface_modify"
    bl_label = "Modify Isosurface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify Isosurface")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this species")

    slice: BoolProperty(name="slice", default=False,
                        )
    boundary: BoolProperty(name="boundary", default=False,
                           )
    distance: FloatProperty(name="distance",
                            description="Distance from origin",
                            default=1)

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'MS' and obj.mode == 'EDIT'
        else:
            return False

    def execute(self, context):
        obj = context.object
        data = obj.data
        bm = bmesh.from_edit_mesh(data)
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        for i in v:
            setattr(batoms.bonds[i], self.key, getattr(self, self.key))
        # batoms.draw()
        return {'FINISHED'}
