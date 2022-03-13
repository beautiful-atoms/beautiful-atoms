import bpy
import bmesh
from bpy.props import (
    StringProperty,
    BoolProperty,
    FloatVectorProperty,
    FloatProperty,
    IntVectorProperty
)
from batoms import Batoms
from batoms.ops.base import OperatorBatoms, OperatorBatomsEdit


class BatomsReplace(OperatorBatomsEdit):
    bl_idname = "batoms.replace"
    bl_label = "Replace"
    bl_description = "Replace selected atoms by new species"

    species: StringProperty(
        name="species", default='O',
        description="Replaced by this species")

    def execute(self, context):
        obj = context.object
        data = obj.data
        bm = bmesh.from_edit_mesh(data)
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        batoms.replace(v, self.species)
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class BatomModify(OperatorBatomsEdit):
    bl_idname = "batoms.batom_modify"
    bl_label = "Modify batom"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify batom")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this species")

    scale: FloatProperty(name="scale", default=0.6,
                         min=0, soft_max=2, description="scale",
                         )
    bond: BoolProperty(name="Bond", default=True,
                       )
    
    show: BoolProperty(name="Show", default=True,
                       )

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'BATOMS' and obj.mode == 'EDIT'
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
            setattr(batoms[i], self.key, getattr(self, self.key))
        context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode="EDIT")
        return {'FINISHED'}


class ApplyCell(OperatorBatoms):
    bl_idname = "batoms.apply_cell"
    bl_label = "Apply Cell"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new cell parameters")

    a: FloatVectorProperty(
        name="a", default=(1, 0, 0),
        # subtype = "XYZ",
        description="Cell in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1, 0),
        # subtype = "XYZ",
        description="Cell in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 0, 1),
        # subtype = "XYZ",
        description="Cell in c axis")

    def execute(self, context):
        cell = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms.cell = cell
        batoms.obj.select_set(True)
        return {'FINISHED'}


class ApplyTransform(OperatorBatoms):
    bl_idname = "batoms.apply_transform"
    bl_label = "Apply Transform"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new transform parameters")

    a: IntVectorProperty(
        name="a", default=(1, 0, 0, 10), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "MATRIX",
        description="Transform matrix in a axis")
    b: IntVectorProperty(
        name="b", default=(0, 1, 0, 10), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "XYZ",
        description="Transform matrix in b axis")
    c: IntVectorProperty(
        name="c", default=(0, 0, 1, 0), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "XYZ",
        description="Transform matrix in c axis")

    def execute(self, context):
        transform = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms.transform(transform)
        batoms.obj.select_set(True)
        return {'FINISHED'}


class ApplyBoundary(OperatorBatoms):
    bl_idname = "batoms.apply_boundary"
    bl_label = "Apply Boundary"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new boundary parameters")

    a: FloatVectorProperty(
        name="a", default=(0, 1), size=2,
        description="boundary  in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1), size=2,
        description="boundary  in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 1), size=2,
        description="boundary  in c axis")

    def execute(self, context):
        boundary = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms.boundary = boundary
        batoms.obj.select_set(True)
        return {'FINISHED'}
