import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy_extras.object_utils import AddObjectHelper
from bpy.props import (BoolProperty,
                       FloatVectorProperty,
                       IntVectorProperty
                       )
from batoms import Batoms

class ApplyCell(Operator, AddObjectHelper):
    bl_idname = "batoms.apply_cell"
    bl_label = "Apply Cell"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new cell parameters")

    a: FloatVectorProperty(
        name="a", default=(1, 0, 0),
        # subtype = "XYZ",
        description = "Cell in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1, 0),
        # subtype = "XYZ",
        description = "Cell in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 0, 1),
        # subtype = "XYZ",
        description = "Cell in c axis")

    # def draw(self, context):
    #     layout = self.layout
    #     box = layout.box()
    #     box.label(text="Basic Parameters")
    #     box.prop(self, 'a')
    #     box.prop(self, 'b')
    #     box.prop(self, 'c')

    def execute(self, context):
        cell = [self.a, self.b, self.c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.cell = cell
        return {'FINISHED'}

class ApplyTransform(Operator, AddObjectHelper):
    bl_idname = "batoms.apply_transform"
    bl_label = "Apply Transform"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new transform parameters")

    a: IntVectorProperty(
        name="a", default=(1, 0, 0, 10), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "MATRIX",
        description = "Transform matrix in a axis")
    b: IntVectorProperty(
        name="b", default=(0, 1, 0, 10), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "XYZ",
        description = "Transform matrix in b axis")
    c: IntVectorProperty(
        name="c", default=(0, 0, 1, 0), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "XYZ",
        description = "Transform matrix in c axis")

    def execute(self, context):
        transform = [self.a, self.b, self.c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.transform(transform)

        return {'FINISHED'}


class ApplyBoundary(Operator, AddObjectHelper):
    bl_idname = "batoms.apply_boundary"
    bl_label = "Apply Boundary"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new boundary parameters")

    a: FloatVectorProperty(
        name="a", default=(0, 1), size = 2,
        description = "boundary  in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1), size = 2,
        description = "boundary  in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 1), size = 2,
        description = "boundary  in c axis")
        
    def execute(self, context):
        boundary = [self.a, self.b, self.c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.boundary = boundary
        return {'FINISHED'}