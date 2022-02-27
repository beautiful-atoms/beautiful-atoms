import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (BoolProperty,
                       FloatVectorProperty,
                       IntVectorProperty
                       )
from batoms import Batoms

class ApplyCell(Operator):
    bl_idname = "batoms.apply_cell"
    bl_label = "Apply Cell"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new cell parameters")

    cell_a: FloatVectorProperty(
        name="a", default=(1, 0, 0),
        # subtype = "XYZ",
        description = "Cell in a axis")
    cell_b: FloatVectorProperty(
        name="b", default=(0, 1, 0),
        # subtype = "XYZ",
        description = "Cell in b axis")
    cell_c: FloatVectorProperty(
        name="c", default=(0, 0, 1),
        # subtype = "XYZ",
        description = "Cell in c axis")

    def execute(self, context):
        cell = [self.cell_a, self.cell_b, self.cell_c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.cell = cell
        return {'FINISHED'}

class ApplyTransform(Operator):
    bl_idname = "batoms.apply_transform"
    bl_label = "Apply Transform"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new transform parameters")

    transform_a: IntVectorProperty(
        name="a", default=(1, 0, 0, 0), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "XYZ",
        description = "Transform matrix in a axis")
    transform_b: IntVectorProperty(
        name="b", default=(0, 1, 0, 0), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "XYZ",
        description = "Transform matrix in b axis")
    transform_c: IntVectorProperty(
        name="c", default=(0, 0, 1, 0), size = 4,
        soft_min = -10, soft_max = 10,
        # subtype = "XYZ",
        description = "Transform matrix in c axis")

    def execute(self, context):
        transform = [self.transform_a, self.transform_b, self.transform_c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.transform(transform)
        return {'FINISHED'}

class ApplyBoundary(Operator):
    bl_idname = "batoms.apply_boundary"
    bl_label = "Apply Boundary"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new boundary parameters")

    boundary_a: FloatVectorProperty(
        name="a", default=(0, 1), size = 2,
        description = "boundary  in a axis")
    boundary_b: FloatVectorProperty(
        name="b", default=(0, 1), size = 2,
        description = "boundary  in b axis")
    boundary_c: FloatVectorProperty(
        name="c", default=(0, 1), size = 2,
        description = "boundary  in c axis")
        
    def execute(self, context):
        boundary = [self.boundary_a, self.boundary_b, self.boundary_c]
        if not context.object.batoms.flag:
            print('Please select a Batoms object.')
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms.boundary = boundary
        return {'FINISHED'}