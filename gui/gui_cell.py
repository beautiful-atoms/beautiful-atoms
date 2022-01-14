import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (BoolProperty,
                       FloatVectorProperty,
                       )

from batoms.butils import get_selected_batoms
from batoms import Batoms

class Cell_PT_prepare(Panel):
    bl_label       = "Cell"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Cell"

  
    def draw(self, context):
        layout = self.layout
        clpanel = context.scene.clpanel

        layout.operator("batoms.apply_cell")
        layout.prop(clpanel, "cell_a")
        layout.prop(clpanel, "cell_b")
        layout.prop(clpanel, "cell_c")
        layout.prop(clpanel, "pbc")

        layout.operator("batoms.apply_transform")
        layout.prop(clpanel, "transform_a")
        layout.prop(clpanel, "transform_b")
        layout.prop(clpanel, "transform_c")
        #
        layout.operator("batoms.apply_boundary")
        layout.prop(clpanel, "boundary_a")
        layout.prop(clpanel, "boundary_b")
        layout.prop(clpanel, "boundary_c")



class CellProperties(bpy.types.PropertyGroup):
    def selected_batoms(self, context):
        if context.object.batoms.batom.flag:
            return context.object.batoms.batom.label
        return None
    def Callback_modify_transform(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = [clpanel.transform_a, clpanel.transform_b, clpanel.transform_c]
        modify_transform(self.selected_batoms(context), transform)
    def Callback_modify_pbc(self, context):
        clpanel = bpy.context.scene.clpanel
        pbc = clpanel.pbc
        modify_batoms_attr(self.selected_batoms(context), 'pbc', pbc)
    def Callback_modify_boundary(self, context):
        clpanel = bpy.context.scene.clpanel
        modify_batoms_attr(self.selected_batoms(context), 'boundary', clpanel.boundary)

    pbc: BoolProperty(
        name = "pbc", default=True,
        description = "pbc", update = Callback_modify_pbc)
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
    transform_a: FloatVectorProperty(
        name="a", default=(1, 0, 0, 0), size = 4,
        # subtype = "XYZ",
        description = "Transform matrix in a axis")
    transform_b: FloatVectorProperty(
        name="b", default=(0, 1, 0, 0), size = 4,
        # subtype = "XYZ",
        description = "Transform matrix in b axis")
    transform_c: FloatVectorProperty(
        name="c", default=(0, 0, 1, 0), size = 4,
        # subtype = "XYZ",
        description = "Transform matrix in c axis")
    
    boundary_a: FloatVectorProperty(
        name="a", default=(0, 1), size = 2,
        description = "boundary  in a axis")
    boundary_b: FloatVectorProperty(
        name="b", default=(0, 1), size = 2,
        description = "boundary  in b axis")
    boundary_c: FloatVectorProperty(
        name="c", default=(0, 1), size = 2,
        description = "boundary  in c axis")

class ApplyCell(Operator):
    bl_idname = "batoms.apply_cell"
    bl_label = "Apply Cell"
    bl_description = ("Apply new cell parameters")
    def execute(self, context):
        clpanel = bpy.context.scene.clpanel
        cell = [clpanel.cell_a, clpanel.cell_b, clpanel.cell_c]
        if not context.object.batoms.batom.flag:
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.batom.label)
        batoms.cell = cell
        return {'FINISHED'}

class ApplyTransform(Operator):
    bl_idname = "batoms.apply_transform"
    bl_label = "Apply Transform"
    bl_description = ("Apply new transform parameters")
    def execute(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = [clpanel.transform_a, clpanel.transform_b, clpanel.transform_c]
        if not context.object.batoms.batom.flag:
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.batom.label)
        batoms.transform(transform)
        return {'FINISHED'}
class ApplyBoundary(Operator):
    bl_idname = "batoms.apply_boundary"
    bl_label = "Apply Boundary"
    bl_description = ("Apply new boundary parameters")
    def execute(self, context):
        clpanel = bpy.context.scene.clpanel
        boundary = [clpanel.boundary_a, clpanel.boundary_b, clpanel.boundary_c]
        if not context.object.batoms.batom.flag:
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.batom.label)
        batoms.boundary = boundary
        return {'FINISHED'}

def modify_batoms_attr(name, key, value):
    if name is None: return
    batoms = Batoms(label = name)
    setattr(batoms, key, value)
    # batoms.select = True
        
def modify_transform(name, transform):
    if name is None: return
    batoms = Batoms(label = name)
    batoms.transform(transform)
    batoms.hide = True
