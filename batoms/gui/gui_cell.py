import bpy
from bpy.types import Panel
from bpy.props import (BoolProperty,
                       FloatVectorProperty,
                       IntVectorProperty
                       )

from batoms import Batoms


class Cell_PT_prepare(Panel):
    bl_label = "Cell"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Cell"

    def draw(self, context):
        layout = self.layout
        clpanel = context.scene.clpanel

        op = layout.operator("batoms.apply_cell")
        op.a = clpanel.cell_a
        op.b = clpanel.cell_b
        op.c = clpanel.cell_c
        layout.prop(clpanel, "cell_a")
        layout.prop(clpanel, "cell_b")
        layout.prop(clpanel, "cell_c")
        layout.prop(clpanel, "pbc")

        op = layout.operator("batoms.apply_transform")
        op.a = clpanel.transform_a
        op.b = clpanel.transform_b
        op.c = clpanel.transform_c
        layout.prop(clpanel, "transform_a")
        layout.prop(clpanel, "transform_b")
        layout.prop(clpanel, "transform_c")
        #
        op = layout.operator("batoms.apply_boundary")
        op.a = clpanel.boundary_a
        op.b = clpanel.boundary_b
        op.c = clpanel.boundary_c
        layout.prop(clpanel, "boundary_a")
        layout.prop(clpanel, "boundary_b")
        layout.prop(clpanel, "boundary_c")


class CellProperties(bpy.types.PropertyGroup):

    def selected_batoms(self, context):
        if context.object.batoms.batom.flag:
            return context.object.batoms.batom.label
        return None

    def Callback_modify_pbc(self, context):
        clpanel = bpy.context.scene.clpanel
        pbc = clpanel.pbc
        modify_batoms_attr(self.selected_batoms(context), 'pbc', pbc)

    pbc: BoolProperty(
        name="pbc", default=True,
        description="pbc", update=Callback_modify_pbc)
    cell_a: FloatVectorProperty(
        name="a", default=(1, 0, 0),
        # subtype = "XYZ",
        description="Cell in a axis")
    cell_b: FloatVectorProperty(
        name="b", default=(0, 1, 0),
        # subtype = "XYZ",
        description="Cell in b axis")
    cell_c: FloatVectorProperty(
        name="c", default=(0, 0, 1),
        # subtype = "XYZ",
        description="Cell in c axis")
    transform_a: IntVectorProperty(
        name="a", default=(1, 0, 0, 0), size=4,
        min=-20, max=20,
        # subtype = "XYZ",
        description="Transform matrix in a axis")
    transform_b: IntVectorProperty(
        name="b", default=(0, 1, 0, 0), size=4,
        min=-20, max=20,
        # subtype = "XYZ",
        description="Transform matrix in b axis")
    transform_c: IntVectorProperty(
        name="c", default=(0, 0, 1, 0), size=4,
        min=-20, max=20,
        # subtype = "XYZ",
        description="Transform matrix in c axis")

    boundary_a: FloatVectorProperty(
        name="a", default=(0, 1), size=2,
        description="boundary  in a axis")
    boundary_b: FloatVectorProperty(
        name="b", default=(0, 1), size=2,
        description="boundary  in b axis")
    boundary_c: FloatVectorProperty(
        name="c", default=(0, 1), size=2,
        description="boundary  in c axis")


def modify_batoms_attr(name, key, value):
    if name is None:
        return
    batoms = Batoms(label=name)
    setattr(batoms, key, value)
    # batoms.select = True


def modify_transform(name, transform):
    if name is None:
        return
    batoms = Batoms(label=name)
    batoms.transform(transform)
    batoms.hide = True
