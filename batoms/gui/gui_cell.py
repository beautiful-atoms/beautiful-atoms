import bpy
from bpy.types import Panel
from bpy.props import (BoolProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       IntVectorProperty
                       )

from batoms import Batoms
from batoms.cell import Bcell
from batoms.gui.gui_batoms import get_active_batoms

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

        # op = layout.operator("batoms.apply_cell")
        # op.a = [clpanel.cell_a0, clpanel.cell_a1, clpanel.cell_a2]
        # op.b = [clpanel.cell_b0, clpanel.cell_b1, clpanel.cell_b2]
        # op.c = [clpanel.cell_c0, clpanel.cell_c1, clpanel.cell_c2]
        box = layout.box()
        row = box.row(align = True)
        row.prop(clpanel, "cell_a0")
        row.prop(clpanel, "cell_a1")
        row.prop(clpanel, "cell_a2")
        row = box.row(align = True)
        row.prop(clpanel, "cell_b0")
        row.prop(clpanel, "cell_b1")
        row.prop(clpanel, "cell_b2")
        row = box.row(align = True)
        row.prop(clpanel, "cell_c0")
        row.prop(clpanel, "cell_c1")
        row.prop(clpanel, "cell_c2")
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

def get_cell(i1, i2):
    """Helper function to easily get cell-property."""

    def getter(self):
        cell = get_active_cell()
        if cell is not None:
            return cell.array[i1, i2]
        else:
            return 0
    return getter

def set_cell(i1, i2):
    """Helper function to easily set cell-property."""

    def setter(self, value):
        cell = get_active_cell()
        if cell is not None:
            cell[i1, i2] = value
            
    return setter


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
    # getter/setter only used for float not float vector
    cell_a0: FloatProperty(
        name="a0", default=0,
        description="Cell in a axis",
        get=get_cell(0, 0),
        set=set_cell(0, 0))
    cell_a1: FloatProperty(
        name="a1", default=0,
        description="Cell in a axis",
        get=get_cell(0, 1),
        set=set_cell(0, 1))
    cell_a2: FloatProperty(
        name="a2", default=0,
        description="Cell in a axis",
        get=get_cell(0, 2),
        set=set_cell(0, 2))
    cell_b0: FloatProperty(
        name="b0", default=0,
        description="Cell in a axis",
        get=get_cell(1, 0),
        set=set_cell(1, 0))
    cell_b1: FloatProperty(
        name="b1", default=0,
        description="Cell in a axis",
        get=get_cell(1, 1),
        set=set_cell(1, 1))
    cell_b2: FloatProperty(
        name="b2", default=0,
        description="Cell in a axis",
        get=get_cell(1, 2),
        set=set_cell(1, 2))
    cell_c0: FloatProperty(
        name="c0", default=0,
        description="Cell in a axis",
        get=get_cell(2, 0),
        set=set_cell(2, 0))
    cell_c1: FloatProperty(
        name="c1", default=0,
        description="Cell in a axis",
        get=get_cell(2, 1),
        set=set_cell(2, 1))
    cell_c2: FloatProperty(
        name="c2", default=0,
        description="Cell in a axis",
        get=get_cell(2, 2),
        set=set_cell(2, 2))
    
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

def get_active_cell():
    context = bpy.context
    if context.object and context.object.batoms.type != 'OTHER':
        cell = Bcell(label=context.object.batoms.label)
        return cell
    return None