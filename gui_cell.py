import bpy
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

from batoms.butils import get_selected_batoms
from batoms import Batoms

class Cell_PT_prepare(Panel):
    bl_label       = "Cell"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Cell"
    bl_idname = "BATOMS_PT_Cell"

  
    def draw(self, context):
        layout = self.layout
        clpanel = context.scene.clpanel

        box = layout.box()
        row = box.row()
        row.prop(clpanel, "pbc")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Cell")
        row = box.row()
        row.prop(clpanel, "cell_a")
        row.prop(clpanel, "cell_b")
        row.prop(clpanel, "cell_c")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="SuperCell")
        row = box.row()
        row.prop(clpanel, "supercell_a")
        row.prop(clpanel, "supercell_b")
        row.prop(clpanel, "supercell_c")
        #
        box = layout.box()
        col = box.column()
        col.prop(clpanel, "boundary")



class CellProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    def Callback_modify_supercell(self, context):
        clpanel = bpy.context.scene.clpanel
        supercell = [clpanel.supercell_a, clpanel.supercell_b, clpanel.supercell_c]
        modify_supercell(self.selected_batoms, supercell)
    def Callback_modify_cell(self, context):
        clpanel = bpy.context.scene.clpanel
        cell = [clpanel.cell_a, clpanel.cell_b, clpanel.cell_c]
        modify_batoms_attr(self.selected_batoms, 'cell', cell)
    def Callback_modify_pbc(self, context):
        clpanel = bpy.context.scene.clpanel
        pbc = clpanel.pbc
        modify_batoms_attr(self.selected_batoms, 'pbc', pbc)
    def Callback_modify_boundary(self, context):
        clpanel = bpy.context.scene.clpanel
        modify_batoms_attr(self.selected_batoms, 'boundary', clpanel.boundary)

    pbc: BoolProperty(
        name = "pbc", default=True,
        description = "pbc", update = Callback_modify_pbc)
    cell_a: FloatProperty(
        name = "a",
        description = "cell a")
    cell_b: FloatProperty(
        name = "b",
        description = "cell b")
    cell_c: FloatProperty(
        name = "c",
        description = "cell c", update = Callback_modify_cell)
    supercell_a: IntProperty(
        name = "a", default=1,
        description = "cell a")
    supercell_b: IntProperty(
        name = "b", default=1,
        description = "cell b")
    supercell_c: IntProperty(
        name = "c", default=1,
        description = "cell c", update = Callback_modify_supercell)
    boundary: FloatVectorProperty(
        name="boundary", default=(0.00, 0.0, 0.0),
        subtype = "XYZ",
        description = "boundary  in a, b, c axis", update = Callback_modify_boundary)
    
def modify_batoms_attr(batoms_name_list, key, value):
    batoms_list = []
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        setattr(batoms, key, value)
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True
        
def modify_supercell(batoms_name_list, supercell):
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        batoms.repeat(supercell)
