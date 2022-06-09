import bpy
from bpy.types import Panel
from bpy.props import (BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       )


# The panel.
class Bond_PT_prepare(Panel):
    bl_label = "Bond"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Bond"

    def draw(self, context):
        layout = self.layout
        bond = context.scene.batoms.bond

        layout.operator("bond.bond_order_auto_set", icon='MODIFIER_ON', text="Auto Set Bond Order")
        layout.separator()
        layout.label(text="Bond style")
        col = layout.column()
        col.prop(bond, "bond_style", expand=True)
        layout.prop(bond, "bondwidth")
        layout.prop(bond, "order")
        layout.prop(bond, "show")


class BondProperties(bpy.types.PropertyGroup):
    def Callback_bond_style(self, context):
        bond = bpy.context.scene.batoms.bond
        bond_style = list(bond.bond_style)[0]
        bpy.ops.bond.bond_modify(key='style', style=bond_style)

    def Callback_modify_bondwidth(self, context):
        bond = bpy.context.scene.batoms.bond
        bondwidth = bond.bondwidth
        bpy.ops.bond.bond_modify(key='width', bondwidth=bondwidth)

    def Callback_modify_show(self, context):
        bond = bpy.context.scene.batoms.bond
        show = bond.show
        bpy.ops.bond.bond_modify(key='show', show=show)

    def Callback_modify_order(self, context):
        bond = bpy.context.scene.batoms.bond
        order = bond.order
        bpy.ops.bond.bond_modify(key='order', order=order)

    bond_style: EnumProperty(
        name="style",
        description="bond style",
        items=(('0', "Unicolor cylinder", ""),
               ('1', "Bicolor cylinder", ""),
               ('2', "Dashed line", ""),
               ('3', "Spring", "")),
        default={'1'},
        update=Callback_bond_style,
        options={'ENUM_FLAG'},
    )
    bondwidth: FloatProperty(
        name="bondwidth", default=0.1,
        min=0, soft_max=1,
        description="bondwidth", update=Callback_modify_bondwidth)
    order: IntProperty(name="Bond order", default=1,
                       min=1, max=3,
                       update=Callback_modify_order)
    show: BoolProperty(name="show", default=False,
                            update=Callback_modify_show)
