"""
"""

import bpy
from bpy.types import Menu, Panel, UIList
from bpy.props import (
    BoolProperty,
    FloatProperty,
    EnumProperty,
    StringProperty,
)



from batoms import Batoms
from batoms.gui.utils import (get_active_bpy_data,
        get_attr, get_enum_attr, set_attr, set_enum_attr,
        get_active_module, set_module_attr
        )


model_style_items = [("Surface", "Surface", "", 0),
                     ("Dot", "Dot", "", 1),
                     ("Wireframe", "Wireframe", "", 2),
                     ]



class BondProperties(bpy.types.PropertyGroup):

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_attr("show", get_active_bpy_data('Bbond')),
                       set=set_attr("show", set_module_attr('bond'))
                       )
    atomRadius: FloatProperty(name="atomRadius",
                default=0.5,
                description="average radius of the atoms",
                get=get_attr("atomRadius", get_active_bpy_data('Bbond')),
                set=set_attr("atomRadius", set_module_attr('bond'))
                )
    minCave: FloatProperty(name="minCave",
                default=3,
                description="minmum radius of the cave",
                get=get_attr("minCave", get_active_bpy_data('Bbond')),
                set=set_attr("minCave", set_module_attr('bond'))
                )
    resolution: FloatProperty(name="resolution",
                default=1,
                description="resolution",
                get=get_attr("resolution", get_active_bpy_data('Bbond')),
                set=set_attr("resolution", set_module_attr('bond'))
                )



class VIEW3D_PT_Batoms_bond(Panel):
    bl_label = "Bond"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_bond"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def draw(self, context):
        name = 'None'
        if context.object:
            if context.object.batoms.type != 'OTHER':
                name = context.object.batoms.label
        layout = self.layout
        # layout.label(text="Active: " + name)
        layout.operator("bond.bond_order_auto_set", icon='MODIFIER_ON', text="Auto Set Bond Order")
        layout.separator()


class BATOMS_MT_bond_pair_context_menu(Menu):
    bl_label = "Bond Pair Specials"
    bl_idname = "BATOMS_MT_bond_pair_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("bond.bond_pair_add",
                             icon='ADD', text="Add Bond Pair")
        layout.separator()
        op = layout.operator("bond.bond_pair_remove",
                             icon='X', text="Delete All Bond Pair")
        op.all = True


class BATOMS_UL_bond_pair(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        bond_pair = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.66, align=False)
            split.prop(bond_pair, "name", text="",
                       emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(bond_pair, "max", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_bond_pair(Panel):
    bl_label = "Bond Pair"
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Pair"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_parent_id = 'VIEW3D_PT_Batoms_bond'
    # bl_options = {'DEFAULT_CLOSED'}

    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def draw(self, context):
        layout = self.layout

        ob = context.object
        ba = bpy.data.collections[ob.batoms.label].Bbond
        if len(ba.settings) > 0:
            kb = ba.settings[ba.ui_list_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_bond_pair", "", ba,
                          "settings", ba, "ui_list_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("bond.bond_pair_add", icon='ADD', text="")
        op = col.operator("bond.bond_pair_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_bond_pair_context_menu",
                 icon='DOWNARROW_HLT', text="")

        if kb:
            col.separator()

            sub = col.column(align=True)

            split = layout.split(factor=0.4)
            row = split.row()

            row = split.row()
            row.alignment = 'RIGHT'

            sub = row.row(align=True)
            sub.label()  # XXX, for alignment only
            subsub = sub.row(align=True)

            sub = row.row()

            layout.use_property_split = True
            row = layout.row()
            col = layout.column()
            sub = col.column(align=True)
            sub.prop(kb, "min", text="Min")
            sub.prop(kb, "max", text="Max")
            sub.prop(kb, "width", text="Width")
            col.prop(kb, "search", text="Search")
            col.prop(kb, "style", text="Style")
            col.prop(kb, "order", text="Order")
            col.prop(kb, "polyhedra",  text="Polyhedra")
            col.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color1",  text="color1")
            col.prop(kb, "color2",  text="color2")
            op = layout.operator("bond.draw", icon='GREASEPENCIL', text="Draw")
