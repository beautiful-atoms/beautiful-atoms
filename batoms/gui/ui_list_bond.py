"""
startup/bl_ui/properties_data_mesh.py

"""

import bpy
from bpy.types import Menu, Panel, UIList
from collections import defaultdict
from batoms import Batoms

class BATOMS_MT_bond_pair_context_menu(Menu):
    bl_label = "Bond Pair Specials"
    bl_idname = "BATOMS_MT_bond_pair_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("bond.bond_pair_add", icon='ADD', text="Add Bond Pair")
        layout.separator()
        # layout.operator("object.shape_key_mirror", icon='ARROW_LEFTRIGHT').use_topology = False
        # layout.operator("object.shape_key_mirror", text="Mirror Shape Key (Topology)").use_topology = True
        # layout.separator()
        # layout.operator("object.join_shapes")
        # layout.operator("object.shape_key_transfer")
        # layout.separator()
        op = layout.operator("bond.bond_pair_remove", icon='X', text="Delete All Bond Pair")
        op.all = True

        # layout.separator()
        # layout.operator("object.shape_key_move", icon='TRIA_UP_BAR', text="Move to Top").type = 'TOP'
        # layout.operator("object.shape_key_move", icon='TRIA_DOWN_BAR', text="Move to Bottom").type = 'BOTTOM'


class BATOMS_UL_bond_pairs(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        bond_pair = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.66, align=False)
            split.prop(bond_pair, "name", text="", emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(bond_pair, "max", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_bond_pairs(Panel):
    bl_label       = "Bond Pair"
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Pairs"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
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
        ba = bpy.data.collections[ob.batoms.label].batoms
        if len(ba.bbond) >0:
            kb = ba.bbond[ba.bond_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_bond_pairs", "", ba, "bbond", ba, "bond_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("bond.bond_pair_add", icon='ADD', text="")
        op = col.operator("bond.bond_pair_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_bond_pair_context_menu", icon='DOWNARROW_HLT', text="")

        if kb:
            col.separator()

            sub = col.column(align=True)
            # sub.operator("object.shape_key_move", icon='TRIA_UP', text="").type = 'UP'
            # sub.operator("object.shape_key_move", icon='TRIA_DOWN', text="").type = 'DOWN'

            split = layout.split(factor=0.4)
            row = split.row()

            row = split.row()
            row.alignment = 'RIGHT'

            sub = row.row(align=True)
            sub.label()  # XXX, for alignment only
            subsub = sub.row(align=True)
            # subsub.active = enable_pin
            # subsub.prop(ob, "show_only_shape_key", text="")
            # sub.prop(ob, "use_shape_key_edit_mode", text="")

            sub = row.row()
            # if key.use_relative:
                # sub.operator("object.shape_key_clear", icon='X', text="")
            # else:
                # sub.operator("object.shape_key_retime", icon='RECOVER_LAST', text="")

            layout.use_property_split = True
            row = layout.row()
            # row.active = enable_edit_value
            # row.prop(kb, "select")
            col = layout.column()
            # sub.active = enable_edit_value
            sub = col.column(align=True)
            sub.prop(kb, "min", text="Min")
            sub.prop(kb, "max", text="Max")
            col.prop(kb, "search", text="Search")
            col.prop(kb, "style", text="Style")
            col.prop(kb, "order", text="Order")
            col.prop(kb, "polyhedra",  text="Polyhedra")
            col.prop(kb, "color1",  text="color1")
            col.prop(kb, "color2",  text="color2")
        #
        
        
