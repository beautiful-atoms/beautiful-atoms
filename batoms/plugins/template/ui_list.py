"""
"""

import bpy
from bpy.types import Menu, Panel, UIList


class BATOMS_MT_template_context_menu(Menu):
    bl_label = "Molecular Surface Specials"
    bl_idname = "BATOMS_MT_template_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("surface.template_add", icon='ADD',
                             text="Add Molecular Surface")
        layout.separator()
        op = layout.operator("surface.template_remove", icon='X',
                             text="Delete All Molecular Surface")
        op.all = True


class BATOMS_UL_template(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        template = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.5, align=False)
            split.prop(template, "name", text="", emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(template, "type", text="")
            row.prop(template, "select", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_template(Panel):
    bl_label = "template"
    bl_category = "Surface"
    bl_idname = "BATOMS_PT_template"
    bl_space_type = "VIEW_3D"
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
        if len(ba.btemplate) > 0:
            kb = ba.btemplate[ba.template_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_template", "", ba, "btemplate",
                          ba, "template_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("surface.template_add", icon='ADD', text="")
        op = col.operator("surface.template_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_template_context_menu", icon='DOWNARROW_HLT', text="")

        if kb:
            col.separator()

            sub = col.column(align=True)

            split = layout.split(factor=0.4)
            row = split.row()

            row = split.row()
            row.alignment = 'RIGHT'

            sub = row.row(align=True)
            sub.label()  # XXX, for alignment only

            sub = row.row()
            layout.use_property_split = True
            row = layout.row()
            col = layout.column()
            sub = col.column(align=True)
            sub.prop(kb, "type", text="Type")
            sub.prop(kb, "scale", text="Scale")
            # sub.prop(kb, "resolution", text="Resolution")
            sub.prop(kb, "select", text="Select")
            col.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color",  text="color")
            col.separator()
            op = layout.operator(
                "surface.template_draw", icon='GREASEPENCIL', text="Draw")
