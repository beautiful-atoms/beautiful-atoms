"""
"""

import bpy
from bpy.types import Menu, Panel, UIList


class BATOMS_MT_molecular_surface_context_menu(Menu):
    bl_label = "Molecular Surface Specials"
    bl_idname = "BATOMS_MT_molecular_surface_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("surface.molecular_surface_add", icon='ADD',
                             text="Add Molecular Surface")
        layout.separator()
        op = layout.operator("surface.molecular_surface_remove", icon='X',
                             text="Delete All Molecular Surface")
        op.all = True


class BATOMS_UL_molecular_surface(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        ms = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.5, align=False)
            split.prop(ms, "name", text="", emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(ms, "type", text="")
            row.prop(ms, "select", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_molecular_surface(Panel):
    bl_label = "Molecular Surface"
    bl_category = "Plugins"
    bl_idname = "BATOMS_PT_molecular_surface"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}

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
        ba = bpy.data.collections[ob.batoms.label].Bmolecularsurface
        if len(ba.settings) > 0:
            kb = ba.settings[ba.ui_list_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_molecular_surface", "", ba, "settings",
                          ba, "ui_list_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("surface.molecular_surface_add", icon='ADD', text="")
        op = col.operator("surface.molecular_surface_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_molecular_surface_context_menu", icon='DOWNARROW_HLT', text="")

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
            sub.prop(kb, "probe", text="Probe")
            sub.prop(kb, "resolution", text="Resolution")
            sub.prop(kb, "select", text="Select")
            col.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color_by",  text="Color by")
            if kb.color_by == 'None':
                col.prop(kb, "color",  text="color")
            else:
                col.prop(kb, "color1",  text="color1")
                col.prop(kb, "color2",  text="color2")
            col.separator()
            op = layout.operator(
                "surface.molecular_surface_draw", icon='GREASEPENCIL', text="Draw")
