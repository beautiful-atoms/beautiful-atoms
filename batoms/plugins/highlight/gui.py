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



class HighlightProperties(bpy.types.PropertyGroup):
    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_attr("show", get_active_bpy_data('Bhighlight')),
                       set=set_attr("show", set_module_attr('highlight'))
                       )
    atomRadius: FloatProperty(name="atomRadius",
                default=0.5,
                description="average radius of the atoms",
                get=get_attr("atomRadius", get_active_bpy_data('Bhighlight')),
                set=set_attr("atomRadius", set_module_attr('highlight'))
                )
    minCave: FloatProperty(name="minCave",
                default=3,
                description="minmum radius of the cave",
                get=get_attr("minCave", get_active_bpy_data('Bhighlight')),
                set=set_attr("minCave", set_module_attr('highlight'))
                )
    resolution: FloatProperty(name="resolution",
                default=1,
                description="resolution",
                get=get_attr("resolution", get_active_bpy_data('Bhighlight')),
                set=set_attr("resolution", set_module_attr('highlight'))
                )



class VIEW3D_PT_Batoms_highlight(Panel):
    bl_label = "Highlight"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Plugins"
    bl_idname = "VIEW3D_PT_Batoms_highlight"
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
        iso = context.scene.Bhighlight

        layout.separator()


class BATOMS_MT_highlight_context_menu(Menu):
    bl_label = "Highlight Specials"
    bl_idname = "BATOMS_MT_highlight_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("batoms.highlight_add",
                             icon='ADD', text="Add highlight")
        layout.separator()
        op = layout.operator("batoms.highlight_remove",
                             icon='X', text="Delete All highlight")
        op.all = True


class BATOMS_UL_highlight(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        highlight = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.5, align=False)
            split.prop(highlight, "name", text="",
                       emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(highlight, "select", text="")
            row.prop(highlight, "scale", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_highlight(Panel):
    bl_label = "Item settings"
    bl_category = "Surface"
    bl_idname = "BATOMS_PT_highlight"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_parent_id = 'VIEW3D_PT_Batoms_highlight'
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
        ba = bpy.data.collections[ob.batoms.label].Bhighlight
        if len(ba.settings) > 0:
            kb = ba.settings[ba.ui_list_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_highlight", "", ba,
                          "settings", ba, "ui_list_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("batoms.highlight_add", icon='ADD', text="")
        op = col.operator("batoms.highlight_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_highlight_context_menu",
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

            sub = row.row()
            layout.use_property_split = True
            row = layout.row()
            col = layout.column()
            sub = col.column(align=True)
            sub.prop(kb, "select", text="select")
            sub.prop(kb, "scale", text="scale")
            sub.prop(kb, "style", text="style")
            sub.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color",  text="color")
            col.separator()
        op = layout.operator("batoms.highlight_draw",
                             icon='GREASEPENCIL', text="Draw")
