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



class TemplateProperties(bpy.types.PropertyGroup):
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=model_style_items,
        get=get_enum_attr("model_style", get_active_bpy_data('Btemplate')),
        set=set_enum_attr("model_style", set_module_attr('template')),
        default=0,
    )

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_attr("show", get_active_bpy_data('Btemplate')),
                       set=set_attr("show", set_module_attr('template'))
                       )



class VIEW3D_PT_Batoms_template(Panel):
    bl_label = "template"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Plugins"
    bl_idname = "VIEW3D_PT_Batoms_template"
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
        iso = context.scene.Btemplate

        layout.label(text="Model style")
        layout.prop(iso, "model_style", expand=True)
        layout.prop(iso, "show", expand=True)
        layout.separator()

class BATOMS_MT_template_context_menu(Menu):
    bl_label = "Molecular Surface Specials"
    bl_idname = "BATOMS_MT_template_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("template.template_add", icon='ADD',
                             text="Add Molecular Surface")
        layout.separator()
        op = layout.operator("template.template_remove", icon='X',
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
            row.prop(template, "prop1", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_template(Panel):
    bl_label = "Item settings"
    bl_category = "Surface"
    bl_idname = "BATOMS_PT_template"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_parent_id = 'VIEW3D_PT_Batoms_template'
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
        ba = bpy.data.collections[ob.batoms.label].Btemplate
        if len(ba.settings) > 0:
            kb = ba.settings[ba.ui_list_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_template", "", ba, "settings",
                          ba, "ui_list_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("template.template_add", icon='ADD', text="")
        op = col.operator("template.template_remove", icon='REMOVE', text="")
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
            sub.prop(kb, "prop1", text="Prop1")
            col.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color",  text="color")
            col.separator()
            op = layout.operator(
                "template.template_draw", icon='GREASEPENCIL', text="Draw")
