"""
"""

import bpy
from bpy.types import Menu, Panel, UIList
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       )

from batoms.internal_data.bpy_data import get_volumetric_data


class VolumetricDataProperties(bpy.types.PropertyGroup):
    # volumetric data
    name: StringProperty(
        name="name",
        description="Name of volumetric data",
        default="create",
        )
    select_data1: EnumProperty(
        name="select_data1",
        description="Select volumetric data",
        items=get_volumetric_data,
        default=0
        )
    select_data2: EnumProperty(
        name="select_data2",
        description="Select volumetric data",
        items=get_volumetric_data,
        default=0
        )
    operator: EnumProperty(
        name="operator",
        items=[('Add', 'Add', '', 0),
                ('Minus', 'Minus', '', 1)
              ],
        description="Volumetric data to be used.",
        default=0,
        )
class VIEW3D_PT_Batoms_volumetric_data(Panel):
    bl_label = "Volumetric Data"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_volumetric_data"
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
        batoms = context.scene.batoms
        op = layout.operator("batoms.volumetric_data_create",
                             icon='GREASEPENCIL', text="Create from others")
        op.name = batoms.volumetric_data.name
        op.select_data1 = batoms.volumetric_data.select_data1
        op.select_data2 = batoms.volumetric_data.select_data2
        op.operator = batoms.volumetric_data.operator
        layout.separator()

class BATOMS_MT_volumetric_data_context_menu(Menu):
    bl_label = "Volumetric Data Specials"
    bl_idname = "BATOMS_MT_volumetric_data_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("batoms.volumetric_data_add",
                             icon='ADD', text="Add Data")
        layout.separator()
        op = layout.operator("batoms.volumetric_data_remove",
                             icon='X', text="Delete All Data")
        op.all = True


class BATOMS_UL_volumetric_data(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        volumetric_data = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.66, align=False)
            split.prop(volumetric_data, "name", text="",
                       emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            # row.prop(volumetric_data, "distance", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_volumetric_data(Panel):
    bl_label = "Volumetric_data"
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_volumetric_data"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_parent_id = 'VIEW3D_PT_Batoms_volumetric_data'
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
        if len(ba.settings_volumetric_data) > 0:
            kb = ba.settings_volumetric_data[ba.ui_list_index_volumetric_data]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5

        row.template_list("BATOMS_UL_volumetric_data", "", ba,
                          "settings_volumetric_data", ba, "ui_list_index_volumetric_data", rows=rows)

        col = row.column(align=True)
        op = col.operator("batoms.volumetric_data_add", icon='ADD', text="")
        op = col.operator("batoms.volumetric_data_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_volumetric_data_context_menu",
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
            # sub.prop(kb, "radius_style", text="Radius_style")
            col.prop(kb, "shape",  text="Shape")
