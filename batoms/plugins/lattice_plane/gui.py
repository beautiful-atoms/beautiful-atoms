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



class LatticePlaneProperties(bpy.types.PropertyGroup):
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=model_style_items,
        get=get_enum_attr("model_style", get_active_bpy_data('Blatticeplane')),
        set=set_enum_attr("model_style", set_module_attr('lattice_plane')),
        default=0,
    )

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_attr("show", get_active_bpy_data('Blatticeplane')),
                       set=set_attr("show", set_module_attr('lattice_plane'))
                       )



class VIEW3D_PT_Batoms_lattice_plane(Panel):
    bl_label = "Lattice Plane"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Plugins"
    bl_idname = "VIEW3D_PT_Batoms_lattice_plane"
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
        iso = context.scene.Blatticeplane

        layout.label(text="Model style")
        layout.prop(iso, "model_style", expand=True)
        layout.prop(iso, "show", expand=True)
        layout.separator()


class BATOMS_MT_lattice_plane_context_menu(Menu):
    bl_label = "Lattice Plane Specials"
    bl_idname = "BATOMS_MT_lattice_plane_context_menu"

    def draw(self, _context):
        layout = self.layout
        op = layout.operator("plane.lattice_plane_add",
                             icon='ADD', text="Add Lattice Plane")
        layout.separator()
        op = layout.operator("plane.lattice_plane_remove",
                             icon='X', text="Delete All Lattice Plane")
        op.all = True


class BATOMS_UL_lattice_plane(UIList):
    def draw_item(self, _context, layout, _data, item, icon, active_data, _active_propname, index):
        lattice_plane = item
        custom_icon = 'OBJECT_DATAMODE'
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.66, align=False)
            split.prop(lattice_plane, "name", text="",
                       emboss=False, icon=custom_icon)
            row = split.row(align=True)
            row.emboss = 'NONE_OR_STATUS'
            row.prop(lattice_plane, "distance", text="")
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon=custom_icon)


class BATOMS_PT_lattice_plane(Panel):
    bl_label = "Item settings"
    bl_category = "Plugins"
    bl_idname = "BATOMS_PT_Lattice_Plane"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_parent_id = 'VIEW3D_PT_Batoms_lattice_plane'
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
        ba = bpy.data.collections[ob.batoms.label].Blatticeplane
        if len(ba.settings) > 0:
            kb = ba.settings[ba.ui_list_index]
        else:
            kb = None

        row = layout.row()

        rows = 3
        if kb:
            rows = 5
        row.template_list("BATOMS_UL_lattice_plane", "", ba,
                          "settings", ba, "ui_list_index", rows=rows)

        col = row.column(align=True)
        op = col.operator("plane.lattice_plane_add", icon='ADD', text="")
        op = col.operator("plane.lattice_plane_remove", icon='REMOVE', text="")
        if kb is not None:
            op.name = kb.name
        col.separator()

        col.menu("BATOMS_MT_lattice_plane_context_menu",
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
            row.prop(kb, "indices", text="Miller index")
            col = layout.column()
            sub = col.column(align=True)
            sub.prop(kb, "distance", text="Distance")
            sub.prop(kb, "scale", text="Scale")
            col.prop(kb, "show_edge",  text="Show edge")
            col.prop(kb, "material_style", text="material_style")
            col.prop(kb, "color_by",  text="Color by")
            if kb.color_by == 'None':
                col.prop(kb, "color",  text="color")
            else:
                col.prop(kb, "color1",  text="color1")
                col.prop(kb, "color2",  text="color2")
            op = layout.operator("plane.lattice_plane_draw",
                                 icon='GREASEPENCIL', text="Draw")
