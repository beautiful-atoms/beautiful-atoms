import bpy
from bpy.types import Panel
from bpy.props import (
    BoolProperty,
    FloatProperty,
    EnumProperty,
)
from batoms import Batoms


class Batoms_PT_prepare(Panel):
    bl_label = "Batoms"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    def draw(self, context):
        name = 'None'
        if context.object:
            if context.object.batoms.type != 'OTHER':
                name = context.object.batoms.label
        layout = self.layout
        layout.label(text="Active: " + name)
        layout.operator("batoms.import")
        layout.operator("batoms.export")

        bapanel = context.scene.bapanel

        layout.label(text="Model style")
        col = layout.column()
        col.prop(bapanel, "model_style", expand=True)
        layout.label(text="Add label")
        layout.prop(bapanel, "show_label", expand=True)
        layout.label(text="Radius style")
        layout.prop(bapanel, "radius_style", expand=True)
        layout.label(text="Color style")
        layout.prop(bapanel, "color_style", expand=True)
        layout.label(text="Polyhedra style")
        layout.prop(bapanel, "polyhedra_style", expand=True)

        layout.prop(bapanel, "show", expand=True)
        layout.prop(bapanel, "wrap", expand=True)
        layout.prop(bapanel, "scale")

        layout.operator("batoms.replace")


class BatomsProperties(bpy.types.PropertyGroup):
    def Callback_model_style(self, context):
        bapanel = bpy.context.scene.bapanel
        model_style = list(bapanel.model_style)[0]
        modify_batoms_attr(context, 'model_style', model_style)
    
    def Callback_show_label(self, context):
        bapanel = bpy.context.scene.bapanel
        show_label = list(bapanel.show_label)[0]
        if show_label == 'none':
            show_label = None
        modify_batoms_attr(context, 'show_label', show_label)

    def Callback_radius_style(self, context):
        bapanel = bpy.context.scene.bapanel
        radius_style = list(bapanel.radius_style)[0]
        modify_batoms_attr(context, 'radius_style', radius_style)

    def Callback_color_style(self, context):
        bapanel = bpy.context.scene.bapanel
        color_style = list(bapanel.color_style)[0]
        modify_batoms_attr(context, 'color_style', color_style)

    def Callback_polyhedra_style(self, context):
        bapanel = bpy.context.scene.bapanel
        polyhedra_style = list(bapanel.polyhedra_style)[0]
        modify_batoms_attr(context, 'polyhedra_style', polyhedra_style)

    def Callback_modify_show(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(context, 'show', bapanel.show)

    def Callback_modify_wrap(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(context, 'wrap', bapanel.wrap)

    def Callback_modify_scale(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(context, 'scale', bapanel.scale)

    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0', "Space-filling", "Use ball and stick"),
               ('1', "Ball-and-stick", "Use ball"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Stick", "Use stick")),
        default={'0'},
        update=Callback_model_style,
        options={'ENUM_FLAG'},
    )

    show_label: EnumProperty(
        name="show_label",
        description="Structural models",
        items=(('none', "None", "None"),
               ('index', "Index", "Index"),
               ('species', "Species", "Species"),
               ('charge', "Charge", "charge"),
               ),
        default={'none'},
        update=Callback_show_label,
        options={'ENUM_FLAG'},
    )

    radius_style: EnumProperty(
        name="radius_style",
        description="Structural models",
        items=(('0', "Covalent", "covalent"),
               ('1', "VDW", "van der Waals"),
               ('2', "Ionic", "ionic")),
        default={'0'},
        update=Callback_radius_style,
        options={'ENUM_FLAG'},
    )

    color_style: EnumProperty(
        name="color_style",
        description="Color",
        items=(('0', "JMOL", "JMOL"),
               ('1', "VESTA", "VESTA"),
               ('2', "CPK", "CPK")),
        default={'0'},
        update=Callback_color_style,
        options={'ENUM_FLAG'},
    )

    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(('0', "0", "atoms, bonds and polyhedra"),
               ('1', "1", "atoms, polyhedra"),
               ('2', "2", "central atoms, polyhedra"),
               ('3', "3", "polyhedra")),
        default='0',
        update=Callback_polyhedra_style,
    )

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       update=Callback_modify_show)

    wrap: BoolProperty(name="wrap",
                       default=False,
                       description="wrap all atoms into cell",
                       update=Callback_modify_wrap)

    scale: FloatProperty(
        name="scale", default=1.0,
        min=0.0, soft_max=2.0,
        description="scale", update=Callback_modify_scale)


def modify_batoms_attr(context, key, value):
    """
    """
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        setattr(batoms, key, value)
        # batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
