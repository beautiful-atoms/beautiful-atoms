import bpy
from bpy.types import Panel
from bpy.props import (StringProperty,
                       BoolProperty,
                       FloatProperty,
                       EnumProperty,
                       )
from batoms.utils.butils import get_selected_batoms, get_selected_vertices
from batoms import Batoms


class Batoms_PT_prepare(Panel):
    bl_label = "Batoms"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def draw(self, context):
        obj = context.object
        layout = self.layout
        layout.label(text="Active: " + obj.batoms.label)
        layout.operator("batoms.import")
        layout.operator("batoms.export")

        bapanel = context.scene.bapanel

        layout.label(text="Model style")
        col = layout.column()
        col.prop(bapanel, "model_style", expand=True)
        layout.label(text="Radius style")
        layout.prop(bapanel, "radius_style", expand=True)

        layout.prop(bapanel, "show", expand=True)
        layout.prop(bapanel, "wrap", expand=True)
        layout.prop(bapanel, "scale")

        layout.operator("batoms.replace")

        


class BatomsProperties(bpy.types.PropertyGroup):
    def Callback_model_style(self, context):
        bapanel = bpy.context.scene.bapanel
        model_style = list(bapanel.model_style)[0]
        modify_batoms_attr(context, 'model_style', model_style)

    def Callback_radius_style(self, context):
        bapanel = bpy.context.scene.bapanel
        radius_style = list(bapanel.radius_style)[0]
        modify_batoms_attr(context, 'radius_style', radius_style)

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
        min = 0.0, soft_max = 2.0,
        description="scale", update=Callback_modify_scale)


def modify_batoms_attr(context, key, value):
    """
    """
    batoms = Batoms(label=context.object.batoms.label)
    setattr(batoms, key, value)
    # batoms.obj.select_set(True)
    bpy.context.view_layer.objects.active = batoms.obj
