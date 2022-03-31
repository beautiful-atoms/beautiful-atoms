import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (
    FloatVectorProperty,
    IntVectorProperty,
    FloatProperty,
)
from batoms.render.render import Render
from batoms import Batoms

class Render_PT_prepare(Panel):
    bl_label = "Render"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "RENDER_PT_Tools"

    def draw(self, context):
        layout = self.layout
        repanel = context.scene.repanel

        layout = self.layout
        col = layout.column()
        # row = col.row(align=True)
        col.prop(repanel, "viewport")
        col.prop(repanel, "light_direction")
        layout.prop(repanel, "resolution")
        layout.prop(repanel, "scale")
        layout.prop(repanel, "distance")


class RenderProperties(bpy.types.PropertyGroup):
    """
    """
    def Callback_modify_viewport(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(context, 'viewport', repanel.viewport)
    
    def Callback_modify_light_direction(self, context):
        repanel = bpy.context.scene.repanel
        modify_light_attr(context, 'direction', repanel.light_direction)

    def Callback_modify_resolution(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(context, 'resolution', repanel.resolution)

    def Callback_modify_scale(self, context):
        modify_camera_attr(context, 'scale', context.scene.repanel.scale)

    def Callback_modify_distance(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(context, 'distance', repanel.distance)

    viewport: FloatVectorProperty(
        name="Viewport", size=3, default=(0, 0, 1),
        soft_min = -1, soft_max = 1,
        subtype = "XYZ",
        description="Miller viewport for the render", update=Callback_modify_viewport)
    
    light_direction: FloatVectorProperty(
        name="Light direction", size=3, default=(0, 0, 1),
        soft_min = -1, soft_max = 1,
        subtype = "XYZ",
        description="Light direction for the render", update=Callback_modify_light_direction)
    

    resolution: IntVectorProperty(
        name="resolution", size=2, default=(1000, 1000),
        soft_min = 100, soft_max = 2000,
        description="Miller resolution for the render", update=Callback_modify_resolution)
    
    scale: FloatProperty(
        name="scale", default=10,
        soft_min = 1, soft_max = 100,
        description="scale", update=Callback_modify_scale)


    distance: FloatProperty(
        name="Distance", default=3,
        description="distance from origin", update=Callback_modify_distance)


def modify_render_attr(context, key, value):
    from batoms.batoms import Batoms
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        setattr(batoms.render, key, value)
        batoms.render.init()
        context.space_data.region_3d.view_perspective = 'CAMERA'


def modify_light_attr(context, key, value):
    from batoms.batoms import Batoms
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        setattr(batoms.render.lights['Default'], key, value)
        batoms.render.init()
        context.space_data.region_3d.view_perspective = 'CAMERA'

def modify_camera_attr(context, key, value):
    from batoms.batoms import Batoms
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        if key.upper() == 'SCALE':
            batoms.render.camera.set_ortho_scale(value)
        else:
            setattr(batoms.render.camera, key, value)
            batoms.render.init()
        context.space_data.region_3d.view_perspective = 'CAMERA'
