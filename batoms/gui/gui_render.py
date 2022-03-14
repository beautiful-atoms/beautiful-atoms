import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (
    IntVectorProperty,
    FloatProperty,
)
from batoms.render.render import Render


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
        row = layout.row()
        row.prop(repanel, "viewport")
        layout.prop(repanel, "distance")
        layout.operator("batoms.add_render")


class RenderProperties(bpy.types.PropertyGroup):
    """
    """
    def Callback_modify_viewport(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(context, 'viewport', repanel.viewport)

    def Callback_modify_distance(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(context, 'distance', repanel.distance)

    viewport: IntVectorProperty(
        name="Viewport", size=3, default=(0, 0, 1),
        description="Miller viewport for the render", update=Callback_modify_viewport)
    distance: FloatProperty(
        name="Distance", default=3,
        description="distance from origin", update=Callback_modify_distance)


def modify_render_attr(context, key, value):
    from batoms.batoms import Batoms
    render = Render()
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        render.batoms = batoms
        setattr(render, key, value)


class AddButton(Operator):
    bl_idname = "batoms.add_render"
    bl_label = "Add"
    bl_description = "Add distance, angle and dihedra angle"

    def execute(self, context):
        repanel = context.scene.repanel
        render = Render(viewport=repanel.viewport,
                        distance=repanel.distance)
        return {'FINISHED'}
