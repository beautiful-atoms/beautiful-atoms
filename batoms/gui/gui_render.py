import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (BoolProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       )
from batoms.utils.butils import get_selected_batoms
from batoms.batoms import Batoms
from batoms.render.render import Render

# The panel.
class Render_PT_prepare(Panel):
    bl_label       = "Render"
    bl_space_type  = "VIEW_3D"
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
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    def Callback_modify_viewport(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(self.selected_batoms, 'viewport', repanel.viewport)
    def Callback_modify_distance(self, context):
        repanel = bpy.context.scene.repanel
        modify_render_attr(self.selected_batoms, 'distance', repanel.distance)
    
    
    viewport: IntVectorProperty(
        name="Viewport", size = 3, default=(0, 0, 1),
        description = "Miller viewport for the render", update = Callback_modify_viewport)
    distance: FloatProperty(
        name="Distance", default=3,
        description = "distance from origin", update = Callback_modify_distance)
    
    

def modify_render_attr(selected_batoms, key, value):
    render = Render()
    if len(selected_batoms) == 1:
        render.batoms = Batoms(label = selected_batoms[0])
        setattr(render, key, value)

def add_render(viewport, distance):
    """
    """
    render = Render(viewport = viewport, distance=distance)

class AddButton(Operator):
    bl_idname = "batoms.add_render"
    bl_label = "Add"
    bl_description = "Add distance, angle and dihedra angle"

    def execute(self, context):
        repanel = context.scene.repanel
        add_render(repanel.viewport, 
                    repanel.distance, 
                    )
        return {'FINISHED'}