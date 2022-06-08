import bpy
import bmesh
from bpy.props import (
    StringProperty,
    BoolProperty,
    FloatVectorProperty,
    FloatProperty,
    IntProperty,
    IntVectorProperty
)
from batoms import Batoms
from batoms.ops.base import OperatorBatoms, OperatorBatomsEdit


class RenderAdd(OperatorBatoms):
    bl_idname = "batoms.render_add"
    bl_label = "Attach Render"
    bl_description = "Add a render to a Batoms object."

    viewport: FloatVectorProperty(
        name="viewport", default=(0, 0, 1), size=3,
        soft_min=-1, soft_max=1,
        description="Viewport")

    resolution: FloatVectorProperty(
        name="resolution", default=(1000, 1000), size=2,
        soft_min=-1, soft_max=1,
        description="resolution")

    def execute(self, context):
        from batoms.render.render import Render
        batoms = Batoms(label=context.object.batoms.label)
        render = Render(viewport=self.viewport,
                        resolution=self.resolution,
                        )
        batoms.render = render
        batoms.render.init()
        bpy.context.space_data.region_3d.view_perspective = 'CAMERA'
        return {'FINISHED'}
