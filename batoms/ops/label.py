"""
"""
import bpy
from bpy_extras import view3d_utils
from bpy.props import (StringProperty,
                       IntProperty,
                       IntVectorProperty,
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms
import blf
import gpu
from gpu_extras.batch import batch_for_shader
import numpy as np

def draw_callback_text(self, context):
    # get the context arguments
    if len(self._positions) == 0:
        return
    scene = context.scene
    region = context.region
    rv3d = context.region_data
    # coord = event.mouse_region_x, event.mouse_region_y
    font_id = 0  # XXX, need to find out how best to get this.
    coord = np.mean(self._positions, axis = 0)
    coord_2d = view3d_utils.location_3d_to_region_2d(region, rv3d, coord, default=None)
    # draw some text
    blf.position(font_id, coord_2d[0], coord_2d[1] + 10, 0)
    blf.size(font_id, 20, 72)
    blf.draw(font_id, self.results)


class BatomsLabel(OperatorBatoms):
    """ """
    bl_idname = "batoms.label"
    bl_label = "Label"

    
    label: StringProperty(
        name="Label", default='Index',
        description="Label")
    
    select: StringProperty(
        name="Label", default='Index',
        description="Label")

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.show_label(self.label)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}
