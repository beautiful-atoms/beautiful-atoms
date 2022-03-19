import bpy
from bpy_extras import view3d_utils
import blf


class DrawText:
    def __init__(self, context, label):
        self.label = label

    def draw_callback_text(self, positions, texts):
        # get the context arguments
        context = bpy.context
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        font_id = 0
        n = len(positions)
        for i in range(n):
            coord = positions[i]
            coord_2d = view3d_utils.location_3d_to_region_2d(
                region, rv3d, coord, default=None)
            # draw some text
            blf.position(font_id, coord_2d[0], coord_2d[1], 0)
            blf.size(font_id, 20, 72)
            blf.draw(font_id, str(texts[i]))

    def add_handle(self, positions, texts):
        self._handle_label = bpy.types.SpaceView3D.draw_handler_add(
            self.draw_callback_text, (positions, texts), 'WINDOW', 'POST_PIXEL')

    def remove_handle(self):
        if hasattr(self, "_handle_label"):
            bpy.types.SpaceView3D.draw_handler_remove(
                self._handle_label, 'WINDOW')
            delattr(self, "_handle_label")
