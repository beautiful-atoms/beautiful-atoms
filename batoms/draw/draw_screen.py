import bpy
from bpy_extras import view3d_utils
import blf
import gpu
from gpu_extras.batch import batch_for_shader
import numpy as np


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


class DrawCrystalAxes:
    def __init__(self, context, label):
        self.label = label
        self.length = 50
        self.texts = ['a', 'b', 'c']

    def draw_callback_text(self, cell):
        # get the context arguments
        context = bpy.context
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        font_id = 0
        n = len(cell)
        origin = view3d_utils.location_3d_to_region_2d(
            region, rv3d, [0, 0, 0], default=None)
        for i in range(3):
            color = [0, 0, 0, 1]
            color[i] = 1
            coord = cell[i]
            coord_2d = np.array(view3d_utils.location_3d_to_region_2d(
                region, rv3d, coord, default=None) - origin)
            vector = coord_2d/np.linalg.norm(coord_2d)*self.length + [100, 100]
            positions = [[100, 100], vector]
            # draw some text
            blf.position(font_id, vector[0], vector[1], 0)
            blf.size(font_id, 20, 72)
            blf.draw(font_id, str(self.texts[i]))
            # 50% alpha, 2 pixel width line
            shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
            gpu.state.blend_set('ALPHA')
            gpu.state.line_width_set(4.0)
            batch = batch_for_shader(shader, 'LINE_STRIP', {"pos": positions})
            shader.bind()
            shader.uniform_float("color", tuple(color))
            batch.draw(shader)
            # restore opengl defaults
            gpu.state.line_width_set(1.0)
            gpu.state.blend_set('NONE')

    def add_handle(self, cell):
        self._handle_label = bpy.types.SpaceView3D.draw_handler_add(
            self.draw_callback_text, (cell,), 'WINDOW', 'POST_PIXEL')

    def remove_handle(self):
        if hasattr(self, "_handle_label"):
            bpy.types.SpaceView3D.draw_handler_remove(
                self._handle_label, 'WINDOW')
            delattr(self, "_handle_label")
