"""
"""
import bpy
import bmesh
from batoms import Batoms
import blf
import gpu
from gpu_extras.batch import batch_for_shader


def draw_callback_px(self, context):
    font_id = 0  # XXX, need to find out how best to get this.
    # draw some text
    blf.position(font_id, 100, 300, 0)
    blf.size(font_id, 2, 7)
    blf.draw(font_id, self.results)

    # 50% alpha, 2 pixel width line
    shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
    # gpu.state.blend_set('ALPHA')
    gpu.state.line_width_set(2.0)
    # print(self.positions)
    batch = batch_for_shader(shader, 'LINES', {"pos": self.positions})
    shader.bind()
    shader.uniform_float("color", (0.0, 0.0, 0.0, 1))
    batch.draw(shader)

    # restore opengl defaults
    gpu.state.line_width_set(1.0)
    gpu.state.blend_set('NONE')


class MeasureButton(bpy.types.Operator):
    """Records the selection order while running and
    when finished with ESC """
    bl_idname = "batoms.measure"
    bl_label = "Measurement"

    results = ""

    @classmethod
    def poll(cls, context):
        return context.mode in {'EDIT_MESH'}

    def modal(self, context, event):
        context.area.tag_redraw()
        if event.type in {'RIGHTMOUSE'}:
            obj = context.object
            data = obj.data
            bm = bmesh.from_edit_mesh(data)
            v = [s.index for s in bm.select_history if isinstance(
                s, bmesh.types.BMVert)]
            batoms = Batoms(label=obj.batoms.label)
            results = []
            for i in v:
                self.positions.append(tuple(batoms[i].position))
            if len(v) == 1:
                measurement_type = 'Position: '
                results = batoms[v[0]].position
            elif len(v) == 2:
                results = batoms.get_distances(v[0],
                                               [v[1]])
                measurement_type = 'Bond length: '
            elif len(v) == 3:
                results = batoms.get_angle(v[0], v[1], v[2])
                measurement_type = 'Angle: '
            elif len(v) == 4:
                results = batoms.get_dihedrals(v[0], v[1], v[2], v[3])
                measurement_type = 'Dihedral: '
            else:
                measurement_type = "More than four atoms are selected. Unsupported."
            # update measurement value
            results.shape = (-1,)
            results = [str(round(float(i), 2)) for i in results]
            results = measurement_type + ' '.join(results)
            # self.report({'INFO'}, results)
            bpy.ops.object.mode_set(mode='EDIT')
            self.results = results
            return {'RUNNING_MODAL'}
        elif event.type in {'ESC'}:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'CANCELLED'}
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        if context.area.type == 'VIEW_3D':
            # the arguments we pass the the callback
            args = (self, context)
            # Add the region OpenGL drawing callback
            # draw in view space with 'POST_VIEW' and 'PRE_VIEW'
            self._handle = bpy.types.SpaceView3D.draw_handler_add(
                draw_callback_px, args, 'WINDOW', 'POST_VIEW')
            self.results = ""
            self.positions = []
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "View3D not found, cannot run operator")
            return {'CANCELLED'}
