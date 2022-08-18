"""
Manipulate atoms and molecules with mouse interactively
based on force field methods.
"""

import bpy
import numpy as np
from batoms.utils.butils import (get_selected_batoms, get_selected_vertices_bmesh)
from batoms import Batoms
from .modal_rigid_body import mouse2positions
from bpy.props import (StringProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       )

import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def translate(batoms, selected_vertices, displacement, 
        speed=1, steps=5):
    """
    """
    if len(selected_vertices) == 0:
        return
    displacement = displacement[:3]*speed/100
    # print("Label: ", batoms.label)
    # print("selected_vertices: ", selected_vertices)
    # Move selected atom with mouse
    positions = batoms.positions
    positions[selected_vertices] += displacement
    batoms.positions = positions
    # optimize
    batoms.optmize(steps=steps)


class OB_Force_Field_Operator(bpy.types.Operator):
    """Force field"""
    bl_idname = "batoms.ob_force_field_operator"
    bl_label = "Open Babel Force field translate Operator"

    mouse_position: FloatVectorProperty(size=4, default=(0, 0, 0, 0))
    nframe: IntProperty(default=0)
    previous: StringProperty()

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'BATOMS'
        else:
            return False

    def modal(self, context, event):
        """
        """
        if event.type in {'ESC'}:
            return {'CANCELLED'}
        elif event.type == 'MIDDLEMOUSE':
            self.previous = 'MIDDLEMOUSE'
            return {'PASS_THROUGH'}
        elif event.type == 'MOUSEMOVE':
            mouse_position = np.array([event.mouse_x, event.mouse_y, 0, 0])
            if self.previous == 'MIDDLEMOUSE':
                self.mouse_position = mouse_position
                self.previous = 'MOUSEMOVE'
                return {'RUNNING_MODAL'}
            # read the batoms
            batoms = Batoms(context.object.batoms.label)
            # read the selected vertices
            selected_vertices = get_selected_vertices_bmesh(batoms.obj)
            # print("selected_vertices:", selected_vertices)
            if len(selected_vertices) == 0:
                self.mouse_position = mouse_position
                return {'RUNNING_MODAL'}
            delta = mouse_position - self.mouse_position
            displacement = mouse2positions(
                delta, self.viewports_3D.spaces.active.region_3d.    view_matrix)
            self.mouse_position = mouse_position
            obffpanel = context.scene.obffpanel
            speed = obffpanel.speed
            steps = obffpanel.steps
            translate(batoms, selected_vertices, displacement,
                      speed, steps)
            self.nframe += 1
            self.previous = 'MOUSEMOVE'
            return {'RUNNING_MODAL'}
        else:
            return {'PASS_THROUGH'}

    def invoke(self, context, event):
        """
        """
        if context.object:
            for area in bpy.context.screen.areas:
                if area.type == 'VIEW_3D':
                    self.viewports_3D = area
            self.mouse_position = np.array(
                [event.mouse_x, event.mouse_y, 0, 0])
            self.nframe = 0
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}


class BATOMS_PT_OB_Force_Field_Modal(bpy.types.Panel):
    bl_idname = "BATOMS_PT_OB_Force_Field_Modal"
    bl_label = "Open Babel Force field"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Real"

    def draw(self, context):
        obffpanel = context.scene.obffpanel
        layout = self.layout
        box = layout.box()
        col = box.column()
        col.prop(obffpanel, "speed", expand=True)
        col.prop(obffpanel, "steps", expand=True)
        layout = self.layout
        layout.operator("batoms.ob_force_field_operator",
                        text='Force field', icon='FACESEL')


class OBForceFieldProperties(bpy.types.PropertyGroup):

    speed: FloatProperty(
        name="Speed", default=1,
        description="Speed to move atoms")
    steps: IntProperty(
        name="Steps", default=5,
        description="Optimization steps")
