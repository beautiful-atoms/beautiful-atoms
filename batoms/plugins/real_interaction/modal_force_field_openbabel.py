"""
Manipulate atoms and molecules with mouse interactively
based on force field methods.


Molecule: Rigid molecules, to constrain all internal
degrees of freedom using the RATTLE-type constraints of the The FixBondLengths class to constrain all internal atomic distances.

Intermolecular: atomic radii.


"""

import bpy
import numpy as np
from batoms.utils.butils import (get_selected_batoms,
                                 get_selected_objects, get_selected_vertices_bmesh,
                                 object_mode)
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
              fmax=0.05, steps=5, frame_start=0):
    """
    """
    if len(selected_vertices) == 0:
        return
    displacement = displacement[:3]*0.01
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
            fmax = obffpanel.fmax
            steps = obffpanel.steps
            translate(batoms, selected_vertices, displacement,
                      fmax, steps, frame_start=self.nframe)
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


class OB_Force_Field_Modal_Panel(bpy.types.Panel):
    bl_idname = "ob_force_field_modal_panel"
    bl_label = "Open Babel Force field"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Real"

    def draw(self, context):
        obffpanel = context.scene.obffpanel
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.prop(obffpanel, "fmax", expand=True)
        row = box.row()
        row.prop(obffpanel, "steps", expand=True)
        layout = self.layout
        layout.operator("batoms.ob_force_field_operator",
                        text='Force field', icon='FACESEL')


class OBForceFieldProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()

    def Callback_modify_fmax(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = clpanel.transform
        # modify_transform(self.selected_batoms, transform)

    def Callback_modify_steps(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = clpanel.transform
        # modify_transform(self.selected_batoms, transform)
    fmax: FloatProperty(
        name="Max force", default=0.05,
        description="max force", update=Callback_modify_fmax)
    steps: IntProperty(
        name="Max steps", default=5,
        description="Max steps", update=Callback_modify_steps)
