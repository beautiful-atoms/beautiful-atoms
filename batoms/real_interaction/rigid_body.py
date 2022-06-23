"""
Manipulate atoms and molecules with mouse interactively
based on force field methods.


Molecule: Rigid molecules, to constrain all internal
degrees of freedom using the RATTLE-type constraints of the
The FixBondLengths class to constrain all internal atomic distances.

Intermolecular: atomic radii.

"""

import bpy
import numpy as np
from batoms.utils.butils import get_selected_batoms, read_batoms_list
from batoms import Batoms
from batoms.data import covalent_radii
from bpy.props import (StringProperty,
                       BoolProperty,
                       FloatVectorProperty,
                       )
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def translate(selected_batoms, displacement):
    """
    """
    from ase.geometry import get_distances
    displacement = displacement[:3]*0.1
    # for ba
    batoms = Batoms(selected_batoms)
    atoms = batoms.atoms
    radii = covalent_radii[atoms.get_atomic_numbers()]
    n = len(atoms)
    batoms_list = read_batoms_list()
    flag = True
    for name in batoms_list:
        if name == selected_batoms:
            continue
        batoms1 = Batoms(name)
        batoms1.select = False
        atoms1 = batoms1.atoms
        radii1 = covalent_radii[atoms1.get_atomic_numbers()]
        n1 = len(atoms1)
        radii_mat = np.tile(radii, (n1, 1)).T
        radii1_mat = np.tile(radii1, (n, 1))
        radii_sum = radii_mat + radii1_mat
        p = atoms.positions + displacement
        p1 = atoms1.positions
        dvec, dmat = get_distances(p, p1)
        flag_mat = np.less(dmat, radii_sum)
        if flag_mat.any():
            flag = False
    if flag:
        batoms.translate(displacement)


def mouse2positions(delta, mat):
    """
    using rotation matrix
    """
    displacement = delta@mat
    return displacement


class Rigid_Body_Operator(bpy.types.Operator):
    """Rigid body
    """
    bl_idname = "object.rigid_body_operator"
    bl_label = "Rigi body translate Operator"

    mouse_position: FloatVectorProperty(size=4, default=(0, 0, 0, 0))
    previous: StringProperty()

    @property
    def selected_batoms(self):
        return get_selected_batoms()

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
            if len(self.selected_batoms) != 1:
                self.mouse_position = mouse_position
                return {'RUNNING_MODAL'}
            delta = mouse_position - self.mouse_position
            displacement = mouse2positions(
                delta, self.viewports_3D.spaces.active.region_3d.view_matrix)
            translate(self.selected_batoms[0], displacement)
            self.previous = 'MOUSEMOVE'
            return {'RUNNING_MODAL'}
        else:
            mouse_position = np.array([event.mouse_x, event.mouse_y, 0, 0])
            self.mouse_position = mouse_position
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
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}


class Rigid_Body_Modal_Panel(bpy.types.Panel):
    bl_idname = "rigid_body_modal_panel"
    bl_label = "Real Interactive"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Real"
    bl_idname = "RIGIDBODY_PT_Modal"

    def draw(self, context):
        rbpanel = context.scene.rbpanel
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.prop(rbpanel, "fix", expand=True)

        box = layout.box()
        row = box.row()
        row.operator("object.rigid_body_operator",
                     text='Rigid Body', icon='FACESEL')


class RigidBodyProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()

    def Callback_modify_fix(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = clpanel.transform
        # modify_transform(self.selected_batoms, transform)

    def Callback_modify_cell(self, context):
        clpanel = bpy.context.scene.clpanel
        cell = [clpanel.cell_a, clpanel.cell_b, clpanel.cell_c]
        # modify_batoms_attr(self.selected_batoms, 'cell', cell)

    fix: BoolProperty(
        name="fix", default=True,
        description="fix", update=Callback_modify_fix)
