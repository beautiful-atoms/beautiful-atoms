"""
Manipulate atoms and molecules with mouse interactively
based on force field methods.


Molecule: Rigid molecules, to constrain all internal
degrees of freedom using the RATTLE-type constraints of the The FixBondLengths class to constrain all internal atomic distances.

Intermolecular: atomic radii.


"""

import bpy
import numpy as np
from batoms.utils.butils import get_selected_batoms, get_selected_objects, get_selected_vertices, object_mode
from batoms import Batoms
from .modal_rigid_body import mouse2positions
from bpy.props import (StringProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       )

from ase.optimize import QuasiNewton, MDMin
try:
    from asap3 import EMT
except:
    from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def translate(selected_vertices, displacement,
              fmax=0.05, steps=5, frame_start=0):
    """
    """
    if len(selected_vertices) == 0:
        return
    displacement = displacement[:3]*0.03
    label = selected_vertices[0][0]
    object_mode()
    # bpy.context.scene.frame_set(frame_start)
    batoms = Batoms(label)
    # Move selected atom with mouse
    atoms = batoms.atoms
    fixed = []
    for label, species, name, index in selected_vertices:
        # batom = Batom(name)
        ind = np.where(atoms.arrays['species'] == species)[0]
        # positions = batom.positions
        atoms.positions[ind[index]] += displacement
        fixed.append(FixAtoms(indices=ind[index]))
    # run optimize
    # print(atoms.positions)
    atoms.constraints = fixed
    optimize(atoms, fmax, steps)
    batoms.set_frames([atoms], frame_start=frame_start)
    # set new positions of atoms
    batoms.positions = atoms
    # batoms.model_style = 1
    # batoms.bondsetting.add(['Al', 'Al'])
    # batoms.draw_bonds()
    # batoms.draw()


def optimize(atoms, fmax=0.05, steps=5):
    """
    """
    atoms.calc = EMT()
    qn = QuasiNewton(atoms)
    qn.run(fmax=fmax, steps=steps)


def add_constraint(atoms, mol_index):
    """
    """
    from ase.constraints import FixAtoms, FixBondLengths
    # RATTLE-type constraints on Molecule.
    atoms.constraints = FixBondLengths([(3 * i + j, 3 * i + (j + 1) % 3)
                                        for i in range(3**3)
                                        for j in [0, 1, 2]])
    #
    # Fix others
    mask = set(range(len(atoms))) - mol_index
    fixlayers = FixAtoms(mask=mask)
    return atoms


class Force_Field_Operator(bpy.types.Operator):
    """Force field"""
    bl_idname = "object.force_field_operator"
    bl_label = "Force field translate Operator"

    mouse_position: FloatVectorProperty(size=4, default=(0, 0, 0, 0))
    nframe: IntProperty(default=0)
    previous: StringProperty()

    @property
    def selected_batoms(self):
        return get_selected_batoms()

    @property
    def selected_batom(self):
        return get_selected_objects('batom')

    @property
    def selected_vertices(self):
        return get_selected_vertices()

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
            selected_vertices = self.selected_vertices
            if len(selected_vertices) != 1:
                self.mouse_position = mouse_position
                return {'RUNNING_MODAL'}
            delta = mouse_position - self.mouse_position
            displacement = mouse2positions(
                delta, self.viewports_3D.spaces.active.region_3d.    view_matrix)
            self.mouse_position = mouse_position
            ffpanel = context.scene.ffpanel
            fmax = ffpanel.fmax
            steps = ffpanel.steps
            translate(selected_vertices, displacement,
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


class Force_Field_Modal_Panel(bpy.types.Panel):
    bl_idname = "force_field_modal_panel"
    bl_label = "Force field"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Real"
    bl_idname = "FORCEFIELD_PT_Modal"

    def draw(self, context):
        ffpanel = context.scene.ffpanel
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.prop(ffpanel, "fmax", expand=True)
        row = box.row()
        row.prop(ffpanel, "steps", expand=True)
        layout = self.layout
        layout.operator("object.force_field_operator",
                        text='Force field', icon='FACESEL')


class ForceFieldProperties(bpy.types.PropertyGroup):
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


def register():
    bpy.utils.register_class(Force_Field_Operator)
    bpy.utils.register_class(Force_Field_Modal_Panel)


def unregister():
    bpy.utils.unregister_class(Force_Field_Operator)
    bpy.utils.unregister_class(Force_Field_Modal_Panel)


if __name__ == "__main__":
    register()

    # test call
    bpy.ops.object.modal_operator('INVOKE_DEFAULT')
