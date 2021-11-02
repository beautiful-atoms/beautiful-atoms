"""
Manipulate atoms and molecules with mouse interactively based on force field methods.


Molecule: Rigid molecules, to constrain all internal degrees of freedom using the RATTLE-type constraints of the The FixBondLengths class to constrain all internal atomic distances.

Intermolecular: atomic radii.


"""
from numpy.lib.function_base import select
import bpy
import numpy as np
from batoms.butils import get_selected_batoms, get_selected_objects
from batoms import Batoms
from batoms.data import covalent_radii, vdw_radii
from batoms.modal.rigid_body import mouse2positions
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

def translate(selected_batom, displacement):
    print(selected_batom)
    pass
def update_positions(atoms, index, mouse_position):
    """

    """
    # Move selected atom with mouse
    atoms[-1].x += atoms.get_cell()[0, 0] / 8
    # run optimize
    optimize(atoms)
    # set new positions of atoms
def optimize(atoms):
    """
    """
    from ase.calculators.emt import EMT
    from ase.optimize import QuasiNewton

    atoms.calc = EMT()
    qn = QuasiNewton(atoms)
    qn.run(fmax=0.05, steps=5)

class Force_Field_Operator(bpy.types.Operator):
    """Move an object with the mouse, example"""
    bl_idname = "object.force_field_operator"
    bl_label = "Force field translate Operator"

    mouse_position: FloatVectorProperty(size = 4, default = (0, 0, 0, 0))
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_batom(self):
        return get_selected_objects('batom')
    def modal(self, context, event):
        """
        """
        if event.type in {'ESC'}:
            return {'CANCELLED'}
        
        elif event.type == 'MOUSEMOVE':
            mouse_position = np.array([event.mouse_x, event.mouse_y, 0, 0])
            if len(self.selected_batoms) != 1:
                self.mouse_position = mouse_position
                return {'RUNNING_MODAL'}
            delta = mouse_position - self.mouse_position
            displacement = mouse2positions(delta, self.viewports_3D.spaces.active.region_3d.    view_matrix)
            self.mouse_position = mouse_position
            translate(self.selected_batom, displacement)
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
            self.mouse_position = np.array([event.mouse_x, event.mouse_y, 0, 0])
            self.first_mouse_x = event.mouse_x
            self.first_value = context.object.location.x
            # bpy.ops.transform.translate('INVOKE_DEFAULT') 
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
        row.prop(ffpanel, "maxf", expand  = True)
        layout = self.layout
        layout.operator("object.force_field_operator", text='Force field', icon='FACESEL')

class ForceFieldProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    def Callback_modify_maxf(self, context):
        clpanel = bpy.context.scene.clpanel
        transform = clpanel.transform
        # modify_transform(self.selected_batoms, transform)
    def Callback_modify_cell(self, context):
        clpanel = bpy.context.scene.clpanel
        cell = [clpanel.cell_a, clpanel.cell_b, clpanel.cell_c]
        # modify_batoms_attr(self.selected_batoms, 'cell', cell)
    maxf: FloatProperty(
        name="Max force", default=0.05,
        description = "max force", update = Callback_modify_maxf)
    
    

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
