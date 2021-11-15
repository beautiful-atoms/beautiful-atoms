import bpy
from bpy_extras.object_utils import AddObjectHelper

from bpy.props import (
    StringProperty,
)
from ase.build import molecule
from batoms import Batoms

def add_molecule(label, formula):
    """
    This function takes inputs and returns vertex and face arrays.
    no actual mesh data creation is done here.
    """
    atoms = molecule(formula)
    batoms = Batoms(label = label, atoms = atoms)
    return batoms


class AddMolecule(bpy.types.Operator, AddObjectHelper):
    """Add a simple molecule mesh"""
    bl_idname = "mesh.primitive_molecule_add"
    bl_label = "Add Molecule"
    bl_options = {'REGISTER', 'UNDO'}

    label: StringProperty(
        name = "Label", default='h2o',
        description = "Label")
    formula: StringProperty(
        name = "Formula", default='H2O',
        description = "formula")
    
    
    def execute(self, context):

        batoms = add_molecule(
            self.label,
            self.formula
        )

        
        return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(AddMolecule.bl_idname, icon='MESH_CUBE')


if __name__ == "__main__":
    # test call
    bpy.ops.mesh.primitive_molecule_add()
