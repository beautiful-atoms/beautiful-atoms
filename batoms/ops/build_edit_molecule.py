import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       )
from ase.build import molecule, bulk
from ase import Atoms
from batoms.utils.butils import get_selected_batoms, get_selected_vertices
from batoms import Batoms
import numpy as np

hydrogens = {
    'C': molecule('CH4')[0:4],
    'N': molecule('NH3'),
    'O': molecule('H2O')[0:2],
    # 'S': molecule('H2S')[0:2],
    'F': Atoms('F'),
    'Cl': Atoms('Cl'),
    'Br': Atoms('Br'),
}


def molecule_replace_element(batoms, indices, element):
    """Replace by element

    Args:
        element (str): The target element
    """
    for i in indices:
        batoms.replace(i, element)


class MolecueReplaceElement(Operator):
    bl_idname = "batoms.molecule_replace_element"
    bl_label = "Add Element"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Element")

    element: StringProperty(
        name="Element", default='C',
        description="element")
    bond_order: IntProperty(
        name="Bond order", default=1,
        description="bond order")

    def execute(self, context):
        obj = context.object
        if not obj.batoms.flag:
            print('Please select a Batom.')
            return {'FINISHED'}
        batoms = Batoms(label=obj.batoms.label)
        batoms.model_style = 1
        bpy.context.view_layer.objects.active = batoms.obj
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        count = len(obj.data.vertices)
        sel = np.zeros(count, dtype=np.bool)
        obj.data.vertices.foreach_get('select', sel)
        indices = np.where(sel)[0]
        molecule_replace_element(batoms, indices, self.element)
        return {'FINISHED'}
