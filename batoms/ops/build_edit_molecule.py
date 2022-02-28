import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       )
from ase.build import molecule
from ase import Atoms
from batoms import Batoms
import numpy as np

elements = {
    'C': {'species': ['C', 'H', 'H', 'H', 'H'],
          'positions': [[0, 0, 0],
                        [0, 1.024, 0.373],
                        [0, 0, -1.09],
                        [-0.89, -0.51, 0.36],
                        [0.89, -0.51, 0.36]]},
    'N': {'species': ["N", "H", "H", "H"],
          "positions": [[0, 0, 0],
                        [0, 0.94, 0.388],
                        [-0.814, -0.47, 0.388],
                        [0.814, -0.47, 0.388]]},
    'O': {'species': ["O", "H", "H"],
          "positions": [[0, 0, 0],
                        [0, -0.76, 0.6],
                        [0, 0.76, 0.6]]},
    'S': {'species': ["S", "H", "H"],
          "positions": [[0, 0, 0],
                        [0, -0.962, 0.967],
                        [0, 0.962, 0.967]]},
    'F': {'species': ["F"],
          'positions': [[0, 0, 0]]},
    'Cl': {'species': ["Cl"],
           'positions': [[0, 0, 0]]},
    'Br': {'species': ["Br"],
           'positions': [[0, 0, 0]]},
}


def molecule_replace_element(batoms, indices, element):
    """Replace by element

    Args:
        element (str): The target element
    """
    positions = batoms.positions
    bondlists = batoms.bonds.bondlists
    for i in indices:
        mol = elements[element]
        mol['positions'] += positions[i]
        # bond
        ba1 = np.where(bondlists[:, 0] == i)[0]
        ba2 = np.where(bondlists[:, 1] == i)[0]
        ba = ba1 + ba2
        batoms.add_vertices(**mol)
    batoms.delete(indices)
    batoms.model_style = 1


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
