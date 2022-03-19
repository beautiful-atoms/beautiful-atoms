import bpy
from batoms import Batoms
import pytest


def test_batoms_label_element():
    """Draw label for each atom"""
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label="ch4", formula="CH4")
    ch4 = Batoms('ch4')
    ch4.show_label = 'index'
    ch4.show_label = None
    ch4.show_label = 'species'
    ch4.show_label = 'elements'
    ch4.set_attributes({'charges': np.zeros(5)})
    ch4.show_label = 'charges'



if __name__ == "__main__":
    test_batoms_label_element()
    print("\n Label: All pass! \n")
