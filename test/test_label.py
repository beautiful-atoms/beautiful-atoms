import bpy
from batoms import Batoms
import pytest


def test_batoms_label_element():
    """Create a molecule use GUI ASE"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label="ch4", formula="CH4")
    ch4 = Batoms('ch4')
    ch4.show_label = 1
    ch4.show_label = 0
    ch4.show_label = 2
    ch4.show_label = 3


if __name__ == "__main__":
    test_batoms_label_element()
    print("\n Label: All pass! \n")
