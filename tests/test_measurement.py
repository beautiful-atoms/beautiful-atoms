import bpy
from batoms import Batoms
import pytest


def test_batoms_measure():
    """Create a molecule use GUI ASE"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label="ch4", formula="CH4")
    ch4 = Batoms("ch4")
    bpy.ops.object.mode_set(mode="OBJECT")
    ch4.obj.data.vertices.foreach_set("select", [1, 1, 0, 0, 0])
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.batoms.measure()


if __name__ == "__main__":
    test_batoms_measure()
    print("\n Label: All pass! \n")
