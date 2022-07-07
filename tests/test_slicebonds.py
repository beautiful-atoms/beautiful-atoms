import bpy
import pytest
from ase.build import molecule
from batoms import Batoms
import numpy as np


def test_slicebonds():
    """Setting sliced Bonds"""
    from batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms("CH4")
    # one atoms
    assert ch4.bond[0].order == 1
    ch4.bond[1].order = 2
    assert ch4.bond[1].order == 2
    ch4.bond[1].show = 0
    assert ch4.bond[1].show == 0
    # > one bonds
    assert np.isclose(ch4.bond[0:2].order, np.array([1, 2])).all()
    ch4.bond[0:2].order = 2
    assert np.isclose(ch4.bond[0:2].order, np.array([2, 2])).all()


if __name__ == "__main__":
    test_slicebonds()
    print("\n SliceBatoms: All pass! \n")
