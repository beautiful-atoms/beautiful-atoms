import bpy
import pytest
from ase.build import molecule
from batoms import Batoms
import numpy as np


def test_slicebatoms():
    """Setting sliced Batoms"""
    from batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms("CH4")
    # one atoms
    assert ch4[0].species == 'C'
    ch4[1].scale = 1
    assert ch4[1].scale == 1
    ch4[1].show = 0
    assert ch4[1].show == 0
#     ch4[1].species = "Cl"
#     assert ch4[1].species == "Cl"
    assert np.isclose(ch4[1].position,
            np.array([0.62911803,  0.62911803,  0.62911803])).all()
    # > one atoms
    assert ch4[0:2].species[0] == 'C'
    assert np.isclose(ch4[0:2].scale, np.array([0.4, 1])).all()
    ch4[0:2].scale = 2
    assert np.isclose(ch4[0:2].scale, np.array([2, 2])).all()


if __name__ == "__main__":
    test_slicebatoms()
    print("\n SliceBatoms: All pass! \n")
