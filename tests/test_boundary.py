import bpy
import pytest
from ase.io import read
from ase.build import bulk
from batoms import Batoms
import numpy as np

def test_boundary():
    bpy.ops.batoms.delete()
    o = Batoms('o', ['O'], [[1, 1, 1]],
               pbc=True,
               cell=[2, 2, 2])
    o.boundary = [1, 0, 0]
    assert len(o.boundary) == 2

def test_boundary_off_origin():
    bpy.ops.batoms.delete()
    au = Batoms("au", from_ase=bulk("Au", cubic=True))
    au.boundary = [1, 1, 0]
    assert np.allclose(au.boundary.positions[0],
                       np.array([0, -2.03999996, 2.03999996]))
    # repeat
    au.translate([0, 0, 2])
    au.boundary = [1, 1, 0]
    assert np.allclose(au.boundary.positions[0],
                       np.array([0, -2.03999996, 4.03999996]))


def test_boundary_oxide():
    from ase.io import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif", ":")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    assert len(tio2.boundary.obj.data.vertices) == 9
    assert np.allclose(tio2.boundary.positions[0],
                       np.array([0, 0, 2.969203]))


def test_boundary_animation():
    from ase.io import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2_10.xyz", ":")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    assert len(tio2.boundary.obj.data.vertices) == 9





if __name__ == "__main__":
    test_boundary()
    print("\n Bondsetting: All pass! \n")
