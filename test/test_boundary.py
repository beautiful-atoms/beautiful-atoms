import bpy
import pytest
from ase.io import read
from ase.build import bulk
from batoms import Batoms


def test_boundary():
    bpy.ops.batoms.delete()
    tio2 = read("datas/tio2.cif")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    assert len(tio2.boundary.obj.data.vertices) == 9


def test_boundary_off_origin():
bpy.ops.batoms.delete()
o = Batoms('o', ['O'], [[1, 1, 1]],
            pbc=True,
            cell=[2, 2, 2])
o.boundary = [1, 0, 0]
print(o.boundary.positions)
o.translate([1, 0, 0])
o.boundary = [1, 0, 0]
print(o.boundary.positions)



def test_boundary_positions():
    bpy.ops.batoms.delete()
    au = Batoms("au", from_ase=bulk("Au", cubic=True))
    au.boundary = [1, 1, 0]
    print(au.boundary.positions)
    print(au.boundary.location)
    # repeat
    au.translate([0, 0, 2])
    au.boundary = [1, 1, 0]
    print(au.boundary.positions)
    print(au.boundary.location)


if __name__ == "__main__":
    test_boundary()
    print("\n Bondsetting: All pass! \n")
