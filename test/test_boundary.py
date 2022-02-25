import pytest
from ase.io import read
from batoms.utils.butils import removeAll
from batoms import Batoms


def test_boundary():
    removeAll()
    tio2 = read("datas/tio2.cif")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    assert len(tio2.boundary.obj.data.vertices) == 9


if __name__ == "__main__":
    test_boundary()
    print("\n Bondsetting: All pass! \n")
