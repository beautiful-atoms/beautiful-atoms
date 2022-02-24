import pytest
from batoms.utils.butils import removeAll
from batoms.batoms import Batoms
from ase.io import read
from ase.build import molecule
import numpy as np


def test_polyhedra_molecule():
    removeAll()
    ch4 = molecule("CH4")
    mh4 = ch4.copy()
    mh4.translate([2, 2, 0])
    mh4[0].symbol = "N"
    ch4 = ch4 + mh4
    ch4 = Batoms("ch4", from_ase=ch4)
    ch4.model_style = 2
    ch4.bonds.setting[("C", "H")].polyhedra = True
    ch4.model_style = 1


def test_polyhedra_crystal():
    removeAll()
    tio2 = read("datas/tio2.cif")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.model_style = 2
    tio2 = tio2 * [3, 3, 3]
    tio2.pbc = False
    tio2.model_style = 2


def test_polyhedra_setting():
    removeAll()
    ch4 = Batoms("ch4", from_ase=molecule("CH4"))
    ch4.bonds.setting[("C", "H")].polyhedra = True
    ch4.model_style = 2
    ch4.pbc = True
    ch4.cell = [3, 3, 3]
    ch4 = ch4 * [2, 1, 1]
    sel1 = ch4.selects.add("sel1", [0, 1, 2, 3, 4])
    sel1.model_style = 1
    ch4.polyhedras.setting.remove("C")
    assert len(ch4.polyhedras.setting) == 1
    ch4.polyhedras.setting.add("C")
    assert len(ch4.polyhedras.setting) == 2
    ch4.polyhedras.setting["C"] = {"color": [0.8, 0.1, 0.3, 0.3]}
    ch4.model_style = 2
    ch4.render.engine = "workbench"
    ch4.get_image([1, 1, 0], output="polyhedras.png")


if __name__ == "__main__":
    test_polyhedra_molecule()
    test_polyhedra_crystal()
    test_polyhedra_setting()
    print("\n Polyhedra: All pass! \n")
