import bpy
import pytest
from batoms.utils.butils import removeAll
from batoms.batoms import Batoms
from ase.build import molecule, bulk
from batoms.bio.bio import read
from time import time

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_bonds():
    bpy.ops.batoms.delete()
    c2h6so = molecule("C2H6SO")
    c2h6so = Batoms("c2h6so", from_ase=c2h6so)
    c2h6so.model_style = 1
    c2h6so.show = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    c2h6so.model_style = 1
    c2h6so.bonds[0].order = 2
    c2h6so.bonds.setting["C-H"].color1 = [0, 1, 0, 1]
    c2h6so.bonds.setting["C-H"].color2 = [0, 1, 1, 1]
    c2h6so.bonds.setting["C-H"].order = 2
    c2h6so.model_style = 1
    c2h6so.bonds[0].order = 2


def test_bonds_high_order():
    bpy.ops.batoms.delete()
    c6h6 = Batoms("c6h6", from_ase=molecule("C6H6"))
    c6h6.model_style = 1
    c6h6.bonds[0].order = 2
    c6h6.bonds[2].order = 2
    c6h6.bonds[5].order = 2
    c6h6.bonds.setting["C-C"].order = 2


def test_bonds_performance():
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o = h2o*[10, 10, 10]
    tstart = time()
    h2o.model_style = 1
    t = time() - tstart
    assert t < 5


def test_bonds_add():
    bpy.ops.batoms.delete()
    au = bulk("Au")
    au = Batoms("au", from_ase=au)
    au = au*[2, 2, 2]
    assert len(au.bonds.setting) == 0
    au.bonds.setting.add(("Au", "Au"))
    assert len(au.bonds.setting) == 1


def test_bonds_search_bond_1():
    bpy.ops.batoms.delete()
    mol = read("../tests/datas/anthraquinone.cif")
    mol.boundary = 0.01
    mol.model_style = 1
    mol.bonds.show_search = True
    if use_cycles:
        set_cycles_res(mol)
    mol.get_image([1, -0.3, 0.1], output="anthraquinone.png", **extras)


def test_bonds_search_bond_2():
    bpy.ops.batoms.delete()
    mof = read("../tests/datas/mof-5.cif")
    mof.boundary = 0.01
    mof.bonds.setting[("Zn", "O")].polyhedra = True
    mof.model_style = 1
    if use_cycles:
        set_cycles_res(mof)
    mof.get_image([0, 1, 0], output="mof-5.png", **extras)


def test_bonds_search_bond_3():
    from ase.io import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif", ":")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    tio2.model_style = 1
    # assert len(tio2.bonds.search_bond) < 20


def test_hydrogen_bond():
    bpy.ops.batoms.delete()
    ch3oh = Batoms(label="ch3oh", from_ase=molecule("CH3OH"))
    ch3oh.bonds.setting[("H", "O")].min = 2.0
    ch3oh.bonds.setting[("H", "O")].max = 3.0
    ch3oh.bonds.setting[("H", "O")].style = "2"
    ch3oh.model_style = 1
    if use_cycles:
        set_cycles_res(ch3oh)
    ch3oh.get_image([1, 0, 0], output="bonds-hb.png", **extras)


if __name__ == "__main__":
    test_bonds()
    test_bonds_high_order()
    test_bonds_performance()
    test_bonds_add()
    test_bonds_search_bond_1()
    test_bonds_search_bond_2()
    test_hydrogen_bond()
    print("\n Bonds.setting: All pass! \n")
