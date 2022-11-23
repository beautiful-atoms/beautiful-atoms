import bpy
import pytest
from batoms.batoms import Batoms
from ase.build import molecule, bulk
from batoms.bio.bio import read
from time import time
import numpy as np

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_bond():
    bpy.ops.batoms.delete()
    from batoms.batoms import Batoms
    from ase.build import molecule
    c2h6so = molecule("C2H6SO")
    c2h6so = Batoms("c2h6so", from_ase=c2h6so)
    c2h6so.model_style = 1
    c2h6so.show = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    c2h6so.model_style = 1
    c2h6so.bond[0].order = 2
    c2h6so.bond.settings["C-H"].order = 2
    c2h6so.model_style = 1
    c2h6so.bond[0].order = 2

def test_settings():
    """key search"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label = 'c2h6so', formula="C2H6SO")
    c2h6so = Batoms('c2h6so')
    # species tuple
    assert c2h6so.bond.settings.find(('C', 'H')) is not None
    # species list
    assert c2h6so.bond.settings.find(['C', 'H']) is not None
    # name str
    assert c2h6so.bond.settings.find('C-H') is not None
    # remove by tuple
    c2h6so.bond.settings.remove(('C', 'H'))
    assert c2h6so.bond.settings.find(('C', 'H')) is None
    # add by tuple
    c2h6so.bond.settings.add(('C', 'H'))
    assert c2h6so.bond.settings.find(('C', 'H')) is not None
    # remove by str
    c2h6so.bond.settings.remove('C-H')
    assert c2h6so.bond.settings.find(('C', 'H')) is None


def test_color():
    bpy.ops.batoms.delete()
    from batoms.batoms import Batoms
    bpy.ops.batoms.molecule_add(label="c2h6so", formula="C2H6SO")
    c2h6so = Batoms("c2h6so")
    c2h6so.model_style = 1
    assert np.isclose(c2h6so.bond.settings.\
        instancers['C-H']['1_1'].data.materials[1].\
            diffuse_color[:], np.array([1, 1, 1, 1])).all()
    c2h6so.bond.settings["C-H"].color1 = [0, 1, 0, 1]
    c2h6so.bond.settings["C-H"].color2 = [0, 1, 1, 1]
    assert np.isclose(c2h6so.bond.settings.\
        instancers['C-H']['1_1'].data.materials[1].\
            diffuse_color[:], np.array([0, 1, 1, 1])).all()


def test_bond_high_order():
    from ase.build import molecule, bulk
    bpy.ops.batoms.delete()
    c6h6 = Batoms("c6h6", from_ase=molecule("C6H6"))
    c6h6.model_style = 1
    c6h6.bond[0].order = 2
    c6h6.bond[2].order = 2
    c6h6.bond[5].order = 2
    c6h6.bond.settings["C-C"].order = 2


def test_bond_performance():
    from ase.build import molecule, bulk
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o = h2o*[10, 10, 10]
    tstart = time()
    h2o.model_style = 1
    t = time() - tstart
    assert t < 5


def test_bond_add():
    bpy.ops.batoms.delete()
    au = bulk("Au")
    au = Batoms("au", from_ase=au)
    au = au*[2, 2, 2]
    assert len(au.bond.settings) == 0
    au.bond.settings.add(("Au", "Au"))
    assert len(au.bond.settings) == 1



def test_bond_search_bond_0():
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.bond.settings[("Ti", "O")].search = 0
    tio2.model_style = 1
    tio2.boundary = 0.01
    tio2.model_style = 1
    assert len(tio2.bond) == 14
    tio2.bond.show_search = True
    if use_cycles:
        set_cycles_res(tio2)
    tio2.get_image([1, -0.3, 0.1], output="bond_search_0.png", **extras)


def test_bond_search_bond_1():
    bpy.ops.batoms.delete()
    mol = read("../tests/datas/anthraquinone.cif")
    mol.boundary = 0.01
    mol.model_style = 1
    mol.bond.show_search = True
    if use_cycles:
        set_cycles_res(mol)
    mol.get_image([1, -0.3, 0.1], output="anthraquinone.png", **extras)


def test_bond_search_bond_2():
    bpy.ops.batoms.delete()
    mof = read("../tests/datas/mof-5.cif")
    mof.boundary = 0.01
    mof.bond.settings[("Zn", "O")].polyhedra = True
    mof.model_style = 1
    if use_cycles:
        set_cycles_res(mof)
    mof.get_image([0, 1, 0], output="mof-5.png", **extras)


def test_bond_search_bond_3():
    from ase.io import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif", ":")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    tio2.model_style = 1
    # assert len(tio2.bond.search_bond) < 20


def test_hydrogen_bond():
    bpy.ops.batoms.delete()
    ch3oh = Batoms(label="ch3oh", from_ase=molecule("CH3OH"))
    ch3oh.bond.settings[("H", "O")].min = 2.0
    ch3oh.bond.settings[("H", "O")].max = 3.0
    ch3oh.bond.settings[("H", "O")].style = "2"
    ch3oh.model_style = 1
    if use_cycles:
        set_cycles_res(ch3oh)
    ch3oh.get_image([1, 0, 0], output="bond-hb.png", **extras)

def test_bond_reload():
    """save to blend file and reload
    """
    import os
    from ase.io import read
    from batoms import Batoms
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif", ":")
    tio2 = Batoms("tio2", from_ase=tio2)
    tio2.boundary = 0.01
    tio2.model_style = 1
    tio2.bond.show_search = True
    cwd = os.getcwd()
    filepath = os.path.join(cwd, "test.blend")
    bpy.ops.wm.save_as_mainfile(filepath=filepath)
    bpy.ops.batoms.delete()
    bpy.ops.wm.open_mainfile(filepath=filepath)
    tio2 = Batoms("tio2")
    assert len(tio2.bond.search_bond) == 43


if __name__ == "__main__":
    test_bond()
    test_bond_high_order()
    test_bond_performance()
    test_bond_add()
    test_bond_search_bond_1()
    test_bond_search_bond_2()
    test_hydrogen_bond()
    print("\n Bond.settings: All pass! \n")
