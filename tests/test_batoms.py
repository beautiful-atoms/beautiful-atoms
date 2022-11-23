import bpy
from ase.build import molecule, bulk
from batoms import Batoms
import numpy as np
import pytest

try:
    from pytest_blender.test import pytest_blender_unactive
except ImportError:
    pytest_blender_unactive = False


if '3.1.0' in bpy.app.version_string:
    blender31 = False
else:
    blender31 = True


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_empty():
    """Create an empty Batoms object"""
    from batoms import Batoms
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o")
    assert len(h2o) == 0


def test_batoms_molecule():
    """Create a Batoms object from scratch"""
    bpy.ops.batoms.delete()
    from batoms import Batoms
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    assert len(h2o) == 3


def test_batoms_crystal():
    """Create a Batoms object with cell"""
    bpy.ops.batoms.delete()
    a = 4.08
    positions = [[0, 0, 0], [a / 2, a / 2, 0],
                 [a / 2, 0, a / 2], [0, a / 2, a / 2]]
    au = Batoms(
        label="au",
        species=["Au"] * len(positions),
        positions=positions,
        pbc=True,
        cell=(a, a, a),
    )
    assert au.pbc == [True, True, True]
    assert np.isclose(au.cell[0, 0], a)

def test_batoms_parameters():
    """Create a Batoms object from scratch"""
    bpy.ops.batoms.delete()
    from batoms import Batoms
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
        model_style=1,
        scale=0.5,
    )
    assert h2o.model_style == 1
    assert np.isclose(h2o.scale, 0.5)
    # segments
    assert h2o.segments[0] == 24
    h2o.segments = [6, 6]
    assert h2o.segments[0] == 6


def test_model_style():
    from batoms import Batoms
    from batoms.bio.bio import read
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.1
    tio2.bond.show_search = True
    tio2.model_style = 2
    assert tio2.bond.search_bond.hide == False
    assert tio2.boundary.hide == False
    #
    tio2.polyhedra_style = 2
    assert tio2.bond.hide == True
    assert tio2.bond.search_bond.hide == True
    assert tio2.boundary.hide == False
    #
    tio2.polyhedra_style = 1
    assert tio2.bond.hide == True
    assert tio2.bond.search_bond.hide == False
    assert tio2.boundary.hide == False


def test_batoms_write():
    """Export Batoms to structure file

    Args:
        filename (str): filename and format of output file
    """
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.write("h2o.in")


def test_batoms_transform():
    """Transform: translate, rotate, mirror"""
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.translate([0, 0, 2])
    assert np.isclose(h2o.positions[0], np.array([0, 0, 2.40])).all()
    h2o.mirror("Z")
    assert np.isclose(h2o.positions[0], np.array([0, 0, 1.6])).all()
    # h2o.rotate(90, "Z")
    # assert np.isclose(h2o.positions[0], np.array([0.76, 0, -2.2])).all()


def test_batoms_wrap():
    """Wrap atoms into pbc cell"""
    bpy.ops.batoms.delete()
    a = 4.08
    positions = [[0, 0, 0], [a / 2, a / 2, 0],
                 [a / 2, 0, a / 2], [0, a / 2, a / 2]]
    au = Batoms(
        label="au",
        species=["Au"] * len(positions),
        positions=positions,
        pbc=True,
        cell=(a, a, a),
    )
    assert np.isclose(au[0].position, np.array([0, 0, 0])).all()
    au[0].position = np.array([-0.3, 0, 0])
    # evaluated
    # assert np.isclose(au[0].position, np.array([0.7, 0, 0])).all()


def test_batoms_supercell():
    """make supercell"""
    bpy.ops.batoms.delete()
    au = Batoms("au", from_ase=bulk("Au"))
    # repeat
    au = au*[2, 2, 2]
    # transform
    P = np.array([[2, 3, 0, 5], [0, 1, 0, 5], [0, 0, 1, 0], [0, 0, 0, 1]])
    au = au.transform(P)
    # assert au.cell
    assert len(au) == 16


def test_batoms_occupy():
    """setting occupancy"""
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        species_props={
            "O": {"elements": {"O": 0.8, "N": 0.2}},
            "H": {"elements": {"H": 0.8}},
        },
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.species["O"] = {"elements": {"O": 0.3, "Cl": 0.3}}


def test_batoms_copy():
    """copy batoms"""
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.pbc = True
    h2o.cell = [3, 3, 3]
    h2o_2 = h2o.copy("h2o_2")
    assert isinstance(h2o_2, Batoms)
    assert len(h2o_2) == 3


def test_replace():
    """replace"""
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.replace([1], "C")
    assert len(h2o.species) == 3
    assert h2o[1].species == "C"


def test_batoms_add():
    """Merge two Batoms objects"""
    from batoms import Batoms
    from ase.build import molecule
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    co = Batoms("co", from_ase=molecule("CO"))
    batoms = h2o + co
    assert len(batoms) == 5
    assert len(batoms.species) == 3


def test_from_batoms():
    """Load Batoms"""
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o = Batoms("h2o")
    assert isinstance(h2o, Batoms)
    assert len(h2o.species) == 2
    assert len(h2o) == 3


def test_set_arrays():
    """Set arrays and attributes"""
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.show = [1, 0, 1]
    assert not h2o[1].show
    h2o.scale = 1
    h2o.set_attributes({"scale": np.array([0.3, 0.3, 0.3])})
    # attributes
    h2o.show = 1
    h2o.show = True
    # positions
    positions = h2o.positions
    assert len(positions) == 3
    h2o.positions = h2o.positions + 2
    assert h2o.positions[0][0] == 2
    #
    del h2o[[2]]
    assert len(h2o.arrays["positions"]) == 2

def test_array_attribute():
    from ase.build import bulk
    import numpy as np
    from batoms import Batoms
    from time import time
    bpy.ops.batoms.delete()
    au = bulk('Au')
    # single value
    vel = np.zeros((len(au), 3))
    tensor = np.zeros((len(au), 3, 3))
    au.set_array("vel", vel)
    au.set_array("tensor", tensor)
    au = Batoms('au', from_ase = au)
    au.get_attribute('vel')
    au.get_attribute('tensor')

def test_att_conflict_case1():
    # Case 1: name ending in 0
    from ase.build import bulk
    import numpy as np
    from batoms import Batoms
    from time import time
    bpy.ops.batoms.delete()
    au_ase = bulk('Au') * [2, 2, 2]
    d_arr1 = np.ones((len(au_ase)))
    d_arr2 = np.ones((len(au_ase)))*2
    # delibrately increate additional dimension
    d_arr3 = np.ones((len(au_ase), 2, 2)) * 3
    # Name confusion
    # array with name "d@0" and shape (len(atoms), ) --> attribute name "d@0"
    au_ase.set_array("d@0", d_arr1)
    au_ase.set_array("d@@2", d_arr2)
    # array with name "d" and shape (len(atoms), 2, 2)
    # after that should get d@@@{i} as attribute name
    au_ase.set_array("d", d_arr3)
    au_bl = Batoms('au', from_ase = au_ase)
    assert np.isclose(au_bl.get_attribute('d@0'), d_arr1).all()
    assert np.isclose(au_bl.get_attribute('d@@2'), d_arr2).all()
    assert np.isclose(au_bl.get_attribute('d'), d_arr3).all()

def test_att_conflict_case2():
    # Case 2: huge matrix
    from ase.build import bulk
    import numpy as np
    from batoms import Batoms
    from time import time
    bpy.ops.batoms.delete()
    au_ase = bulk('Au') * [2, 2, 1]
    d_arr1 = np.ones((len(au_ase), 1))
    # delibrately increate additional dimension
    d_arr2 = np.ones((len(au_ase), 15)) * 2
    # Name confusion
    # array with name "d1" and shape (len(atoms), 1) --> attribute name "d10"
    au_ase.set_array("d1", d_arr1)
    # array with name "d" and shape (len(atoms), 15) --> attribute name also "d10"
    # after that should get d@10 as attribute name
    au_ase.set_array("d", d_arr2)
    au_bl = Batoms('au', from_ase = au_ase)
    assert np.isclose(au_bl.get_attribute('d1'), d_arr1).all()
    assert np.isclose(au_bl.get_attribute('d'), d_arr2).all()




###############################
# Patch from TT for Atoms.array
###############################
def test_set_arrays_precision():
    import numpy as np
    from batoms import Batoms
    from batoms.bio.bio import read
    import ase.io
    atoms_ase = ase.io.read("../tests/datas/ch4_int_flag.extxyz")
    # ASE treats additional I-field as np.int32,
    # but on most recent platforms default int is np.int64
    # this may cause error in loading Batoms
    assert atoms_ase.arrays["some_int_flag"].dtype in (np.int32, np.int64)
    assert atoms_ase.arrays["positions"].dtype == np.float64
    # Without patch following will throw error
    atoms_bl = read("../tests/datas/ch4_int_flag.extxyz")
    # Change the array to float32 should also work
    atoms_ase2 = atoms_ase.copy()
    arr = atoms_ase.get_array("some_int_flag")
    atoms_ase2.arrays["some_int_flag"] = np.array([1.0, 2.0, 2.0, 2.0, 2.0], dtype="float32")
    assert atoms_ase2.get_array("some_int_flag").dtype == np.float32
    atoms_bl2 = Batoms("ch4", from_ase=atoms_ase2)
    # get_attribute will return np.int64 and np.float64 explicitly
    assert atoms_bl.attributes["some_int_flag"].dtype == np.int64
    assert atoms_bl2.attributes["some_int_flag"].dtype == np.float64


def test_repeat():
    """Repeat"""
    from ase.build import molecule, bulk
    from batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms(label="h2o", from_ase=h2o)
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.repeat([2, 2, 2])
    assert len(h2o) == 24
    assert np.isclose(h2o.positions[3] - h2o.positions[0], np.array([0, 0, 3])).all()


def test_get_geometry():
    """Test geometry"""
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    angle = h2o.get_angle(1, 0, 2)
    assert np.isclose(angle, 103.9998)
    d = h2o.get_distances(0, 1)
    assert np.isclose(d, 0.96856)
    com = h2o.get_center_of_mass()
    assert np.isclose(com, np.array([0, 0, 0.052531])).all()
    cog = h2o.get_center_of_geometry()
    assert np.isclose(cog, np.array([0, 0, -0.178892])).all()


def test_make_real():
    from ase.build import molecule
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.realize_instances = True
    h2o.realize_instances = False

def test_as_dict():
    from batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(label='au', formula = 'Au')
    au = Batoms("au")
    data = au.as_dict()
    assert 'lattice_plane' not in data.keys()
    # active plugin
    au.lattice_plane
    data = au.as_dict()
    assert 'lattice_plane' in data.keys()



@pytest.mark.skipif(
    blender31,
    reason="Requires Blender >= 3.1.",
)
def test_export_mesh_x3d():
    from ase.build import molecule
    from batoms import Batoms
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    c2h6so = Batoms("c2h6so", from_ase=molecule("C2H6SO"))
    c2h6so.cell = [3, 3, 3]
    c2h6so.model_style = 1
    c2h6so.export_mesh("c2h6so.x3d", with_cell=True, with_bond=True)
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.model_style = 2
    tio2.bond.show_search = True
    tio2.export_mesh("tio2.x3d", with_cell=True,
                     with_polyhedra=True,
                     with_boundary=True,
                     with_search_bond=True,
                     with_bond=True)


def test_export_mesh_obj():
    from ase.build import molecule
    from batoms import Batoms
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    c2h6so = Batoms("c2h6so", from_ase=molecule("C2H6SO"))
    c2h6so.cell = [3, 3, 3]
    c2h6so.model_style = 1
    c2h6so.export_mesh("c2h6so.obj", with_cell=True, with_bond=True)
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.model_style = 2
    tio2.bond.show_search = True
    tio2.export_mesh("tio2.obj", with_cell=True,
                     with_polyhedra=True,
                     with_boundary=True,
                     with_search_bond=True,
                     with_bond=True)


def test_calc_electrostatic_potential():
    from batoms import Batoms
    from ase.io.cube import read_cube_data
    from time import time
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label='c2h6so', formula='C2H6SO')
    c2h6so = Batoms("c2h6so")
    tstart = time()
    c2h6so.auto_assign_charge()
    t = time() - tstart
    print("auto_assign_charge: {:1.2f}".format(t))
    tstart = time()
    c2h6so.calc_electrostatic_potential(c2h6so.positions - 1)
    t = time() - tstart
    print("calc_electrostatic_potential: {:1.2f}".format(t))
    assert t < 5


if __name__ == "__main__":
    test_empty()
    test_batoms_molecule()
    test_batoms_crystal()
    test_batoms_write()
    test_batoms_transform()
    test_batoms_wrap()
    test_batoms_supercell()
    test_batoms_occupy()
    test_batoms_copy()
    test_replace()
    test_batoms_add()
    test_from_batoms()
    test_set_arrays()
    test_att_conflict_case1()
    test_att_conflict_case2()
    test_set_arrays_precision()
    test_repeat()
    test_get_geometry()
    test_make_real()
    print("\n Batoms: All pass! \n")
