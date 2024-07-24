import bpy
from batoms import Batoms
import numpy as np
import pytest

if bpy.app.version_string >= "4.0.0":
    blender40 = False
else:
    blender40 = True


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
    positions = [[0, 0, 0], [a / 2, a / 2, 0], [a / 2, 0, a / 2], [0, a / 2, a / 2]]
    au = Batoms(
        label="au",
        species=["Au"] * len(positions),
        positions=positions,
        pbc=True,
        cell=(a, a, a),
    )
    assert au.pbc == [True, True, True]
    assert np.isclose(au.cell[0, 0], a)


def test_batoms_parameters(h2o):
    """Create a Batoms object from scratch"""
    h2o.model_style = 1
    h2o.scale = 0.5
    assert h2o.model_style == 1
    assert np.isclose(h2o.scale, 0.5)
    # segments
    assert h2o.segments[0] == 24
    h2o.segments = [6, 6]
    assert h2o.segments[0] == 6


def test_model_style(tio2):
    tio2.boundary = 0.1
    tio2.bond.show_search = True
    tio2.model_style = 2
    assert tio2.bond.search_bond.hide is False
    assert tio2.boundary.hide is False
    #
    tio2.polyhedra_style = 2
    assert tio2.bond.hide is True
    assert tio2.bond.search_bond.hide is True
    assert tio2.boundary.hide is False
    #
    tio2.polyhedra_style = 1
    assert tio2.bond.hide is True
    assert tio2.bond.search_bond.hide is False
    assert tio2.boundary.hide is False


def test_batoms_write(h2o):
    """Export Batoms to structure file"""
    from ase.io import read

    h2o.write("h2o.in")
    new_h2o = read("h2o.in")
    assert len(h2o) == len(new_h2o)


def test_batoms_transform(h2o):
    """Transform: translate, rotate, mirror"""
    positions = h2o.positions
    h2o.translate([0, 0, 2])
    # local positions
    assert np.isclose(h2o.positions[0], positions[0]).all()
    # global positions
    assert np.isclose(h2o.global_positions[0], positions[0] + np.array([0, 0, 2])).all()
    # h2o.mirror("Z")
    # assert np.isclose(h2o.global_positions[0], np.array([0, 0, 1.6])).all()
    # h2o.rotate(90, "Z")
    # assert np.isclose(h2o.positions[0], np.array([0.76, 0, -2.2])).all()


def test_batoms_wrap(au):
    """Wrap atoms into pbc cell"""
    assert np.isclose(au[0].position, np.array([0, 0, 0])).all()
    au[0].position = np.array([-0.3, 0, 0])
    # evaluated
    # assert np.isclose(au[0].position, np.array([0.7, 0, 0])).all()


def test_batoms_supercell(au):
    """make supercell"""
    # repeat
    au = au * [2, 2, 2]
    assert len(au) == 32
    # transform
    P = np.array([[2, 3, 0, 5], [0, 1, 0, 5], [0, 0, 1, 0], [0, 0, 0, 1]])
    au = au.transform(P)
    # assert au.cell
    assert len(au) == 64


def test_batoms_occupy(h2o):
    """setting occupancy"""
    h2o.species["O"] = {"elements": {"O": 0.3, "Cl": 0.3}}


def test_batoms_copy(h2o):
    """copy batoms"""
    h2o.pbc = True
    h2o.cell = [3, 3, 3]
    h2o_2 = h2o.copy("h2o_2")
    assert isinstance(h2o_2, Batoms)
    assert len(h2o_2) == 3


def test_replace(h2o):
    """replace"""
    h2o.replace([1], "C")
    assert len(h2o.species) == 3
    assert h2o[1].species[0] == "C"


def test_batoms_add(h2o):
    """Merge two Batoms objects"""
    from batoms import Batoms
    from ase.build import molecule

    co = Batoms("co", from_ase=molecule("CO"))
    batoms = h2o + co
    assert len(batoms) == 5
    assert len(batoms.species) == 3


def test_from_batoms(h2o):
    """Load Batoms"""
    new_h2o = Batoms("h2o")
    assert isinstance(new_h2o, Batoms)
    assert len(new_h2o.species) == 2
    assert len(new_h2o) == 3


def test_set_arrays(h2o):
    """Set arrays and attributes"""
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

    # (nx2) array
    au = bulk("Au")
    array2d = np.zeros((len(au), 2))
    au.set_array("array2d", array2d)
    #
    # (nx3) array
    array3d = np.zeros((len(au), 3))
    au.set_array("array3d", array3d)
    au = Batoms("au", from_ase=au)
    au.get_attribute("array2d")
    au.get_attribute("array3d")
    bpy.ops.batoms.delete()


@pytest.mark.skipif(
    blender40,
    reason="Requires Blender >= 4.0.0",
)
def test_array_attribute_4d():
    from ase.build import bulk
    import numpy as np
    from batoms import Batoms

    # (nx4) array
    au = bulk("Au")
    array4d = np.zeros((len(au), 4))
    au.set_array("array4d", array4d)
    au = Batoms("au", from_ase=au)
    au.get_attribute("array4d")
    bpy.ops.batoms.delete()


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
    atoms_ase2.arrays["some_int_flag"] = np.array(
        [1.0, 2.0, 2.0, 2.0, 2.0], dtype="float32"
    )
    assert atoms_ase2.get_array("some_int_flag").dtype == np.float32
    atoms_bl2 = Batoms("ch4", from_ase=atoms_ase2)
    # get_attribute will return np.int64 and np.float64 explicitly
    assert atoms_bl.attributes["some_int_flag"].dtype == np.int64
    assert atoms_bl2.attributes["some_int_flag"].dtype == np.float64


def test_repeat(h2o):
    """Repeat"""
    import numpy as np

    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.repeat([2, 2, 2])
    assert len(h2o) == 24
    print(h2o.positions)
    assert np.isclose(h2o.positions[3] - h2o.positions[0], np.array([0, 0, 3])).all()


def test_get_geometry(h2o):
    """Test geometry"""
    angle = h2o.get_angle(1, 0, 2)
    assert np.isclose(angle, 103.9998)
    d = h2o.get_distances(0, 1)
    assert np.isclose(d, 0.96856)
    com = h2o.get_center_of_mass()
    assert np.isclose(com, np.array([0, 0, 0.052531])).all()
    cog = h2o.get_center_of_geometry()
    assert np.isclose(cog, np.array([0, 0, -0.178892])).all()


def test_make_real(h2o):
    h2o.realize_instances = True
    h2o.realize_instances = False


def test_as_dict(au):
    data = au.as_dict()
    assert "lattice_plane" not in data.keys()
    # active plugin
    au.lattice_plane
    data = au.as_dict()
    assert "lattice_plane" in data.keys()


def test_batoms_ops():
    from batoms import Batoms

    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(label="au", formula="Au")
    au = Batoms("au")
    assert len(au) == 1


@pytest.mark.skipif(
    not blender40, reason="In Blender>=4.0, export_scene.x3d is not working."
)
def test_export_mesh_x3d(c2h6so):
    from batoms.bio.bio import read

    c2h6so.cell = [3, 3, 3]
    c2h6so.model_style = 1
    c2h6so.export_mesh("c2h6so.x3d", with_cell=True, with_bond=True)
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.model_style = 2
    tio2.bond.show_search = True
    tio2.export_mesh(
        "tio2.x3d",
        with_cell=True,
        with_polyhedra=True,
        with_boundary=True,
        with_search_bond=True,
        with_bond=True,
    )


@pytest.mark.skipif(
    not blender40, reason="In Blender 4.0, export_scene.obj is not working."
)
def test_export_mesh_obj(c2h6so, tio2):
    c2h6so.cell = [3, 3, 3]
    c2h6so.model_style = 1
    c2h6so.export_mesh("c2h6so.obj", with_cell=True, with_bond=True)
    tio2.boundary = 0.01
    tio2.model_style = 2
    tio2.bond.show_search = True
    tio2.export_mesh(
        "tio2.obj",
        with_cell=True,
        with_polyhedra=True,
        with_boundary=True,
        with_search_bond=True,
        with_bond=True,
    )


def test_calc_electrostatic_potential(c2h6so):
    from time import time

    # pytest skip if openbabel is not installed
    pytest.importorskip("openbabel")
    tstart = time()
    c2h6so.auto_assign_charge()
    t = time() - tstart
    print("auto_assign_charge: {:1.2f}".format(t))
    tstart = time()
    c2h6so.calc_electrostatic_potential(c2h6so.positions - 1)
    t = time() - tstart
    print("calc_electrostatic_potential: {:1.2f}".format(t))
    assert t < 5
