from ase.build import molecule, bulk
from batoms import Batoms
import numpy as np
from batoms.utils.butils import removeAll
import pytest


def test_empty():
    """Create an empty Batoms object"""
    removeAll()
    h2o = Batoms("h2o")
    assert len(h2o) == 0


def test_batoms_molecule():
    """Create a Batoms object from scratch"""
    removeAll()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    assert len(h2o) == 3


def test_batoms_crystal():
    """Create a Batoms object with cell"""
    removeAll()
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


def test_batoms_species():
    """Setting properties of species"""
    removeAll()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    assert len(h2o.species) == 2
    # default covalent radius
    assert np.isclose(h2o.radius["H"], 0.31)
    # vdw radius
    h2o.species["H"].radius_style = 1
    assert np.isclose(h2o.radius["H"], 1.2)
    # default Jmol color
    assert np.isclose(h2o["H"].color, np.array([1, 1, 1, 1])).all()
    # VESTA color
    h2o.species["H"].color_style = 2
    assert np.isclose(h2o["H"].color, np.array([1, 0.8, 0.8, 1])).all()
    # materials
    h2o.species["H"].materials = {
        "Metallic": 0.9, "Specular": 1.0, "Roughness": 0.01}
    # elements
    h2o.species["H"].elements = {"H": 0.5, "O": 0.4}
    assert len(h2o["H"].elements) == 3
    # species X
    h2o.replace([0], "X")
    h2o["X"].color = [0.8, 0.8, 0.0, 0.3]


def test_batoms_write():
    """Export Batoms to structure file

    Args:
        filename (str): filename and format of output file
    """
    removeAll()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.write("h2o.in")


def test_batoms_transform():
    """Transform: translate, rotate, mirror"""
    removeAll()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    h2o.translate([0, 0, 2])
    assert np.isclose(h2o.positions[0], np.array([0, 0, 2.40])).all()
    h2o.mirror("Z")
    assert np.isclose(h2o.positions[0], np.array([0, 0, 1.6])).all()
    # h2o.rotate(90, 'Z')
    # assert np.isclose(h2o.positions[0], np.array([0.76, 0, -2.2])).all()


def test_batoms_wrap():
    """Wrap atoms into pbc cell"""
    removeAll()
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
    removeAll()
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
    removeAll()
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
    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.pbc = True
    h2o.cell = [3, 3, 3]
    h2o_2 = h2o.copy("h2o_2")
    assert isinstance(h2o_2, Batoms)
    assert len(h2o_2) == 3


def test_replace():
    """replace"""
    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.replace([1], "C")
    assert len(h2o.species) == 3
    assert h2o[1].species == "C"


def test_batoms_add():
    """Merge two Batoms objects"""
    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    co = Batoms("co", from_ase=molecule("CO"))
    batoms = h2o + co
    assert len(batoms) == 5
    assert len(batoms.species) == 3


def test_from_batoms():
    """Load Batoms"""
    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o = Batoms("h2o")
    assert isinstance(h2o, Batoms)
    assert len(h2o.species) == 2
    assert len(h2o) == 3


def test_set_arrays():
    """Set arrays and attributes"""
    removeAll()
    mol = molecule("H2O")
    h2o = Batoms("h2o", from_ase=mol)
    h2o.show = [1, 0, 1]
    assert not h2o[1].show
    h2o.scale = [1, 1, 1]
    h2o.set_attributes({"scale": np.array([0.3, 0.3, 0.3])})
    positions = h2o.positions
    assert len(positions) == 3
    del h2o[[2]]
    assert len(h2o.arrays["positions"]) == 2


def test_repeat():
    """Repeat"""
    removeAll()
    h2o = molecule("H2O")
    h2o = Batoms(label="h2o", from_ase=h2o)
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.repeat([2, 2, 2])
    assert len(h2o) == 24


def test_get_geometry():
    """Test geometry"""
    removeAll()
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
    from batoms.utils.butils import removeAll

    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    h2o.make_real()


if __name__ == "__main__":
    test_empty()
    test_batoms_molecule()
    test_batoms_crystal()
    test_batoms_species()
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
    test_repeat()
    test_get_geometry()
    test_make_real()
    print("\n Batoms: All pass! \n")
