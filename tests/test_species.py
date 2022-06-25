import bpy
from batoms import Batoms
import numpy as np
import pytest

try:
    from pytest_blender.test import pytest_blender_unactive
except ImportError:
    pytest_blender_unactive = False


def test_batoms_species():
    """Setting properties of species"""
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    assert len(h2o.species) == 2
    # default covalent radius
    assert np.isclose(h2o.radius["H"], 0.31)
    # vdw radius
    h2o.radius_style = 1
    assert np.isclose(h2o.radius["H"], 1.2)
    # default Jmol color
    assert np.isclose(h2o["H"].color, np.array([1, 1, 1, 1])).all()
    # VESTA color
    h2o.color_style = 2
    assert np.isclose(h2o["H"].color, np.array([1, 0.8, 0.8, 1])).all()
    # materials
    h2o.species["H"].materials = {
        "Metallic": 0.9, "Specular": 1.0, "Roughness": 0.01}
    # elements
    h2o.species["H"].elements = {"H": 0.5, "O": 0.4}
    assert len(h2o["H"].elements) == 3
    h2o["H"].material_style = "mirror"
    # species X
    h2o.replace([0], "X")
    h2o["X"].color = [0.8, 0.8, 0.0, 0.3]


def test_species_color():
    from batoms import Batoms
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.species["C"].color = [1, 0, 0, 1]
    ch4.species["C"].material_style = "metallic"
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()
    ch4.species.update()
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()


def test_auto_build_species():
    """auto build species use spglib"""
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    magnetite = read("../tests/datas/magnetite.cif")
    magnetite.auto_build_species()
    assert len(magnetite.species) == 5
    assert len(magnetite.bonds.setting) == 10
    assert len(magnetite.polyhedras.setting) == 3


if __name__ == "__main__":
    test_batoms_species()
    test_species_color()
    test_auto_build_species()