import bpy
import pytest
from batoms import Batoms
from batoms.database.mp import mp_search


def test_database_materials_project():
    bpy.ops.batoms.delete()
    key = None
    if key:
        mp_search(key, "mp-2815")


def test_database_pubchem():
    from batoms.database.pubchem import pubchem_search

    bpy.ops.batoms.delete()
    batoms = pubchem_search("31423")
    assert batoms.label == "cid_31423"


def test_database_pymatgen():
    pytest.importorskip("pymatgen")
    from pymatgen.core.structure import Molecule
    from pymatgen.core import Lattice, Structure

    bpy.ops.batoms.delete()
    # molecule
    co = Molecule(["C", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
    co = Batoms(label="co", from_pymatgen=co)
    # lattice
    fe = Structure(Lattice.cubic(2.8), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    fe = Batoms(label="fe", from_pymatgen=fe)
    assert fe.pbc == [True, True, True]
