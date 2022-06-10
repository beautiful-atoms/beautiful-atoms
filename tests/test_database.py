import bpy
import pytest
from batoms import Batoms
from batoms.database.pymatgen import pymatgen_search
from pymatgen.core.structure import Molecule
from pymatgen.core import Lattice, Structure


def test_database_materials_project():
    bpy.ops.batoms.delete()
    # bas = pymatgen_search('your key', "mp-2815")


def test_database_pubchem():
    from batoms.database.pubchem import pubchem_search

    bpy.ops.batoms.delete()
    bas = pubchem_search("31423")


def test_database_pymatgen():
    bpy.ops.batoms.delete()
    # molecule
    co = Molecule(["C", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
    co = Batoms(label="co", from_pymatgen=co)
    # lattice
    fe = Structure(Lattice.cubic(2.8), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    fe = Batoms(label="fe", from_pymatgen=fe)
    assert fe.pbc == [True, True, True]


if __name__ == "__main__":
    test_database_materials_project()
    test_database_pubchem()
    test_database_pymatgen()
    print("\n Database: All pass! \n")
