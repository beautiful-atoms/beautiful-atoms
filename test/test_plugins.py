from batoms import Batoms
from batoms.plugins.pymatgen import pymatgen_search
from pymatgen.core.structure import Molecule
from pymatgen.core import Lattice, Structure
from batoms.utils.butils import removeAll


def test_plugins_materials_project():
    removeAll()
    # bas = pymatgen_search('your key', "mp-2815")


def test_plugins_pubchem():
    from batoms.plugins.pubchem import pubchem_search

    removeAll()
    bas = pubchem_search("31423")


def test_plugins_pymatgen():
    removeAll()
    # molecule
    co = Molecule(["C", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
    co = Batoms(label="co", from_pymatgen=co)
    # lattice
    fe = Structure(Lattice.cubic(2.8), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    fe = Batoms(label="fe", from_pymatgen=fe)
    assert fe.pbc == [True, True, True]


if __name__ == "__main__":
    test_plugins_materials_project()
    test_plugins_pubchem()
    test_plugins_pymatgen()
    print("\n Plugins: All pass! \n")
