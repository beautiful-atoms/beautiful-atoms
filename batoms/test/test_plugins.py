
def test_pymatgen():
    from batoms.bio import read
    from batoms.plugins.pymatgen import pymatgen_search
    from batoms.butils import removeAll
    removeAll()
    bas = pymatgen_search('your key', "mp-2815")


def test_pubchem():
    from batoms.bio import read
    from batoms.plugins.pubchem import pubchem_search
    from batoms.butils import removeAll
    removeAll()
    bas = pubchem_search("74236")

