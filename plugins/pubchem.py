"""
"""
from batoms import Batoms

def pubchem_search(cid):
    from ase.data.pubchem import pubchem_search
    # import ssl
    # ssl._create_default_https_context = ssl._create_unverified_context
    data = pubchem_search(cid = cid)
    batoms = Batoms(label = 'cid_%s'%cid, atoms = data.atoms)
    return batoms