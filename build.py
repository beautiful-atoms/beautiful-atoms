"""

Wrapper functions for ASE build

"""

from ase import atom, build
from batoms.batoms import Batoms


def molecule(label, symbol, **kwargs):
    atoms = build.molecule(symbol, **kwargs)
    batoms = Batoms(label = label, atoms = atoms)
    return batoms

def bulk(label, symbol, **kwargs):
    atoms = build.bulk(symbol, **kwargs)
    batoms = Batoms(label = label, atoms = atoms)
    return batoms

def surface(label, lattice, indices, layers, **kwargs):
    if isinstance(lattice, Batoms):
        lattice = lattice.atoms
    atoms = build.surface(lattice, indices, layers, **kwargs)
    atoms.info.pop('species', None)
    batoms = Batoms(label = label, atoms = atoms)
    return batoms
