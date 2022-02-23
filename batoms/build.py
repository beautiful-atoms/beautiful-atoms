"""

Wrapper functions for ASE build

"""

from ase import build
from batoms.batoms import Batoms


def molecule(label, symbol, **kwargs):
    """_summary_

    Args:
        label (_type_): _description_
        symbol (_type_): _description_

    Returns:
        _type_: _description_
    """
    atoms = build.molecule(symbol, **kwargs)
    batoms = Batoms(label=label, from_ase=atoms)
    return batoms


def bulk(label, symbol, **kwargs):
    """_summary_

    Args:
        label (_type_): _description_
        symbol (_type_): _description_

    Returns:
        _type_: _description_
    """
    atoms = build.bulk(symbol, **kwargs)
    batoms = Batoms(label=label, from_ase=atoms)
    return batoms


def surface(label, lattice, indices, layers, **kwargs):
    """_summary_

    Args:
        label (_type_): _description_
        lattice (_type_): _description_
        indices (_type_): _description_
        layers (_type_): _description_

    Returns:
        _type_: _description_
    """
    if isinstance(lattice, Batoms):
        lattice = lattice.atoms
    atoms = build.surface(lattice, indices, layers, **kwargs)
    atoms.info.pop('species', None)
    batoms = Batoms(label=label, from_ase=atoms)
    return batoms
