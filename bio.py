import os
from ase import io
from batoms import Batoms

def read(filename, **kwargs):
    """
    wrapper function for ase.io.read
    """
    atoms = io.read(filename=filename, **kwargs)
    batoms = Batoms(label = os.path.split(filename)[1], atoms=atoms)
    return batoms
