import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np



def test_polyhedra():
    from ase.build import molecule
    from batoms.batoms import Batoms
    ch4 = Batoms('ch4', atoms = molecule('CH4'))
    ch4.bondsetting['C-H'].polyhedra = True
    ch4.polyhedrasetting['C'].color = [0.8, 0.1, 0.3, 1.0]
    ch4.model_type = 2


if __name__ == '__main__':
    test_polyhedra()
    print('\n Bondsetting: All pass! \n')