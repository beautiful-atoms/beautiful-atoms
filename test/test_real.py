import pytest
from batoms.cell import Bcell
import numpy as np

def test_force_field():
    """
    """
    from ase.build import graphene_nanoribbon
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    gnr = graphene_nanoribbon(2, 2, type='armchair', saturated=True,
                                vacuum=3.5)
    gnr.pbc = False
    gnr = Batoms('gnr', atoms = gnr)


def test_force_field_al():
    """
    """
    from batoms.build import bulk
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    al = bulk('au', 'Al')
    al = al*[1, 20, 1]
    al.pbc = False
    

if __name__ == '__main__':
    test_force_field()
    print('\n Bcell: All pass! \n')