import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_voronoi():
from ase.build import bulk
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
au = bulk('Au', cubic = True)
au = Batoms('au', atoms = au)
au = au*[4, 4, 4]
au.scale = 0.1
au.draw_voronoi()

def test_voronoi():
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
c2h6 = molecule('C2H6')
c2h6.center(3.0)
c2h6 = Batoms('c2h6', atoms = c2h6)
c2h6.scale = 0.1
c2h6.draw_voronoi()
    # au.get_image([0, 0, 1], engine = 'eevee')



if __name__ == '__main__':
    test_diff()
    print('\n Bondsetting: All pass! \n')