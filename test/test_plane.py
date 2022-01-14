import pytest
import numpy as np

def test_lattice_plane():
    """
    """
from ase.build import bulk
from batoms import Batoms
from batoms.butils import removeAll
removeAll()
au = bulk('Au', cubic = True)
au = Batoms('au', aseAtoms = au)
au.planesetting[(1, 1, 1)] = {'distance': 3, 'scale': 1.1}
au.planesetting[(1, 1, 0)] = {'distance': 3, 'color': [0.8, 0.1, 0, 0.8]}
au.planesetting.draw_lattice_plane()
au.draw_cell()
au.get_image([0, 0, 1], output = 'plane-lattice-plane.png')

def test_crystal_shape():
    """
    """
from ase.build import bulk
from batoms import Batoms
from batoms.butils import removeAll
removeAll()
au = bulk('Au', cubic = True)
au = Batoms('au', aseAtoms = au)
au.planesetting[(1, 1, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
# au.planesetting[(0, 0, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
au.planesetting.draw_crystal_shape(origin = au.cell.center)
    au.get_image([0, 0, 1], output = 'plane-crystal.png')

def test_boundary():
    from batoms.batoms import Batoms
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    h2o = read('/home/xing/ase/batoms/h2o-homo.cube')
    h2o.planesetting[(0, 0, 1)] = {'distance': 6, 'boundary': True}
    h2o.planesetting[(0, 0, -1)] = {'distance': -5, 'boundary': True}
    h2o.draw_lattice_plane()
    h2o.get_image([1, 0, 0], output = 'plane-boundary.png')


if __name__ == '__main__':
    test_lattice_plane()
    test_crystal_shape()
    test_boundary()
    print('\n Bcell: All pass! \n')