import pytest
import numpy as np

def test_lattice_plane():
    """
    """
    from batoms.build import bulk
    from batoms.butils import removeAll
    removeAll()
    au = bulk('au', 'Au', cubic = True)
    au.planesetting[(1, 1, 1)] = {'distance': 3}
    au.planesetting[(1, 1, 0)] = {'distance': 3, 'color': [0.8, 0.1, 0, 0.8]}
    au.draw_lattice_plane()

def test_crystal_shape():
    """
    """
    from batoms.build import bulk
    from batoms.butils import removeAll
    removeAll()
    au = bulk('au', 'Au', cubic = True)
    au.planesetting[(1, 1, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
    au.planesetting[(0, 0, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
    au.draw_crystal_shape()



if __name__ == '__main__':
    test_lattice_plane()
    test_crystal_shape()
    print('\n Bcell: All pass! \n')