import pytest
import numpy as np

def test_lattice_plane():
    """
    """
    from batoms.build import bulk
    from batoms.butils import removeAll
    removeAll()
    au = bulk('au', 'Au', cubic = True)
    au.planesetting[(1, 1, 1)] = {'distance': 3, 'scale': 1.1}
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
    # au.planesetting[(0, 0, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
    au.draw_crystal_shape(origin = au.cell.center)

def test_boundary():
    from batoms.batoms import Batoms
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    h2o = read('/home/xing/ase/batoms/h2o-homo.cube')
    # h2o = read('/home/xing/ase/batoms/fe.cube')
    h2o.planesetting[(0, 0, 1)] = {'distance': 5, 'boundary': True}
    h2o.planesetting[(0, 0, -1)] = {'distance': -4, 'boundary': True}
    h2o.draw_lattice_plane()
    h2o.render.run([0, 0, 1], engine = 'eevee')




if __name__ == '__main__':
    test_lattice_plane()
    test_crystal_shape()
    test_boundary()
    print('\n Bcell: All pass! \n')