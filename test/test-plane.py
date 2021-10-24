import pytest
import numpy as np

def test_bplane():
    """
    """
    from batoms.build import bulk
    from batoms.butils import removeAll
    removeAll()
    pt = bulk('pt', 'Pt', cubic = True)
    pt.planesetting[(1, 1, 1)] = {'indices': (1, 1, 1), 'distance': 0.3}
    pt.planesetting[(1, 1, 0)] = {'indices': (1, 1, 0), 'distance': 0.3, 'color': [0.8, 0.1, 0, 0.8]}
    pt.draw_lattice_plane()

def test_crystal_shape():
    """
    """
from batoms.build import bulk
from batoms.butils import removeAll
removeAll()
pt = bulk('pt', 'Pt', cubic = True)
pt.planesetting[(1, 1, 1)] = {'indices': (1, 1, 1), 'distance': 0.5}
pt.draw_crystal_shape()


if __name__ == '__main__':
    test_bplane()
    print('\n Bcell: All pass! \n')