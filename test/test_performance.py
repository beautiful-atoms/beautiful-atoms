import pytest


def test_position():
    from ase.build import bulk
    from batoms import Batoms
    from time import time
    tstart = time()
    au = bulk('Au', cubic = True)
    au = Batoms(label = 'au', atoms = au, segments = [6, 6])
    au.repeat([10, 10, 20])
    au.repeat([5, 5, 5])
    t = time() - tstart
    assert t < 5
    print('Repeat time: {:1.2f}'.format(t))
    # get position
    tstart = time()
    positions = au.positions
    t = time() - tstart
    assert t < 2
    print('Get positions time: {:1.2f}'.format(t))
    # set position
    tstart = time()
    au['Au'].positions = positions
    t = time() - tstart
    assert t < 2
    print('Set positions time: {:1.2f}'.format(t))


if __name__ == '__main__':
    test_position()
    print('\n Performance: All pass! \n')
