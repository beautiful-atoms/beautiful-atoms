import pytest
from ase.build import bulk
import bpy

def test_position():
    from batoms import Batoms
    from time import time

    bpy.ops.batoms.delete()
    tstart = time()
    au = bulk("Au", cubic=True)
    au = Batoms(label="au", from_ase=au, segments=[6, 6])
    au.repeat([10, 10, 10])
    au.repeat([5, 5, 5])
    t = time() - tstart
    assert t < 100
    print("Repeat time: {:1.2f}".format(t))
    # get position
    tstart = time()
    positions = au.positions
    t = time() - tstart
    assert t < 4
    print("Get positions time: {:1.2f}".format(t))
    # set position
    tstart = time()
    au["Au"].positions = positions
    t = time() - tstart
    assert t < 4
    print("Set positions time: {:1.2f}".format(t))


def test_scatter_and_gather_attribute():
    from ase.build import bulk
    import numpy as np
    from batoms import Batoms
    from time import time
    bpy.ops.batoms.delete()
    au = bulk('Au')
    au *= [5, 5, 5]
    au *= [10, 10, 10]
    # single value
    d0 = np.zeros((len(au)))
    d2 = np.zeros((len(au), 2))
    d22 = np.zeros((len(au), 2, 2))
    au.set_array("d0d", d0)
    au.set_array("d1d", d2)
    au.set_array("d2d", d22)
    au = Batoms('au', from_ase = au)
    tstart = time()
    au.get_attribute('d0d')
    t = time() - tstart
    print("Gatther data for data (1,): {:1.2f}".format(t))
    assert t < 1
    tstart = time()
    au.get_attribute('d1d')
    t = time() - tstart
    print("Gatther data for data (2,): {:1.2f}".format(t))
    assert t < 2
    tstart = time()
    au.get_attribute('d2d')
    t = time() - tstart
    print("Gatther data for data (2, 2): {:1.2f}".format(t))
    assert t < 4


if __name__ == "__main__":
    test_position()
    print("\n Performance: All pass! \n")
