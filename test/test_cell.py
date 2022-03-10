import bpy
import pytest
from batoms.cell import Bcell
import numpy as np


def test_cell():
    """ """
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    assert isinstance(cell, Bcell)
    cell[2, 2] += 5
    assert cell[2, 2] == 7
    cell.repeat([2, 2, 2])
    assert np.allclose(cell.array, np.array([[4, 0, 0], [0, 4, 0], [0, 0, 14]]))
    cell2 = cell.copy("pt-2")
    assert isinstance(cell2, Bcell)
    reciprocal = cell.reciprocal
    assert np.isclose(reciprocal[0, 0], 1.570796)


if __name__ == "__main__":
    test_cell()
    print("\n Cell: All pass! \n")
