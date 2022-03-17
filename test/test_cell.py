import bpy
import pytest
from batoms.cell import Bcell
import numpy as np


def test_cell():
    """ """
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    assert isinstance(cell, Bcell)


def test_cell_transform():
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    cell[2, 2] += 5
    assert cell[2, 2] == 7
    cell.translate([0, 0, 2])
    assert cell[2, 2] == 7
    assert cell.positions[3, 2] == 9
    # cell.rotate(np.pi/4, 'Z')
    # assert cell[2, 2] == 7
    # assert cell.positions[3, 2] == 9


def test_repeat():
    from batoms.cell import Bcell
    import numpy as np
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    # cell.rotate(np.pi/4, 'Z')
    # cell.repeat([2, 2, 2])
    # assert np.allclose(cell.array, np.array(
    #     [[4, 0, 0], [0, 4, 0], [0, 0, 4]]))


def test_translate_repeat():
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    cell.translate([0, 0, 2])
    cell.repeat([2, 2, 2])
    assert np.allclose(cell.array, np.array(
        [[4, 0, 0], [0, 4, 0], [0, 0, 4]]))


def test_copy():
    bpy.ops.batoms.delete()
    cell = Bcell(label="pt", array=[2, 2, 2])
    cell2 = cell.copy("pt-2")
    assert isinstance(cell2, Bcell)
    reciprocal = cell.reciprocal
    assert np.isclose(reciprocal[0, 0], 3.14159)


if __name__ == "__main__":
    test_cell()
    test_cell_transform()
    test_repeat()
    test_translate_repeat()
    test_copy()
    print("\n Cell: All pass! \n")
