import pytest
from ase.build import molecule
from batoms import Batoms
from batoms.utils.butils import removeAll
import numpy as np


def test_batoms_batom():
    """Setting individual Batom"""
    removeAll()
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    print(h2o[0])
    h2o[1].scale = 1
    h2o[1].show = 0
    h2o[1].show = 1
    h2o[1].species = "Cl"
    assert np.isclose(h2o[1].position, np.array([0.0, 0.76323903, -0.477047])).all()
    assert h2o[1].species == "Cl"


if __name__ == "__main__":
    test_batoms_batom()
    print("\n Batom: All pass! \n")
