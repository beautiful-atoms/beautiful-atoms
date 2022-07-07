import bpy
import pytest
from time import time
from ase.build import molecule
from batoms.pdbparser import read_pdb
from batoms import Batoms


def test_SAS():
    """ """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.ms.draw()
    print(h2o.ms.setting)
    # area = h2o.ms.get_psasa()


def test_SAS_protein():
    """ """
    import numpy as np

    bpy.ops.batoms.delete()
    prot = read_pdb("../tests/datas/1ema.pdb")  # 1tim
    prot = Batoms("1ema", from_ase=prot)
    prot.ms.setting["1"] = {"resolution": 0.4}
    prot.ms.draw()
    area = prot.ms.get_sasa("1")
    # area = prot.ms.get_psasa()
    assert abs(area[0] - 14461) < 10000
    prot.selects.add("A", "chain A")
    prot.ms.setting.add("2", {"select": "A", "color": [0.8, 0.1, 0.1, 1.0]})
    prot.ms.draw()


def test_SES():
    """ """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.ms.setting["1"].type = "SES"
    h2o.ms.draw()
    # area = h2o.ms.get_sasa(partial=True)[0]


def test_SES_protein():
    """
    compare with jmol, QuickSES, EDTSurf
    isosurface resolution 3 molecular 1.4
    isosurface area set 0
            area                    volume
    """
    bpy.ops.batoms.delete()
    prot = read_pdb("../tests/datas/1ema.pdb")
    prot = Batoms("prot", from_ase=prot)
    tstart = time()
    prot.ms.setting["1"].type = "SES"
    prot.ms.draw()
    t = time() - tstart
    assert t < 5
    area = prot.ms.get_sesa("1")[0]
    assert abs(area - 8011) < 1000


if __name__ == "__main__":
    test_SAS()
    test_SAS_protein()
    test_SES()
    test_SES_protein()
    print("\n MSsetting: All pass! \n")
