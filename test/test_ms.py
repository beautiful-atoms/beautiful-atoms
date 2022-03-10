import bpy
import pytest
from time import time
from ase.build import bulk, molecule
from ase.io import read
from batoms.pdbparser import read_pdb
from batoms import Batoms


def test_SAS():
    """
    isosurface sasurface 1.4
    isosurface area set 0
    36.107
    """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.mssetting.draw_SAS()
    # area = h2o.mssetting.get_psasa()


def test_SAS_protein():
    """

    isosurface resolution 3 sasurface 1.4

    isosurface area set 0

    Jmol    650.3
    Batoms:  655.6
    freesasa: 657.80
    """
    import numpy as np

    bpy.ops.batoms.delete()
    prot = read_pdb("datas/1ema.pdb")  # 1tim
    prot = Batoms("1ema", from_ase=prot)
    prot.mssetting["1"] = {"resolution": 0.4}
    prot.mssetting.draw_SAS()
    area = prot.mssetting.get_sasa("1")
    # area = prot.mssetting.get_psasa()
    assert abs(area[0] - 14461) < 10000
    prot.selects.add("A", "chain A")
    prot.mssetting.add("2", {"select": "A", "color": [0.8, 0.1, 0.1, 1.0]})
    prot.mssetting.draw_SAS()


def test_SES():
    """
    isosurface sasurface 1.4
    isosurface area set 0
    36.107
    """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.mssetting.draw_SES()
    # area = h2o.mssetting.get_sasa(partial=True)[0]


def test_SES_protein():
    """
    isosurface resolution 3 molecular 1.4

    isosurface area set 0
            area                    volume
    Jmol    7054.5  (esolution 1)
    Jmol    17479.0  (esolution 3)
    Jmol    18518.3  (esolution 5)
    Jmol    19180.0  (esolution 8)
    Batoms: 6811   (resolution 0.4) 24134
    Batoms: 6836   (resolution 0.2) 24065
    Batoms: 7249   (resolution 0.1)  23789
    QuickSES: 6783 (resolution 0.4)
    QuickSES: 7069 (resolution 0.2)
    QuickSES: 7298 (resolution 0.1)  23211
    EDTSurf 7338 (Use slightly different vdw radii)
    """
    bpy.ops.batoms.delete()
    prot = read_pdb("datas/1ema.pdb")
    prot = Batoms("prot", from_ase=prot)
    tstart = time()
    prot.mssetting.draw_SES(parallel=2)
    t = time() - tstart
    assert t < 5
    area = prot.mssetting.get_sesa("1")[0]
    assert abs(area - 8011) < 1000


if __name__ == "__main__":
    test_SAS()
    test_SAS_protein()
    test_SES()
    test_SES_protein()
    print("\n MSsetting: All pass! \n")
