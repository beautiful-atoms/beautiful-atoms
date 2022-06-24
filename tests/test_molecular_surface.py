import bpy
import pytest
from time import time
from ase.build import molecule
from batoms.pdbparser import read_pdb
from batoms import Batoms

def test_SAS():
    """
    """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.molecular_surface.draw()
    print(h2o.molecular_surface.settings)
    # area = h2o.molecular_surface.get_psasa()

def test_SES():
    """
    """
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.molecular_surface.settings["1"].type = "SES"
    h2o.molecular_surface.draw()
    # area = h2o.molecular_surface.get_sasa(partial=True)[0]


def test_molecular_surface_ops():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(formula="C2H6SO")
    c2h6so = Batoms('C2H6SO')
    c2h6so.selects.add("1", [2, 4, 6, 7])
    c2h6so.selects.add("2", [3, 5, 8, 9])
    bpy.context.view_layer.objects.active = c2h6so.obj
    bpy.ops.surface.molecular_surface_add(name="2")
    assert len(c2h6so.molecular_surface.settings) == 2
    bpy.ops.surface.molecular_surface_remove(name="1")
    assert len(c2h6so.molecular_surface.settings) == 1
    print(c2h6so.molecular_surface.settings)
    bpy.ops.surface.molecular_surface_add(name="1")
    assert len(c2h6so.molecular_surface.settings) == 2
    c2h6so.molecular_surface.settings["1"].select = "1"
    c2h6so.molecular_surface.settings["2"].select = "2"
    c2h6so.molecular_surface.settings["2"].type = "SES"
    bpy.ops.surface.molecular_surface_draw()


def test_molecule_surface_uilist():
    """molecule_surface_uilist panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    assert ch4.coll.Bmolecularsurface.ui_list_index==0
    bpy.ops.surface.molecular_surface_add(name='2')
    assert ch4.coll.Bmolecularsurface.ui_list_index==1



def test_SAS_protein():
    """
    """
    import numpy as np

    bpy.ops.batoms.delete()
    prot = read_pdb("../tests/datas/1ema.pdb")  # 1tim
    prot = Batoms("1ema", from_ase=prot)
    prot.molecular_surface.settings["1"] = {"resolution": 0.4}
    prot.molecular_surface.draw()
    area = prot.molecular_surface.get_sasa("1")
    # area = prot.molecular_surface.get_psasa()
    assert abs(area[0] - 14461) < 10000
    prot.selects.add("A", "chain A")
    prot.molecular_surface.settings.add("2", {"select": "A", "color": [0.8, 0.1, 0.1, 1.0]})
    prot.molecular_surface.draw()


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
    prot.molecular_surface.settings["1"].type = "SES"
    prot.molecular_surface.draw()
    t = time() - tstart
    assert t < 5
    area = prot.molecular_surface.get_sesa("1")[0]
    assert abs(area - 8011) < 1000



if __name__ == "__main__":
    test_SAS()
    test_SAS_protein()
    test_SES()
    test_SES_protein()
    print("\n MSsetting: All pass! \n")
