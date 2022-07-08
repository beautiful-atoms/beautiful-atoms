import bpy
import pytest
from time import time
from ase.build import molecule
from batoms.pdbparser import read_pdb
from batoms import Batoms
try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}

def test_SAS():
    """
    """
    from batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label = 'h2o', formula = 'H2O')
    h2o = Batoms("h2o")
    h2o.molecular_surface.draw()
    print(h2o.molecular_surface.settings)
    # area = h2o.molecular_surface.get_psasa()

def test_EPM():
    """Electrostatic Potential Maps
    """
    from batoms import Batoms
    from ase.io.cube import read_cube_data
    hartree, atoms = read_cube_data('../tests/datas/h2o-hartree.cube')
    bpy.ops.batoms.delete()
    h2o = Batoms("h2o", from_ase=atoms)
    h2o.volumetric_data['hartree'] = -hartree
    h2o.molecular_surface.draw()
    assert "Vertex Color" not in bpy.data.objects['h2o_1_sas'].data.materials[0].node_tree.nodes
    # color by potential
    h2o.molecular_surface.settings['1'].color_by = 'hartree'
    h2o.molecular_surface.draw()
    if bpy.app.version_string >= '3.1.0':
        assert "Color Attribute" in bpy.data.objects['h2o_1_sas'].data.materials[0].node_tree.nodes
    else:
        assert "Vertex Color" in bpy.data.objects['h2o_1_sas'].data.materials[0].node_tree.nodes
    if use_cycles:
        set_cycles_res(h2o)
    else:
        h2o.render.resolution = [200, 200]
    h2o.get_image([1, 0, 0], output="h2o-EMP.png", **extras)


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
