import bpy
from batoms import Batoms


def test_batoms_delete():
    """ """
    from batoms import Batoms

    bpy.ops.batoms.delete()
    # delete all
    bpy.ops.batoms.molecule_add(label="nh3", formula="NH3")
    bpy.ops.batoms.molecule_add(label="h2o", formula="H2O")
    bpy.ops.batoms.delete()
    assert "nh3" not in bpy.data.collections
    assert "h2o" not in bpy.data.collections
    # delete by label
    bpy.ops.batoms.molecule_add(label="nh3", formula="NH3")
    bpy.ops.batoms.molecule_add(label="h2o", formula="H2O")
    bpy.ops.batoms.delete(label="h2o")
    assert "h2o" not in bpy.data.collections
    # delete by selection
    bpy.ops.batoms.molecule_add(label="nh3", formula="NH3")
    bpy.ops.batoms.molecule_add(label="h2o", formula="H2O")
    nh3 = Batoms("nh3")
    h2o = Batoms("h2o")
    nh3.obj.select_set(False)
    h2o.obj.select_set(True)
    bpy.ops.batoms.delete_selected_batoms()
    assert "h2o" not in bpy.data.collections
    bpy.ops.batoms.delete()


def test_batoms_apply_model_style(h2o):
    bpy.context.view_layer.objects.active = h2o.obj
    assert h2o.model_style == 0
    # model_style
    bpy.ops.batoms.apply_model_style(model_style="1")
    assert h2o.model_style == 1
    # radius_style
    bpy.ops.batoms.apply_radius_style(radius_style="1")
    assert h2o.radius_style == "1"
    # color_style
    bpy.ops.batoms.apply_color_style(color_style="2")
    assert h2o.color_style == "2"


def test_batoms_join_seperate(h2o, ch4):
    """ """
    ch4.obj.select_set(True)
    h2o.obj.select_set(True)
    # join
    bpy.ops.batoms.join(label="ch4")
    assert len(ch4) == 8
    # separate
    ch4.separate()
    assert len(ch4) == 5
    assert len(h2o) == 3
    # join with another name
    ch4.obj.select_set(True)
    h2o.obj.select_set(True)
    bpy.ops.batoms.join(label="a")
    a = Batoms("a")
    assert len(a) == 8
    # separate
    a.separate()
    assert len(ch4) == 5
    assert len(h2o) == 3


def test_batoms_apply_label(h2o):
    bpy.context.view_layer.objects.active = h2o.obj
    bpy.ops.batoms.apply_label(label="elements")
    assert h2o.show_label == "elements"
    bpy.ops.batoms.apply_label(label="")
    assert h2o.show_label == ""


def test_ase_molecule():
    """Create a molecule use GUI ASE"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label="nh3", formula="NH3")
    nh3 = Batoms("nh3")
    assert len(nh3) == 4
    bpy.ops.batoms.delete()


def test_ase_bulk():
    """Create an bulk use GUI ASE"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(label="au", formula="Au")
    au = Batoms("au")
    assert len(au) == 1
    # surface
    au.obj.select_set(True)
    bpy.ops.batoms.surface_add()
    bpy.ops.batoms.delete()


def test_ase_surface():
    """Create an surface use GUI ASE"""
    bpy.ops.batoms.delete()
    # fcc
    bpy.ops.surface.fcc100_add(label="au100")
    au100 = Batoms("au100")
    assert len(au100) == 4
    bpy.ops.surface.fcc110_add(label="au110")
    au110 = Batoms("au110")
    assert len(au110) == 4
    bpy.ops.surface.fcc111_add(label="au111")
    au111 = Batoms("au111")
    assert len(au111) == 4
    bpy.ops.surface.fcc211_add(label="au211")
    au211 = Batoms("au211")
    assert len(au211) == 12
    bpy.ops.surface.fcc111_root_add(label="au111_root")
    au111_root = Batoms("au111_root")
    assert len(au111_root) == 12
    # bcc
    bpy.ops.surface.bcc100_add(label="fe100")
    fe100 = Batoms("fe100")
    assert len(fe100) == 4
    bpy.ops.surface.bcc110_add(label="fe110")
    fe110 = Batoms("fe110")
    assert len(fe110) == 4
    bpy.ops.surface.bcc111_add(label="fe111")
    fe111 = Batoms("fe111")
    assert len(fe111) == 4
    bpy.ops.surface.bcc111_root_add(label="fe111_root")
    fe111_root = Batoms("fe111_root")
    assert len(fe111_root) == 12
    # hcp
    bpy.ops.surface.hcp0001_add(label="ti0001")
    ti0001 = Batoms("ti0001")
    assert len(ti0001) == 4
    bpy.ops.surface.hcp10m10_add(label="ti10m10")
    ti10m10 = Batoms("ti10m10")
    assert len(ti10m10) == 16
    bpy.ops.surface.hcp0001_root_add(label="ti0001_root")
    ti0001_root = Batoms("ti0001_root")
    assert len(ti0001_root) == 12
    # diamond
    bpy.ops.surface.diamond100_add(label="c100")
    c100 = Batoms("c100")
    assert len(c100) == 4
    bpy.ops.surface.diamond111_add(label="c111")
    c111 = Batoms("c111")
    assert len(c111) == 4
    bpy.ops.batoms.delete()


# ==============================================
# Below for edit mode
# ==============================================
def test_batoms_apply_model_style_selected(h2o):
    import numpy as np
    from batoms import Batoms

    bpy.ops.surface.fcc111_add(label="au111", symbol="Au", size=(1, 1, 4))
    au111 = Batoms("au111")
    au111 += h2o
    au111.model_style = 0
    bpy.context.view_layer.objects.active = au111.obj
    # only select h2o
    au111.obj.data.vertices.foreach_set("select", [0, 0, 0, 0, 1, 1, 1])
    # change model_style for selected atoms
    bpy.ops.batoms.apply_model_style_selected(model_style="1")
    assert au111.get_attribute("model_style")[0] == 0
    assert au111.get_attribute("model_style")[-1] == 1
    assert np.isclose(au111.get_attribute("scale")[0], 1)
    assert np.isclose(au111.get_attribute("scale")[-1], 0.4)
    bpy.ops.batoms.delete()
