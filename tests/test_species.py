import bpy
from batoms import Batoms
import numpy as np
import pytest

try:
    from pytest_blender.test import pytest_blender_unactive
except ImportError:
    pytest_blender_unactive = False


def test_settings():
    """species panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    assert ch4.coll.batoms.ui_list_index_species==1
    # add
    ch4.species.add('H_1')
    assert ch4.coll.batoms.ui_list_index_species==2
    # remove
    ch4.species.remove('H_1')
    assert ch4.coll.batoms.ui_list_index_species==1
    # ops add
    bpy.context.view_layer.objects.active = ch4.obj
    bpy.ops.batoms.species_add(species='H_1')
    assert ch4.coll.batoms.ui_list_index_species==2
    # ops remove
    bpy.ops.batoms.species_remove(species='H_1')
    assert ch4.coll.batoms.ui_list_index_species==1


def test_batoms_species():
    """Setting properties of species"""
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    assert len(h2o.species) == 2
    # default covalent radius
    assert np.isclose(h2o.radius["H"], 0.31)
    # vdw radius
    h2o.radius_style = 1
    assert np.isclose(h2o.radius["H"], 1.2)
    # default Jmol color
    assert np.isclose(h2o["H"].color, np.array([1, 1, 1, 1])).all()
    # VESTA color
    h2o.color_style = 2
    assert np.isclose(h2o["H"].color, np.array([1, 0.8, 0.8, 1])).all()
    # materials
    h2o.species["H"].materials = {
        "Metallic": 0.9, "Specular": 1.0, "Roughness": 0.01}
    # elements
    h2o.species["H"].elements = {"H": 0.5, "O": 0.4}
    assert len(h2o["H"].elements) == 3
    h2o["H"].material_style = "mirror"
    # species X
    h2o.replace([0], "X")
    h2o["X"].color = [0.8, 0.8, 0.0, 0.3]


def test_species_color():
    from batoms import Batoms
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.species["C"].color = [1, 0, 0, 1]
    ch4.species["C"].material_style = "metallic"
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()
    ch4.species.update()
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()


def test_auto_build_species():
    """auto build species use spglib"""
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    magnetite = read("../tests/datas/magnetite.cif")
    magnetite.auto_build_species()
    assert len(magnetite.species) == 5
    assert len(magnetite.bonds.setting) == 10
    assert len(magnetite.polyhedras.setting) == 3

def test_geometry_node_object():
    """species instances are used in geometry node,
    in batoms, boundary, search_bond. When instances are
    re-build, we need also update the geometry node."""
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.bond.show_search = True
    tio2.model_style = 1
    tio2.species["Ti"].color = [1, 1, 0, 1]
    tio2.species.update()
    assert tio2.gnodes.node_group.nodes['ObjectInfo_tio2_Ti'].inputs['Object'].default_value is not None
    assert tio2.boundary.gnodes.node_group.nodes['ObjectInfo_tio2_Ti'].inputs['Object'].default_value is not None
    assert tio2.bond.search_bond.gnodes.node_group.nodes['ObjectInfo_tio2_Ti'].inputs['Object'].default_value is not None

def test_color_by_attribute():
    from ase.build import bulk
    from batoms import Batoms
    import numpy as np
    au = bulk("Au", cubic=True)*[5, 5, 5]
    au = Batoms(label = "au", from_ase = au)
    z = au.positions[:, 2]
    au.set_attributes({"z_coor": (z-np.min(z))/(np.max(z)-np.min(z)) })
    au.species.color_by_attribute("z_coor")
    assert au.species['Au'].materials['Au'].node_tree.nodes['Attribute'].attribute_name == "z_coor"
    assert len(au.species['Au'].materials['Au'].node_tree.links) > 2

if __name__ == "__main__":
    test_batoms_species()
    test_species_color()
    test_geometry_node_object()
    test_auto_build_species()
    print("\n Batoms: All pass! \n")
