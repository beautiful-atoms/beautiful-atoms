import bpy
import numpy as np
import pytest


def test_settings(ch4):
    """species panel"""
    assert ch4.coll.batoms.ui_list_index_species == 1
    # add
    ch4.species.add("H_1")
    assert ch4.coll.batoms.ui_list_index_species == 2
    # remove
    ch4.species.remove("H_1")
    assert ch4.coll.batoms.ui_list_index_species == 1
    # ops add
    bpy.context.view_layer.objects.active = ch4.obj
    bpy.ops.batoms.species_add(species="H_1")
    assert ch4.coll.batoms.ui_list_index_species == 2
    # ops remove
    bpy.ops.batoms.species_remove(species="H_1")
    assert ch4.coll.batoms.ui_list_index_species == 1


def test_batoms_species(h2o):
    """Setting properties of species"""
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
    h2o.species["H"].materials = {"Metallic": 0.9, "IOR": 1.0, "Roughness": 0.01}
    # elements
    h2o.species["H"].elements = {"H": 0.5, "O": 0.4}
    assert len(h2o["H"].elements) == 3
    h2o["H"].material_style = "mirror"
    # species X
    h2o.replace([0], "X")
    h2o["X"].color = [0.8, 0.8, 0.0, 0.3]


def test_species_color(ch4):
    ch4.species["C"].color = [1, 0, 0, 1]
    ch4.species["C"].material_style = "metallic"
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()
    ch4.species.update()
    assert np.isclose(ch4.species["C"].color, np.array([1, 0, 0, 1])).all()


def test_auto_build_species():
    """auto build species use spglib"""
    pytest.importorskip("spglib")
    from batoms.bio.bio import read

    bpy.ops.batoms.delete()
    magnetite = read("../tests/datas/magnetite.cif")
    magnetite.auto_build_species()
    assert len(magnetite.species) == 5
    assert len(magnetite.bonds.setting) == 10
    assert len(magnetite.polyhedras.setting) == 3


def test_geometry_node_object(tio2):
    """species instances are used in geometry node,
    in batoms, boundary, search_bond. When instances are
    re-build, we need also update the geometry node."""
    tio2.boundary = 0.01
    tio2.bond.show_search = True
    tio2.model_style = 1
    tio2.species["Ti"].color = [1, 1, 0, 1]
    tio2.species.update()
    species_node = tio2.gn_node_group.nodes["Atoms_tio2"].node_tree.nodes[
        "Atoms_tio2_Ti"
    ]
    assert species_node.inputs["Instancer"].default_value is not None
    assert (
        tio2.boundary.gn_node_group.nodes["ObjectInfo_tio2_Ti"]
        .inputs["Object"]
        .default_value
        is not None
    )
    assert (
        tio2.bond.search_bond.gn_node_group.nodes["ObjectInfo_tio2_Ti"]
        .inputs["Object"]
        .default_value
        is not None
    )


def test_color_by_attribute():
    from ase.build import bulk
    from batoms import Batoms
    import numpy as np

    au = bulk("Au", cubic=True) * [5, 5, 5]
    au = Batoms(label="au", from_ase=au)
    z = au.positions[:, 2]
    au.set_attributes({"z_coor": (z - np.min(z)) / (np.max(z) - np.min(z))})
    au.species.color_by_attribute("z_coor")
    assert (
        au.species["Au"].materials["Au"].node_tree.nodes["Attribute"].attribute_name
        == "z_coor"  # noqa W503
    )
    assert len(au.species["Au"].materials["Au"].node_tree.links) > 2
