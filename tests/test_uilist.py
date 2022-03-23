from numpy import indices
import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest


def test_species():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.species) == 2
    bpy.ops.batoms.species_update()
    bpy.ops.batoms.species_add(species = "O")
    assert len(ch4.species) == 3
    bpy.ops.batoms.species_remove(species = "O")
    assert len(ch4.species) == 2


def test_bond():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    bpy.ops.bond.bond_pair_remove(name = "C-H")
    assert len(ch4.bonds.setting) == 1
    bpy.ops.bond.bond_pair_add(species1 = "C", species2 = "H")
    assert len(ch4.bonds.setting) == 2
    bpy.ops.bond.draw()


def test_polyhedra():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.polyhedras.setting) == 2
    bpy.ops.batoms.polyhedra_draw()
    bpy.ops.batoms.polyhedra_remove(species = "C")
    assert len(ch4.polyhedras.setting) == 1
    bpy.ops.batoms.polyhedra_add(species = "C")
    assert len(ch4.polyhedras.setting) == 2

def test_isosurface():
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube")
    bpy.context.view_layer.objects.active = h2o.obj
    bpy.ops.surface.isosurface_draw()
    assert len(h2o.isosurfaces.setting) == 1
    bpy.ops.surface.isosurface_remove(name = "1")
    print(h2o.isosurfaces.setting)
    assert len(h2o.isosurfaces.setting) == 0
    bpy.ops.surface.isosurface_add(name = "1")
    assert len(h2o.isosurfaces.setting) == 1
    bpy.ops.surface.isosurface_draw()

def test_ms():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(formula = "C2H6SO")
    c2h6so = Batoms('C2H6SO')
    c2h6so.selects.add("1", [2, 4, 6, 7])
    c2h6so.selects.add("2", [3, 5, 8, 9])
    bpy.context.view_layer.objects.active = c2h6so.obj
    bpy.ops.surface.ms_add(name = "2")
    assert len(c2h6so.ms.setting) == 2
    bpy.ops.surface.ms_remove(name = "1")
    assert len(c2h6so.ms.setting) == 1
    print(c2h6so.ms.setting)
    bpy.ops.surface.ms_add(name = "1")
    assert len(c2h6so.ms.setting) == 2
    c2h6so.ms.setting["1"].select = "1"
    c2h6so.ms.setting["2"].select = "2"
    c2h6so.ms.setting["2"].type = "SES"
    bpy.ops.surface.ms_draw()

def test_lattice_plane():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula = "Au")
    au = Batoms('Au')
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.plane.lattice_plane_add(indices = (1, 1, 1))
    au.lattice_plane.setting["1-1-1"].distance = 3.0
    bpy.ops.plane.lattice_plane_add(indices = (1, 0, 0))
    assert len(au.lattice_plane.setting) == 2
    bpy.ops.plane.lattice_plane_draw()
    bpy.ops.plane.lattice_plane_remove(name = "1-0-0")
    assert len(au.lattice_plane.setting) == 1
    print(au.lattice_plane.setting)
    bpy.ops.plane.lattice_plane_draw()


def test_crystal_shape():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula = "Au", cubic = True)
    au = Batoms('Au')
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.plane.crystal_shape_add(indices = (1, 1, 1))
    au.crystal_shape.setting["1-1-1"].symmetry = True
    au.crystal_shape.setting["1-1-1"].distance = 3.0
    bpy.ops.plane.crystal_shape_add(indices = (1, 0, 0))
    assert len(au.crystal_shape.setting) == 2
    bpy.ops.plane.crystal_shape_draw()
    bpy.ops.plane.crystal_shape_remove(name = "1-0-0")
    assert len(au.crystal_shape.setting) == 8
    print(au.crystal_shape.setting)
    bpy.ops.plane.crystal_shape_draw()

if __name__ == "__main__":
    test_species()
    test_bond()
    test_isosurface()
    test_ms()
    test_lattice_plane()
    test_crystal_shape()
    print("\n UIList: All pass! \n")
