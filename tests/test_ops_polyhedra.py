import bpy


def test_polyhedra(ch4):
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.polyhedra.settings) == 2
    bpy.ops.batoms.polyhedra_draw()
    bpy.ops.batoms.polyhedra_remove(species="C")
    assert len(ch4.polyhedra.settings) == 1
    bpy.ops.batoms.polyhedra_add(species="C")
    assert len(ch4.polyhedra.settings) == 2
