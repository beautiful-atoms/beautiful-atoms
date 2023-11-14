import bpy


def test_species(ch4):
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.species) == 2
    bpy.ops.batoms.species_update()
    bpy.ops.batoms.species_add(species="O")
    assert len(ch4.species) == 3
    bpy.ops.batoms.species_remove(species="O")
    assert len(ch4.species) == 2
