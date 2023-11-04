import bpy


def test_batoms_measure(ch4):
    """Create a molecule use GUI ASE"""
    bpy.ops.object.mode_set(mode="OBJECT")
    ch4.obj.data.vertices.foreach_set("select", [1, 1, 0, 0, 0])
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.batoms.measure()
