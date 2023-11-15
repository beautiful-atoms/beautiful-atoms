import bpy


def test_batoms(ch4):
    """Batoms panel"""
    bpy.ops.object.select_all(action="DESELECT")
    bpy.context.view_layer.objects.active = ch4.obj
    # model_style
    assert bpy.context.scene.batoms.batoms.model_style.upper() == "SPACE-FILLING"
    bpy.context.scene.batoms.batoms.model_style = "Ball-and-stick"
    assert ch4.model_style == 1
    # radius_style
    assert bpy.context.scene.batoms.batoms.radius_style.upper() == "COVALENT"
    bpy.context.scene.batoms.batoms.radius_style = "VDW"
    assert ch4.radius_style == "1"
    # color_style
    assert bpy.context.scene.batoms.batoms.color_style.upper() == "JMOL"
    bpy.context.scene.batoms.batoms.color_style = "CPK"
    assert ch4.color_style == "2"
    # show
    assert bpy.context.scene.batoms.batoms.show is True
    bpy.context.scene.batoms.batoms.show = False
    assert ch4.show[0] == 0
