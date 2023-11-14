import bpy


def test_cell(ch4):
    """Cell panel"""
    import numpy as np

    bpy.context.view_layer.objects.active = ch4.obj
    # model_style
    assert np.isclose(bpy.context.scene.batoms.cell.cell_a0, 0)
    bpy.context.scene.batoms.cell.cell_a0 = 3
    assert np.isclose(bpy.context.scene.batoms.cell.cell_a0, 3)
    assert np.isclose(ch4.cell[0, 0], 3)
