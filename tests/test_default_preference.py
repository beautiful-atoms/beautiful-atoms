"""Test if batoms' default preferences and startup file are correctly loaded
This test assumes the user has initialized the Blender environment with default settings
"""
from filecmp import DEFAULT_IGNORES
import bpy
import pytest
import numpy as np


def test_pref_load():
    """Test several preference settings that are significantly
    from factory Blender
    """
    # Default gradiant for theme
    DEFAULT_GRAD = 0.18823
    DEFAULT_HI_GRAD = 0.2392
    CUSTOM_GRAD = 0.50
    CUSTOM_HI_GRAD = 0.90
    theme = bpy.context.preferences.themes[0]
    gradient = theme.view_3d.space.gradients.gradient[0]
    hi_gradient = theme.view_3d.space.gradients.high_gradient[0]
    # Checking if a lighter theme is set
    assert not np.isclose(hi_gradient, DEFAULT_HI_GRAD, atol=0.08)
    assert not np.isclose(gradient, DEFAULT_GRAD, atol=0.08)
    assert np.isclose(hi_gradient, CUSTOM_HI_GRAD, atol=0.08)
    assert np.isclose(gradient, CUSTOM_GRAD, atol=0.08)
    # Batoms-related
    assert "batoms" in bpy.context.preferences.addons
    addon = bpy.context.preferences.addons["batoms"]
    assert addon.preferences.logging_level == "WARNING"
    return


def test_startup_load():
    # Check Layout
    # Check Scene unit system
    scene = bpy.context.scene
    assert scene.unit_settings.system == "NONE"
    assert scene.unit_settings.scale_length == 1.0
    # Check no default objects exists
    assert len(bpy.data.objects) == 0
    assert len(bpy.data.collections) == 1
    return


if __name__ == "__main__":
    test_pref_load()
    test_startup_load()
    print("\n Batoms: All pass! \n")
