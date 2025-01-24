# import bpy
import importlib
from ..utils.butils import get_preferences_addon

plugin_info = {
    "highlight": ["highlight", "Highlight", "Bhighlight"],
    "cavity": ["cavity", "Cavity", "Bcavity"],
    "crystal_shape": ["crystal_shape", "CrystalShape", "Bcrystalshape"],
    "lattice_plane": ["lattice_plane", "LatticePlane", "Blatticeplane"],
    "isosurface": ["isosurface", "Isosurface", "Bisosurface"],
    "molecular_surface": ["molecular_surface", "MolecularSurface", "Bmolecularsurface"],
    "magres": ["magres", "Magres", "Bmagres"],
    "template": ["template", "Template", "Btemplate"],
}


def enable_plugin():
    # Check if "batoms" add-on is enabled in Blender preferences
    if not get_preferences_addon():
        return
    for key in plugin_info.keys():
        if getattr(get_preferences_addon().preferences, key):
            plugin = importlib.import_module(f".{key}", package=__name__)
            plugin.register_class()


def disable_plugin():
    if not get_preferences_addon():
        return
    for key in plugin_info.keys():
        if getattr(get_preferences_addon().preferences, key):
            plugin = importlib.import_module(f".{key}", package=__name__)
            plugin.unregister_class()
