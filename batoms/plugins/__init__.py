# import bpy
import importlib

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
    # if "batoms" not in bpy.context.preferences.addons:
    #     return
    for key in plugin_info.keys():
        # if getattr(bpy.context.preferences.addons["batoms"].preferences, key):
        plugin = importlib.import_module(f".{key}", package=__name__)
        plugin.register_class()


def disable_plugin():
    # if "batoms" not in bpy.context.preferences.addons:
    #     return
    for key in plugin_info.keys():
        # if getattr(bpy.context.preferences.addons["batoms"].preferences, key):
        plugin = importlib.import_module(f".{key}", package=__name__)
        plugin.unregister_class()
