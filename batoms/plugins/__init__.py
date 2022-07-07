import bpy
import importlib

plugins = [
    'isosurface',
    'molecular_surface',
    'lattice_plane',
    'crystal_shape',
    'cavity',
    'magres',
    'real_interaction',
]

def enable_plugin():
    for key in plugins:
        if getattr(bpy.context.preferences.addons['batoms'].preferences, key):    
            plugin = importlib.import_module("batoms.plugins.{}".format(key))
            plugin.register_class()   

def disable_plugin():
    for key in plugins:
        if getattr(bpy.context.preferences.addons['batoms'].preferences, key):    
            plugin = importlib.import_module("batoms.plugins.{}".format(key))
            plugin.unregister_class()   