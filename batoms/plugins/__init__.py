import bpy
import importlib


"""
import pathlib
import os
from time import time
cwd = os.path.dirname(pathlib.Path(__file__))
plugin_folders = []

tstart = time()
for f in os.listdir(cwd):
     if os.path.isdir(os.path.join(cwd, f)) and f != '__pycache__':
         plugin_folders.append(f)

plugin_info = {}
for name in plugin_folders:
    info = getattr(importlib.import_module("batoms.plugins.{}".format(name)), 'plugin_info')
    plugin_info[name] = info
print("Read plugin info: {:1.2f}".format(time() - tstart))
"""


plugin_info = {
    'cavity': ['cavity', 'Cavity', 'Bcavity'],
    'crystal_shape': ['crystal_shape', 'CrystalShape', 'Bcrystalshape'],
    'lattice_plane': ['lattice_plane', 'LatticePlane', 'Blatticeplane'],
    'isosurface': ['isosurface', 'Isosurface', 'Bisosurface'],
    'molecular_surface': ['molecular_surface', 'MolecularSurface', 'Bmolecularsurface'],
    'magres': ['magres', 'Magres', 'Bmagres'],
    }

def enable_plugin():
    for key in plugin_info.keys():
        if getattr(bpy.context.preferences.addons['batoms'].preferences, key):    
            plugin = importlib.import_module("batoms.plugins.{}".format(key))
            plugin.register_class()   

def disable_plugin():
    for key in plugin_info.keys():
        if getattr(bpy.context.preferences.addons['batoms'].preferences, key):    
            plugin = importlib.import_module("batoms.plugins.{}".format(key))
            plugin.unregister_class()   