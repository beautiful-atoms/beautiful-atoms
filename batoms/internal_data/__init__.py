from . import bpy_data
from .bpy_data import Base

import bpy
from bpy.props import (
                       PointerProperty,
                       CollectionProperty,
                       )

classes = [
    bpy_data.Belement,
    bpy_data.Bspecies,
    bpy_data.Batom,
    bpy_data.Battribute,
    bpy_data.Bcell,
    bpy_data.Bvolume,
    bpy_data.Bselect,
    bpy_data.Bboundary,
    bpy_data.BatomsCollection,
    bpy_data.BatomsObject,
]


import importlib

modules = ['bond',
            'polyhedra',
            'render',
            'ribbon']

def enable_module():
    for key in modules:
        module = importlib.import_module("batoms.{}".format(key))
        module.register_class()   

def disable_module():
    for key in modules:
        module = importlib.import_module("batoms.{}".format(key))
        module.unregister_class()   
        
from bpy.types import Collection, Object

def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    
    Collection.batoms = PointerProperty(name='Batoms',
                                        type=bpy_data.BatomsCollection)
    Object.batoms = PointerProperty(name='Batoms',
                                    type=bpy_data.BatomsObject)
    enable_module()
    

def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    disable_module()
    
    del Collection.batoms
    del Object.batoms
    
