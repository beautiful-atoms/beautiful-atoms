from bpy.props import (
    PointerProperty,
)
from . import bpy_data
from bpy.types import Collection, Object
from .bpy_data import Base

classes = [
    bpy_data.Belement,
    bpy_data.Bspecies,
    bpy_data.Batom,
    bpy_data.Battribute,
    bpy_data.Bcell,
    bpy_data.Bvolumetric_data,
    bpy_data.Bselect,
    bpy_data.Bboundary,
    bpy_data.BatomsCollection,
    bpy_data.BatomsObject,
]


__all__ = ["register_class", "unregister_class", "Base"]


def register_class():
    from bpy.utils import register_class

    for cls in classes:
        register_class(cls)

    Collection.batoms = PointerProperty(name="Batoms", type=bpy_data.BatomsCollection)
    Object.batoms = PointerProperty(name="Batoms", type=bpy_data.BatomsObject)


def unregister_class():
    from bpy.utils import unregister_class

    for cls in reversed(classes):
        unregister_class(cls)

    del Collection.batoms
    del Object.batoms
