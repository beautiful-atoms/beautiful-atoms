# Batoms TOOLBAR  - Addon - Blender 2.9x



###
bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (0, 1, 0),
    "blender": (2, 80, 0),
    "location": "File -> Import -> Batoms (xyz, cif, pdb, ...)",
    "description": "Python module for drawing and rendering ASE (Atomic Simulation Environment) atoms and molecules objects using blender.",
    "warning": "",
    "category": "Import-Export",
}

from batoms.batoms import Batoms
from batoms.batom import Batom


import bpy
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from . import (
        custom_property,
        gui_io,
        gui_batoms
        )

# Register

def menu_func_import_batoms(self, context):
    lay = self.layout
    lay.operator(gui_io.IMPORT_OT_batoms.bl_idname,text="batoms file (xyz, cif, pdb, ...)")


classes = [
        custom_property.Batoms,
        custom_property.Batom,
        custom_property.Bcell,
        custom_property.BBond,
        custom_property.BPolyhedra,
        custom_property.BIsosurface,
        custom_property.BVolume,
        gui_io.IMPORT_OT_batoms,
        gui_batoms.Batoms_PT_prepare,
        gui_batoms.BatomsProperties,
        gui_batoms.AddMolecule,
        gui_batoms.AddBulk,
        gui_batoms.AddAtoms,
    ]
#
def register():
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import_batoms)
    
    for cls in classes:
        bpy.utils.register_class(cls)
    scene = bpy.types.Scene
    scene.blpanel = bpy.props.PointerProperty(type=gui_batoms.BatomsProperties)
    bpy.types.Collection.batoms = bpy.props.PointerProperty(name = 'Batoms', type = custom_property.Batoms)
    bpy.types.Collection.bbond = bpy.props.CollectionProperty(name = 'BBond', type = custom_property.BBond)
    bpy.types.Collection.bpolyhedra = bpy.props.CollectionProperty(name = 'BPolyhedra', type = custom_property.BPolyhedra)
    bpy.types.Collection.bisosurface = bpy.props.CollectionProperty(name = 'BIsosurface', type = custom_property.BIsosurface)
    bpy.types.Object.batom = bpy.props.PointerProperty(name = 'Batom', type = custom_property.Batom)
    bpy.types.Object.bcell = bpy.props.PointerProperty(name = 'Bcell', type = custom_property.Bcell)
    bpy.types.Object.bvolume = bpy.props.PointerProperty(name = 'BVolume', type = custom_property.BVolume)
    

def unregister():

    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import_batoms)

    for cls in classes:
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":

    register()
