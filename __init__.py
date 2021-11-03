# Batoms TOOLBAR  - Addon - Blender 2.9x



###
bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (0, 1, 0),
    "blender": (2, 83, 0),
    "location": "File -> Import -> Batoms (xyz, cif, pdb, ...)",
    "description": """Python module for drawing and 
rendering ASE (Atomic Simulation Environment) atoms and 
molecules objects using blender.""",
    "warning": "",
    "category": "Import-Export",
}

from batoms.batoms import Batoms
from batoms.batom import Batom


import bpy
from bpy.types import (
    Collection,
    Object,
)
from bpy.props import (
    PointerProperty,
    CollectionProperty,
)

from . import (
    custom_property,
)
from .gui import (
        gui_io,
        gui_batoms,
        gui_batom,
        gui_volume,
        gui_cell,
        gui_bond,
        gui_polyhedra,
        gui_build,
        gui_plane,
        ops_add_molecule,
        )

from .modal import (
    record_selection,
    rigid_body,
    force_field,
)

# Register



classes = [
        custom_property.Batoms,
        custom_property.Batom,
        custom_property.Bcell,
        custom_property.BBond,
        custom_property.BPolyhedra,
        custom_property.BIsosurface,
        custom_property.BVolume,
        custom_property.BPlane,
        gui_io.IMPORT_OT_batoms,
        gui_batoms.Batoms_PT_prepare,
        gui_batoms.BatomsProperties,
        gui_batoms.ReplaceButton,
        gui_batoms.MeasureButton,
        gui_batoms.FragmentateButton,
        ops_add_molecule.AddMolecule,
        gui_build.Build_PT_prepare,
        gui_build.BuildProperties,
        gui_build.AddMolecule,
        gui_build.AddBulk,
        gui_build.AddAtoms,
        gui_batom.Batom_PT_prepare,
        gui_batom.BatomProperties,
        gui_volume.Volume_PT_prepare,
        gui_volume.VolumeProperties,
        gui_volume.AddButton,
        gui_plane.Plane_PT_prepare,
        gui_plane.PlaneProperties,
        gui_plane.AddButton,
        gui_cell.Cell_PT_prepare,
        gui_cell.CellProperties,
        gui_bond.Bond_PT_prepare,
        gui_bond.BondProperties,
        gui_bond.RemoveButton,
        gui_bond.AddButton,
        gui_polyhedra.Polyhedra_PT_prepare,
        gui_polyhedra.PolyhedraProperties,
        record_selection.EDIT_MESH_OT_record_selection,
        rigid_body.Rigid_Body_Operator,
        rigid_body.Rigid_Body_Modal_Panel,
        rigid_body.RigidBodyProperties,
        force_field.Force_Field_Operator,
        force_field.Force_Field_Modal_Panel,
        force_field.ForceFieldProperties
    ]
#
def register():
    bpy.types.TOPBAR_MT_file_import.append(gui_io.menu_func_import_batoms)
    bpy.types.VIEW3D_MT_mesh_add.append(ops_add_molecule.menu_func)
    bpy.types.VIEW3D_MT_select_edit_mesh.append(record_selection.add_edit_mesh_button)
    
    for cls in classes:
        bpy.utils.register_class(cls)
    scene = bpy.types.Scene
    scene.bapanel = PointerProperty(type=gui_batoms.BatomsProperties)
    scene.bupanel = PointerProperty(type=gui_build.BuildProperties)
    scene.btpanel = PointerProperty(type=gui_batom.BatomProperties)
    scene.clpanel = PointerProperty(type=gui_cell.CellProperties)
    scene.bbpanel = PointerProperty(type=gui_bond.BondProperties)
    scene.popanel = PointerProperty(type=gui_polyhedra.PolyhedraProperties)
    scene.plpanel = PointerProperty(type=gui_plane.PlaneProperties)
    scene.vopanel = PointerProperty(type=gui_volume.VolumeProperties)
    scene.rbpanel = PointerProperty(type=rigid_body.RigidBodyProperties)
    scene.ffpanel = PointerProperty(type=force_field.ForceFieldProperties)
    Collection.batoms = PointerProperty(name = 'Batoms', 
                            type = custom_property.Batoms)
    Collection.bbond = CollectionProperty(name = 'BBond', 
                            type = custom_property.BBond)
    Collection.bpolyhedra = CollectionProperty(name = 'BPolyhedra', 
                            type = custom_property.BPolyhedra)
    Collection.bisosurface = CollectionProperty(name = 'BIsosurface', 
                            type = custom_property.BIsosurface)
    Collection.bplane = CollectionProperty(name = 'BPlane', 
                            type = custom_property.BPlane)
    Object.batom = PointerProperty(name = 'Batom', 
                            type = custom_property.Batom)
    Object.bcell = PointerProperty(name = 'Bcell', 
                            type = custom_property.Bcell)
    Object.bbond = PointerProperty(name = 'BBond', 
                            type = custom_property.BBond)
    Object.bpolyhedra = PointerProperty(name = 'BPolyhedra', 
                            type = custom_property.BPolyhedra)
    Object.bisosurface = PointerProperty(name = 'BIsosurface', 
                            type = custom_property.BIsosurface)
    Object.bvolume = PointerProperty(name = 'BVolume', 
                            type = custom_property.BVolume)
    Object.bplane = PointerProperty(name = 'BPlane', 
                            type = custom_property.BPlane)
    

def unregister():

    bpy.types.TOPBAR_MT_file_import.remove(gui_io.menu_func_import_batoms)
    bpy.types.VIEW3D_MT_mesh_add.remove(ops_add_molecule.menu_func)
    bpy.types.VIEW3D_MT_select_edit_mesh.remove(record_selection.add_edit_mesh_button)

    for cls in classes:
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":

    register()
