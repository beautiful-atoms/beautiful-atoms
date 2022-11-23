from . import (
    ops_io,
    add_nanoparticle,
    add_nanotube,
    add_nanoribbon,
    add_object,
    add_surface,
    molecule_edit_atom,
    molecule_edit_bond,
    ops_batoms,
    ops_species,
    ops_volumetric_data,
    measure,
    manual_mapping,
    ops_render,
    )

classes = [
    add_object.AddMolecule,
    add_object.AddBulk,
    add_object.AddAtoms,
    add_object.AddSmiles,
    add_surface.BuildSurfaceFCC100,
    add_surface.BuildSurfaceFCC110,
    add_surface.BuildSurfaceFCC111,
    add_surface.BuildSurfaceFCC211,
    add_surface.BuildSurfaceFCC111Root,
    add_surface.BuildSurfaceBCC100,
    add_surface.BuildSurfaceBCC110,
    add_surface.BuildSurfaceBCC111,
    add_surface.BuildSurfaceBCC111Root,
    add_surface.BuildSurfaceHCP0001,
    add_surface.BuildSurfaceHCP10m10,
    add_surface.BuildSurfaceHCP0001Root,
    add_surface.BuildSurfaceDiamond100,
    add_surface.BuildSurfaceDiamond111,
    add_nanotube.BuildNanotube,
    add_nanoribbon.BuildNanoribbon,
    add_nanoparticle.BuildDecahedron,
    add_nanoparticle.BuildIcosahedron,
    add_nanoparticle.BuildOctahedron,
    molecule_edit_atom.MolecueEditElement,
    molecule_edit_bond.MolecueEditBond,
    ops_batoms.BatomsReplace,
    ops_batoms.BatomModify,
    ops_batoms.ApplyCell,
    ops_batoms.ApplyTransform,
    ops_batoms.ApplyBoundary,
    ops_batoms.ApplyModelStyle,
    ops_batoms.ApplyRadiusStyle,
    ops_batoms.ApplyColorStyle,
    ops_batoms.BatomsCopy,
    ops_batoms.ApplyModelStyleSelected,
    ops_batoms.ApplyLabel,
    ops_batoms.BatomsAutoBuildSpecies,
    ops_batoms.AddSurface,
    ops_batoms.AddRootSurface,
    ops_batoms.BatomsJoin,
    ops_batoms.BatomsSeparate,
    ops_batoms.deleteBatoms,
    ops_batoms.deleteSelectedBatoms,
    ops_batoms.BatomsDeleteSelectedAtoms,
    ops_io.IMPORT_OT_batoms,
    ops_io.EXPORT_OT_batoms,
    ops_species.SpeciesAdd,
    ops_species.SpeciesRemove,
    ops_species.SpeciesUpdate,
    ops_species.SpeciesModify,
    ops_volumetric_data.VolumetricDataAdd,
    ops_volumetric_data.VolumetricDataRemove,
    ops_volumetric_data.VolumetricDataCreate,
    ops_render.RenderAdd,
    measure.MeasureButton,
]

manuals = [manual_mapping.batoms_ase_manual_map,manual_mapping.batoms_manual_map]



def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)


def register_manual_map():
    from bpy.utils import register_manual_map
    for cls in manuals:
        register_manual_map(cls)


def unregister_manual_map():
    from bpy.utils import unregister_manual_map
    for cls in reversed(manuals):
        unregister_manual_map(cls)


def register_menu():
    import bpy
    bpy.types.TOPBAR_MT_file_import.append(ops_io.menu_func_import_batoms)
    bpy.types.TOPBAR_MT_file_export.append(ops_io.menu_func_export_batoms)

def unregister_menu():
    import bpy
    bpy.types.TOPBAR_MT_file_import.remove(ops_io.menu_func_import_batoms)
    bpy.types.TOPBAR_MT_file_export.remove(ops_io.menu_func_export_batoms)
