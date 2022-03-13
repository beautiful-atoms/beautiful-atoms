bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (2, 0, 0),
    "blender": (3, 0, 0),
    "location": "File -> Import -> Batoms (xyz, cif, pdb, ...)",
    "description": """Python module for drawing and
rendering atoms and molecules objects using blender.""",
    "warning": "",
    "category": "Import-Export",
}

from batoms.batoms import Batoms


import bpy
from bpy.types import Collection, Object
from bpy.props import PointerProperty
from . import custom_property
from .gui import (
    gui_batoms,
    gui_batom,
    gui_toolbar,
    gui_cell,
    gui_bond,
    gui_polyhedra,
    gui_plane,
    gui_render,
    gui_pymatgen,
    gui_pubchem,
    gui_rscb,
    ui_list_species,
    ui_list_bond,
    ui_list_lattice_plane,
    ui_list_crystal_shape,
    ui_list_isosurface,
    ui_list_ms,
    view3d_mt_batoms_add,
)

from .ops import (
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
    ops_surface,
    measure,
    manual_mapping,
    ops_bond,
    ops_plane,
    )

from .ops import classes_ops

from .modal import (
            rigid_body,
            force_field,
            )

# manul

manuals = [manual_mapping.batoms_ase_manual_map,manual_mapping.batoms_manual_map]

classes_prop = [
    custom_property.Belement,
    custom_property.Bspecies,
    custom_property.Batom,
    custom_property.Bcell,
    custom_property.Bbond,
    custom_property.Bpolyhedra,
    custom_property.Bisosurface,
    custom_property.Bvolume,
    custom_property.Bplane,
    custom_property.Blight,
    custom_property.Bcamera,
    custom_property.Brender,
    custom_property.Bsheet,
    custom_property.Bhelix,
    custom_property.Bturn,
    custom_property.Bselect,
    custom_property.Bms,
    custom_property.BatomsCollection,
    custom_property.BatomsObject
]

classes_ops = [
    add_object.deleteBatoms,
    add_object.AddMolecule,
    add_object.AddBulk,
    add_object.AddAtoms,
    add_object.AddSurface,
    add_object.AddRootSurface,
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
    ops_batoms.ApplyCell,
    ops_batoms.ApplyTransform,
    ops_batoms.ApplyBoundary,
    ops_io.IMPORT_OT_batoms,
    ops_io.EXPORT_OT_batoms,
    gui_batoms.Batoms_PT_prepare,
    gui_batoms.BatomsProperties,
    measure.MeasureButton,
    gui_batom.Batom_PT_prepare,
    gui_batom.BatomProperties,
    gui_bond.Bond_PT_prepare,
    gui_bond.BondProperties,
    gui_polyhedra.Polyhedra_PT_prepare,
    gui_polyhedra.PolyhedraProperties,
    gui_cell.Cell_PT_prepare,
    gui_cell.CellProperties,
    gui_plane.Plane_PT_prepare,
    gui_plane.PlaneProperties,
    gui_plane.AddButton,
    gui_render.Render_PT_prepare,
    gui_render.RenderProperties,
    gui_render.AddButton,
    gui_pymatgen.Pymatgen_PT_prepare,
    gui_pymatgen.PymatgenProperties,
    gui_pymatgen.Search,
    gui_pubchem.Pubchem_PT_prepare,
    gui_pubchem.PubchemProperties,
    gui_pubchem.Search,
    gui_rscb.RSCB_PT_prepare,
    gui_rscb.RSCBProperties,
    gui_rscb.RSCB_Import,
    rigid_body.Rigid_Body_Operator,
    rigid_body.Rigid_Body_Modal_Panel,
    rigid_body.RigidBodyProperties,
    force_field.Force_Field_Operator,
    force_field.Force_Field_Modal_Panel,
    force_field.ForceFieldProperties,
    view3d_mt_batoms_add.VIEW3D_MT_batoms_add,
    view3d_mt_batoms_add.VIEW3D_MT_surface_add,
    view3d_mt_batoms_add.VIEW3D_MT_nanotube_add,
    view3d_mt_batoms_add.VIEW3D_MT_nanoparticle_add,
    ui_list_species.BATOMS_MT_species_context_menu,
    ui_list_species.BATOMS_UL_species,
    ui_list_species.BATOMS_PT_species,
    ui_list_bond.BATOMS_MT_bond_pair_context_menu,
    ui_list_bond.BATOMS_UL_bond_pairs,
    ui_list_bond.BATOMS_PT_bond_pairs,
    ui_list_lattice_plane.BATOMS_MT_lattice_plane_context_menu,
    ui_list_lattice_plane.BATOMS_UL_lattice_planes,
    ui_list_lattice_plane.BATOMS_PT_lattice_planes,
    ui_list_crystal_shape.BATOMS_MT_crystal_shape_context_menu,
    ui_list_crystal_shape.BATOMS_UL_crystal_shapes,
    ui_list_crystal_shape.BATOMS_PT_crystal_shapes,
    ui_list_isosurface.BATOMS_MT_isosurface_context_menu,
    ui_list_isosurface.BATOMS_UL_isosurface,
    ui_list_isosurface.BATOMS_PT_isosurface,
    ui_list_ms.BATOMS_MT_ms_context_menu,
    ui_list_ms.BATOMS_UL_ms,
    ui_list_ms.BATOMS_PT_ms,
    ops_species.SpeciesAdd,
    ops_species.SpeciesRemove,
    ops_species.SpeciesUpdate,
    ops_species.SpeciesModify,
    ops_bond.BondPairAdd,
    ops_bond.BondPairRemove,
    ops_bond.BondModify,
    ops_bond.BondDraw,
    ops_plane.LatticePlaneAdd,
    ops_plane.LatticePlaneRemove,
    ops_plane.LatticePlaneDraw,
    ops_plane.LatticePlaneModify,
    ops_plane.CrystalShapeAdd,
    ops_plane.CrystalShapeRemove,
    ops_plane.CrystalShapeDraw,
    ops_plane.CrystalShapeModify,
    ops_surface.MSAdd,
    ops_surface.MSRemove,
    ops_surface.MSDraw,
    ops_surface.MSModify,
    ops_surface.IsosurfaceAdd,
    ops_surface.IsosurfaceRemove,
    ops_surface.IsosurfaceDraw,
    ops_surface.IsosurfaceModify,
]
# classes_ops.extend(classes_ops)
#

addon_tools = []


# handle the keymap
wm = bpy.context.window_manager
# Note that in background mode (no GUI available), keyconfigs are not available either,
# so we have to check this to avoid nasty errors in background case.
kc = wm.keyconfigs.addon

def register():

    # class
    for cls in classes_prop:
        bpy.utils.register_class(cls)

    Collection.batoms = PointerProperty(name='Batoms',
                                        type=custom_property.BatomsCollection)
    Object.batoms = PointerProperty(name='Batoms',
                                    type=custom_property.BatomsObject)
    # class
    for cls in classes_ops:
        bpy.utils.register_class(cls)

    # manual
    for manual in manuals:
        bpy.utils.register_manual_map(manual)
    # menu
    bpy.types.TOPBAR_MT_file_import.append(ops_io.menu_func_import_batoms)
    bpy.types.TOPBAR_MT_file_export.append(ops_io.menu_func_export_batoms)
    bpy.types.VIEW3D_MT_add.prepend(view3d_mt_batoms_add.menu_func)
    
    # in background mode, we don't need tool and we can not regester keymap
    if kc:
        # tool
        bpy.utils.register_tool(gui_toolbar.BatomsTransform, after={"builtin.cursor"}, separator=True, group=True)
        bpy.utils.register_tool(gui_toolbar.BatomsBoundary, after={gui_toolbar.BatomsTransform.bl_idname})
        bpy.utils.register_tool(gui_toolbar.BatomsCell, after={gui_toolbar.BatomsBoundary.bl_idname})
        # bpy.utils.register_tool(gui_toolbar.AddAtoms, after={gui_toolbar.AddSurface.bl_idname})
        bpy.utils.register_tool(gui_toolbar.MoleculeEditElement,
            after={"builtin.scale_cage"}, separator=True, group=True)
        bpy.utils.register_tool(gui_toolbar.MolecueEditBond,
            after={gui_toolbar.MoleculeEditElement.bl_idname}, separator=True, group=True)

    scene = bpy.types.Scene
    scene.bapanel = PointerProperty(type=gui_batoms.BatomsProperties)
    scene.btpanel = PointerProperty(type=gui_batom.BatomProperties)
    scene.clpanel = PointerProperty(type=gui_cell.CellProperties)
    scene.bbpanel = PointerProperty(type=gui_bond.BondProperties)
    scene.popanel = PointerProperty(type=gui_polyhedra.PolyhedraProperties)
    scene.plpanel = PointerProperty(type=gui_plane.PlaneProperties)
    scene.repanel = PointerProperty(type=gui_render.RenderProperties)
    scene.rbpanel = PointerProperty(type=rigid_body.RigidBodyProperties)
    scene.ffpanel = PointerProperty(type=force_field.ForceFieldProperties)
    scene.pmgpanel = PointerProperty(type=gui_pymatgen.PymatgenProperties)
    scene.pubcpanel = PointerProperty(type=gui_pubchem.PubchemProperties)
    


def unregister():
    
    # class
    for cls in classes_prop:
        bpy.utils.unregister_class(cls)
    for cls in classes_ops:
        bpy.utils.unregister_class(cls)
    # manual
    for manual in manuals:
        bpy.utils.unregister_manual_map(manual)
    # menu
    bpy.types.TOPBAR_MT_file_import.remove(ops_io.menu_func_import_batoms)
    bpy.types.TOPBAR_MT_file_export.remove(ops_io.menu_func_export_batoms)
    bpy.types.VIEW3D_MT_add.remove(view3d_mt_batoms_add.menu_func)
    
    # tool
    if kc:
        bpy.utils.unregister_tool(gui_toolbar.BatomsTransform)
        bpy.utils.unregister_tool(gui_toolbar.BatomsBoundary)
        bpy.utils.unregister_tool(gui_toolbar.BatomsCell)
        bpy.utils.unregister_tool(gui_toolbar.MoleculeEditElement)
        bpy.utils.unregister_tool(gui_toolbar.MolecueEditBond)


if __name__ == "__main__":

    register()
