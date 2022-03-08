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
    gui_io,
    gui_batoms,
    gui_batom,
    gui_toolbar,
    gui_volume,
    gui_cell,
    gui_bond,
    gui_polyhedra,
    gui_plane,
    gui_render,
    gui_pymatgen,
    gui_pubchem,
    gui_rscb,
    view3d_mt_batoms_add,
)

from .ops import (add_nanoparticle,
    add_nanotube,
    add_nanoribbon,
    add_object,
    add_surface,
    molecule_edit_atom,
    molecule_edit_bond,
    transform
    )

from .modal import (record_selection,
            rigid_body,
            force_field,
            )

# Register

classes = [
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
    custom_property.Bmssetting,
    custom_property.Batoms_coll,
    custom_property.Batoms_obj,
    add_object.AddMolecule,
    add_object.AddBulk,
    add_object.AddAtoms,
    add_object.AddSurface,
    add_object.AddRootSurface,
    add_object.deleteBatoms,
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
    transform.ApplyCell,
    transform.ApplyTransform,
    transform.ApplyBoundary,
    gui_io.IMPORT_OT_batoms,
    gui_batoms.Batoms_PT_prepare,
    gui_batoms.BatomsProperties,
    gui_batoms.ReplaceButton,
    gui_batoms.MeasureButton,
    gui_batoms.FragmentateButton,
    gui_batom.Batom_PT_prepare,
    gui_batom.BatomProperties,
    gui_bond.Bond_PT_prepare,
    gui_bond.BondProperties,
    gui_bond.RemoveButton,
    gui_bond.AddButton,
    gui_polyhedra.Polyhedra_PT_prepare,
    gui_polyhedra.PolyhedraProperties,
    gui_cell.Cell_PT_prepare,
    gui_cell.CellProperties,
    gui_volume.Volume_PT_prepare,
    gui_volume.VolumeProperties,
    gui_volume.AddButton,
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
    record_selection.EDIT_MESH_OT_record_selection,
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
]
#
addon_keymaps = []


def register():

    # class
    for cls in classes:
        bpy.utils.register_class(cls)

    # menu
    bpy.types.TOPBAR_MT_file_import.append(gui_io.menu_func_import_batoms)
    bpy.types.VIEW3D_MT_add.prepend(view3d_mt_batoms_add.menu_func)

    # tool
    bpy.utils.register_tool(gui_toolbar.AddMolecule, after={"builtin.scale_cage"}, separator=True, group=True)
    bpy.utils.register_tool(gui_toolbar.AddBulk, after={gui_toolbar.AddMolecule.bl_idname})
    bpy.utils.register_tool(gui_toolbar.AddSurface, after={gui_toolbar.AddBulk.bl_idname})
    bpy.utils.register_tool(gui_toolbar.AddAtoms, after={gui_toolbar.AddSurface.bl_idname})
    bpy.utils.register_tool(gui_toolbar.MoleculeEditElement,
        after={"builtin.scale_cage"}, separator=True, group=True)
    bpy.utils.register_tool(gui_toolbar.MolecueEditBond,
        after={gui_toolbar.MoleculeEditElement.bl_idname}, separator=True, group=True)

    # handle the keymap
    wm = bpy.context.window_manager
    # Note that in background mode (no GUI available), keyconfigs are not available either,
    # so we have to check this to avoid nasty errors in background case.
    kc = wm.keyconfigs.addon
    if kc:
        km = wm.keyconfigs.addon.keymaps.new(name='Object Mode', space_type='EMPTY')
        kmi = km.keymap_items.new(gui_toolbar.MoleculeEditElement.bl_idname,
                'RIGHTMOUSE', 'PRESS',
                )
        addon_keymaps.append((km, kmi))

    scene = bpy.types.Scene
    scene.bapanel = PointerProperty(type=gui_batoms.BatomsProperties)
    scene.btpanel = PointerProperty(type=gui_batom.BatomProperties)
    scene.clpanel = PointerProperty(type=gui_cell.CellProperties)
    scene.bbpanel = PointerProperty(type=gui_bond.BondProperties)
    scene.popanel = PointerProperty(type=gui_polyhedra.PolyhedraProperties)
    scene.plpanel = PointerProperty(type=gui_plane.PlaneProperties)
    scene.repanel = PointerProperty(type=gui_render.RenderProperties)
    scene.vopanel = PointerProperty(type=gui_volume.VolumeProperties)
    scene.rbpanel = PointerProperty(type=rigid_body.RigidBodyProperties)
    scene.ffpanel = PointerProperty(type=force_field.ForceFieldProperties)
    scene.pmgpanel = PointerProperty(type=gui_pymatgen.PymatgenProperties)
    scene.pubcpanel = PointerProperty(type=gui_pubchem.PubchemProperties)
    Collection.batoms = PointerProperty(name='Batoms',
                                        type=custom_property.Batoms_coll)
    Object.batoms = PointerProperty(name='Batom',
                                    type=custom_property.Batoms_obj)


def unregister():
    
    # class
    for cls in classes:
        bpy.utils.unregister_class(cls)

    # menu
    bpy.types.TOPBAR_MT_file_import.remove(gui_io.menu_func_import_batoms)
    bpy.types.VIEW3D_MT_add.remove(view3d_mt_batoms_add.menu_func)

    # tool
    bpy.utils.unregister_tool(gui_toolbar.AddMolecule)
    bpy.utils.unregister_tool(gui_toolbar.AddBulk)
    bpy.utils.unregister_tool(gui_toolbar.AddSurface)
    bpy.utils.unregister_tool(gui_toolbar.AddAtoms)
    bpy.utils.unregister_tool(gui_toolbar.MoleculeEditElement)
    bpy.utils.unregister_tool(gui_toolbar.MolecueEditBond)


if __name__ == "__main__":

    register()
