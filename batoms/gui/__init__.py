import bpy
from bpy.props import PointerProperty

from . import (
    gui_batoms,
    gui_batom,
    gui_toolbar,
    gui_cell,
    gui_bond,
    gui_plane,
    gui_render,
    gui_pymatgen,
    gui_pubchem,
    gui_rscb,
    ui_list_species,
    ui_list_bond,
    ui_list_polyhedra,
    ui_list_lattice_plane,
    ui_list_crystal_shape,
    ui_list_isosurface,
    ui_list_cavity,
    ui_list_ms,
    view3d_mt_batoms_add,
)

classes = [
    gui_batoms.Batoms_PT_prepare,
    gui_batoms.BatomsProperties,
    gui_batom.Batom_PT_prepare,
    gui_batom.BatomProperties,
    gui_bond.Bond_PT_prepare,
    gui_bond.BondProperties,
    gui_cell.Cell_PT_prepare,
    gui_cell.CellProperties,
    gui_plane.Plane_PT_prepare,
    gui_plane.PlaneProperties,
    gui_plane.AddButton,
    gui_render.Render_PT_prepare,
    gui_render.RenderProperties,
    gui_pymatgen.Pymatgen_PT_prepare,
    gui_pymatgen.PymatgenProperties,
    gui_pymatgen.Search,
    gui_pubchem.Pubchem_PT_prepare,
    gui_pubchem.PubchemProperties,
    gui_pubchem.Search,
    gui_rscb.RSCB_PT_prepare,
    gui_rscb.RSCBProperties,
    gui_rscb.RSCB_Import,
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
    ui_list_polyhedra.BATOMS_MT_polyhedra_context_menu,
    ui_list_polyhedra.BATOMS_UL_polyhedra,
    ui_list_polyhedra.BATOMS_PT_polyhedra,
    ui_list_lattice_plane.BATOMS_MT_lattice_plane_context_menu,
    ui_list_lattice_plane.BATOMS_UL_lattice_planes,
    ui_list_lattice_plane.BATOMS_PT_lattice_planes,
    ui_list_crystal_shape.BATOMS_MT_crystal_shape_context_menu,
    ui_list_crystal_shape.BATOMS_UL_crystal_shapes,
    ui_list_crystal_shape.BATOMS_PT_crystal_shapes,
    ui_list_isosurface.BATOMS_MT_isosurface_context_menu,
    ui_list_isosurface.BATOMS_UL_isosurface,
    ui_list_isosurface.BATOMS_PT_isosurface,
    ui_list_cavity.BATOMS_MT_cavity_context_menu,
    ui_list_cavity.BATOMS_UL_cavity,
    ui_list_cavity.BATOMS_PT_cavity,
    ui_list_ms.BATOMS_MT_ms_context_menu,
    ui_list_ms.BATOMS_UL_ms,
    ui_list_ms.BATOMS_PT_ms,
]

def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    scene = bpy.types.Scene
    scene.bapanel = PointerProperty(type=gui_batoms.BatomsProperties)
    scene.btpanel = PointerProperty(type=gui_batom.BatomProperties)
    scene.clpanel = PointerProperty(type=gui_cell.CellProperties)
    scene.bbpanel = PointerProperty(type=gui_bond.BondProperties)
    scene.plpanel = PointerProperty(type=gui_plane.PlaneProperties)
    scene.repanel = PointerProperty(type=gui_render.RenderProperties)
    scene.pmgpanel = PointerProperty(type=gui_pymatgen.PymatgenProperties)
    scene.pubcpanel = PointerProperty(type=gui_pubchem.PubchemProperties)
    


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    scene = bpy.types.Scene
    del scene.bapanel
    del scene.btpanel
    del scene.clpanel
    del scene.bbpanel
    del scene.plpanel
    del scene.repanel
    del scene.pmgpanel
    del scene.pubcpanel

def register_menu():
    bpy.types.VIEW3D_MT_add.prepend(view3d_mt_batoms_add.menu_func)

def unregister_menu():
    bpy.types.VIEW3D_MT_add.remove(view3d_mt_batoms_add.menu_func)
    


# handle the keymap
wm = bpy.context.window_manager
# Note that in background mode (no GUI available), keyconfigs are not available either,
# so we have to check this to avoid nasty errors in background case.
kc = wm.keyconfigs.addon

def register_keymap():
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


def unregister_keymap():
    if kc:
        bpy.utils.unregister_tool(gui_toolbar.BatomsTransform)
        bpy.utils.unregister_tool(gui_toolbar.BatomsBoundary)
        bpy.utils.unregister_tool(gui_toolbar.BatomsCell)
        bpy.utils.unregister_tool(gui_toolbar.MoleculeEditElement)
        bpy.utils.unregister_tool(gui_toolbar.MolecueEditBond)

