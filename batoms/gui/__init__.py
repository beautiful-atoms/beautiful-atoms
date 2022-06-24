import bpy
from bpy.props import PointerProperty

from . import (
    gui_batoms,
    gui_slicebatoms,
    gui_toolbar,
    gui_cell,
    gui_slicebonds,
    gui_plane,
    gui_render,
    gui_pymatgen,
    gui_pubchem,
    gui_rscb,
    ui_list_species,
    view3d_mt_batoms_add,
    view3d_mt_object_context_menu,
    view3d_mt_edit_mesh_context_menu,
)

class BatomsCollection(bpy.types.PropertyGroup):
    """
    Collection properties of all panel properties.
    """
    batoms: PointerProperty(type=gui_batoms.BatomsProperties)
    batom: PointerProperty(type=gui_slicebatoms.BatomProperties)
    cell: PointerProperty(type=gui_cell.CellProperties)
    bond: PointerProperty(type=gui_slicebonds.BondProperties)
    plane: PointerProperty(type=gui_plane.PlaneProperties)
    render: PointerProperty(type=gui_render.RenderProperties)
    pymatgen: PointerProperty(type=gui_pymatgen.PymatgenProperties)
    pubchem: PointerProperty(type=gui_pubchem.PubchemProperties)
    rscb: PointerProperty(type=gui_rscb.RSCBProperties)


classes = [
    gui_batoms.Batoms_PT_prepare,
    gui_batoms.BatomsProperties,
    gui_slicebatoms.Batom_PT_prepare,
    gui_slicebatoms.BatomProperties,
    gui_slicebonds.Bond_PT_prepare,
    gui_slicebonds.BondProperties,
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
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_batoms_model_style,
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_batoms_radius_style,
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_batoms_color_style,
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_batoms_label,
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_batoms,
    view3d_mt_object_context_menu.VIEW3D_MT_object_context_bonds,
    view3d_mt_edit_mesh_context_menu.VIEW3D_MT_edit_mesh_context_batoms_model_style,
    view3d_mt_edit_mesh_context_menu.VIEW3D_MT_edit_mesh_context_batoms,
    ui_list_species.BATOMS_MT_species_context_menu,
    ui_list_species.BATOMS_UL_species,
    ui_list_species.BATOMS_PT_species,
    BatomsCollection,
]



def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    scene = bpy.types.Scene
    scene.batoms = PointerProperty(type=BatomsCollection)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    scene = bpy.types.Scene
    del scene.batoms


def register_menu():
    bpy.types.VIEW3D_MT_add.prepend(view3d_mt_batoms_add.menu_func)
    bpy.types.VIEW3D_MT_object_context_menu.prepend(view3d_mt_object_context_menu.menu_func)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.prepend(view3d_mt_edit_mesh_context_menu.menu_func)

def unregister_menu():
    bpy.types.VIEW3D_MT_add.remove(view3d_mt_batoms_add.menu_func)
    bpy.types.VIEW3D_MT_object_context_menu.remove(view3d_mt_object_context_menu.menu_func)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.remove(view3d_mt_edit_mesh_context_menu.menu_func)
    


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

