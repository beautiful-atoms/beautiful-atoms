import bpy
from bpy.props import PointerProperty

from . import (
    gui_batoms,
    gui_io,
    gui_slicebatoms,
    gui_toolbar,
    gui_cell,
    gui_plane,
    gui_render,
    gui_volumetric_data,
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
    plane: PointerProperty(type=gui_plane.PlaneProperties)
    render: PointerProperty(type=gui_render.RenderProperties)
    io: PointerProperty(type=gui_io.BatomsPropertiesIO)
    volumetric_data: PointerProperty(type=gui_volumetric_data.VolumetricDataProperties)


classes = [
    gui_batoms.Batoms_PT_prepare,
    gui_batoms.BatomsProperties,
    gui_slicebatoms.Batom_PT_prepare,
    gui_slicebatoms.BatomProperties,
    gui_cell.Cell_PT_prepare,
    gui_cell.CellProperties,
    gui_plane.Plane_PT_prepare,
    gui_plane.PlaneProperties,
    gui_plane.AddButton,
    gui_render.Render_PT_prepare,
    gui_render.RenderProperties,
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
    gui_volumetric_data.VolumetricDataProperties,
    gui_volumetric_data.VIEW3D_PT_Batoms_volumetric_data,
    gui_volumetric_data.BATOMS_MT_volumetric_data_context_menu,
    gui_volumetric_data.BATOMS_UL_volumetric_data,
    gui_volumetric_data.BATOMS_PT_volumetric_data,
    gui_io.BatomsPropertiesIO,
    gui_io.VIEW3D_PT_Batoms_io,
    gui_io.VIEW3D_PT_Batoms_io_materials_project,
    gui_io.Batoms_IO_search_material_project,
    gui_io.VIEW3D_PT_Batoms_io_pubchem,
    gui_io.Batoms_IO_search_pubchem,
    gui_io.VIEW3D_PT_Batoms_io_RSCB,
    gui_io.Batoms_IO_search_RSCB,
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
