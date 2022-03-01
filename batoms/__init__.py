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
    gui_volume,
    gui_cell,
    gui_bond,
    gui_polyhedra,
    gui_ase,
    gui_plane,
    gui_render,
    gui_pymatgen,
    gui_pubchem,
    gui_rscb,
    gui_tool,
    view3d_mt_batoms_add,
)

from .ops import (build_object, 
    build_surface, 
    build_nanotube,
    build_nanoparticle,
    build_edit_molecule, 
    transform
    )

from .modal import (record_selection,
            rigid_body,
            force_field,
            replace_atoms,
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
    build_object.AddMolecule,
    build_object.AddBulk,
    build_object.AddAtoms,
    build_object.AddSurface,
    build_object.AddRootSurface,
    build_object.deleteBatoms,
    build_surface.BuildSurfaceFCC100,
    build_surface.BuildSurfaceFCC110,
    build_surface.BuildSurfaceFCC111,
    build_surface.BuildSurfaceFCC111Root,
    build_surface.BuildSurfaceBCC100,
    build_surface.BuildSurfaceBCC110,
    build_surface.BuildSurfaceBCC111,
    build_surface.BuildSurfaceBCC111Root,
    build_surface.BuildSurfaceHCP0001,
    build_surface.BuildSurfaceHCP0001Root,
    build_nanotube.BuildNanotube,
    build_nanoparticle.BuildDecahedron,
    build_nanoparticle.BuildIcosahedron,
    build_nanoparticle.BuildOctahedron,
    build_edit_molecule.MolecueReplaceElement,
    replace_atoms.EDIT_MESH_OT_replace_atoms,
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
    gui_ase.ASE_PT_prepare,
    gui_ase.ASEProperties,
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


def register():

    # class
    for cls in classes:
        bpy.utils.register_class(cls)

    # menu
    bpy.types.TOPBAR_MT_file_import.append(gui_io.menu_func_import_batoms)
    bpy.types.VIEW3D_MT_add.prepend(view3d_mt_batoms_add.menu_func)

    # tool
    bpy.utils.register_tool(gui_tool.AddMolecule, after={"builtin.scale_cage"}, separator=True, group=True)
    bpy.utils.register_tool(gui_tool.AddBulk, after={gui_tool.AddMolecule.bl_idname})
    bpy.utils.register_tool(gui_tool.AddSurface, after={gui_tool.AddBulk.bl_idname})
    bpy.utils.register_tool(gui_tool.AddAtoms, after={gui_tool.AddSurface.bl_idname})
    bpy.utils.register_tool(gui_tool.MoleculeReplaceElement, after={"builtin.scale_cage"}, separator=True, group=True)


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
    scene.asepanel = PointerProperty(type=gui_ase.ASEProperties)
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
    bpy.utils.unregister_tool(gui_tool.AddMolecule)
    bpy.utils.unregister_tool(gui_tool.AddBulk)
    bpy.utils.unregister_tool(gui_tool.AddSurface)
    bpy.utils.unregister_tool(gui_tool.AddAtoms)
    bpy.utils.unregister_tool(gui_tool.MoleculeReplaceElement)


if __name__ == "__main__":

    register()
