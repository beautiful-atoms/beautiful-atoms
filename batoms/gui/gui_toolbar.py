"""
Adds an object mode tool to the toolbar.
"""

from bpy.types import WorkSpaceTool
import os


class AddAtoms(WorkSpaceTool):
    bl_space_type = 'VIEW_3D'
    bl_context_mode = 'OBJECT'

    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = (
        "Create an element structure."
    )
    # tools require the icon to be a binary dat file.
    # which can be generated using the python script: blender_icons_geom.py
    bl_icon = "ops.generic.select_circle"
    # bl_icon = os.path.join("../icons" , "add_atoms")
    bl_widget = None
    
    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.add_atoms")
        layout.prop(props, "formula")
        layout.prop(props, "label")


class AddMolecule(WorkSpaceTool):
    bl_space_type = 'VIEW_3D'
    bl_context_mode = 'OBJECT'

    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = (
        "Create an molecule structure from a database."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "batoms_add_molecule"
    bl_widget = None
    
    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.add_molecule")
        layout.prop(props, "formula")
        layout.prop(props, "label")


class AddBulk(WorkSpaceTool):
    bl_space_type = 'VIEW_3D'
    bl_context_mode = 'OBJECT'

    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = (
        "Create an bulk structure from a database."
    )
    bl_icon = "ops.generic.select_lasso"
    bl_widget = None

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.add_bulk")
        layout.prop(props, "formula")
        layout.prop(props, "label")


class AddSurface(WorkSpaceTool):
    bl_space_type = 'VIEW_3D'
    bl_context_mode = 'OBJECT'

    bl_idname = "batoms.add_surface"
    bl_label = "Add surface"
    bl_description = (
        "Create an surface structure from a bulk."
    )
    bl_icon = "ops.transform.translate"

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.add_surface")
        layout.prop(props, "label")
        layout.prop(props, "indices")
        layout.prop(props, "size")
        layout.prop(props, "vacuum")
        layout.prop(props, "nlayer")
        layout.prop(props, "termination")



class MoleculeEditElement(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "EDIT_MESH" # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms.molecule_edit_atom"
    bl_label = "Edit atoms"
    bl_description = (
        "Create an molecule structure from a database."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "batoms_add_molecule"
    bl_widget = None
    

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.molecule_edit_atom")
        layout.prop(props, "element")
        layout.prop(props, "bond_order")


class MolecueEditBond(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "EDIT_MESH" # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms.molecule_edit_bond"
    bl_label = "Edit bond"
    bl_description = (
        "Create an molecule structure from a database."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "batoms_add_molecule"
    bl_widget = None
    

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.molecule_edit_bond")
        layout.prop(props, "order")
        layout.prop(props, "rotate")
        layout.prop(props, "hydrogen")

