"""
Adds an object mode tool to the toolbar.
# TODO Not working in  background mode,
fialed to register the keymaps.
"""

from bpy.types import WorkSpaceTool


class BatomsCell(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "OBJECT"  # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms.apply_cell"
    bl_label = "Cell"
    bl_description = (
        "Apply new cell parameters."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "ops.transform.resize"
    bl_widget = None

    bl_keymap = (
        ("batoms.apply_cell", {"type": "RIGHTMOUSE", "value": "PRESS"},
         {"properties": []}),
        # ("batoms.apply_cell", {"type": "LEFTMOUSE", "value": "PRESS", "ctrl": True},
        #  {"properties": []}),
    )

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.apply_cell")
        layout.prop(props, "a")
        layout.prop(props, "b")
        layout.prop(props, "c")


class BatomsTransform(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "OBJECT"  # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms.apply_transform"
    bl_label = "Transform"
    bl_description = (
        "Apply new transform parameters."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "ops.transform.shear"
    bl_widget = None

    bl_keymap = (
        ("batoms.apply_transform", {"type": "RIGHTMOUSE", "value": "PRESS"},
         {"properties": []}),
        # ("batoms.apply_transform", {"type": "LEFTMOUSE", "value": "PRESS", "ctrl": True},
        #  {"properties": []}),
    )

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.apply_transform")
        layout.prop(props, "a")
        layout.prop(props, "b")
        layout.prop(props, "c")


class BatomsBoundary(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "OBJECT"  # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms.apply_boundary"
    bl_label = "Boundary"
    bl_description = (
        "Apply new Boundary parameters."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "ops.transform.shrink_fatten"
    bl_widget = None

    bl_keymap = (
        ("batoms.apply_boundary", {"type": "RIGHTMOUSE", "value": "PRESS"},
         {"properties": []}),
        # ("batoms.apply_boundary", {"type": "LEFTMOUSE", "value": "PRESS", "ctrl": True},
        #  {"properties": []}),
    )

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.apply_boundary")
        layout.prop(props, "a")
        layout.prop(props, "b")
        layout.prop(props, "c")


class AddSurface(WorkSpaceTool):
    bl_space_type = 'VIEW_3D'
    bl_context_mode = 'OBJECT'

    bl_idname = "batoms.surface_add"
    bl_label = "Add surface"
    bl_description = (
        "Create an surface structure from a bulk."
    )
    bl_icon = "ops.transform.translate"

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.surface_add")
        layout.prop(props, "label")
        layout.prop(props, "indices")
        layout.prop(props, "size")
        layout.prop(props, "vacuum")
        layout.prop(props, "nlayer")
        layout.prop(props, "termination")


class MoleculeEditElement(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "EDIT_MESH"  # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms_tool.molecule_edit_atom"
    bl_label = "Edit atoms"
    bl_description = (
        "Create an molecule structure from a database."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "batoms_molecule_edit_atom"
    bl_widget = None

    bl_keymap = (
        ("batoms.molecule_edit_atom", {"type": "RIGHTMOUSE", "value": "PRESS"},
         {"properties": []}),
        ("batoms.molecule_edit_atom", {"type": "LEFTMOUSE", "value": "PRESS", "ctrl": True},
         {"properties": []}),
    )

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.molecule_edit_atom")
        layout.prop(props, "element")
        layout.prop(props, "bond_order")


class MolecueEditBond(WorkSpaceTool):
    bl_space_type = "VIEW_3D"
    bl_context_mode = "EDIT_MESH"  # "OBJECT" or "EDIT_MESH"

    bl_idname = "batoms_tool.molecule_edit_bond"
    bl_label = "Edit bond"
    bl_description = (
        "Create an molecule structure from a database."
    )
    # bl_icon = "ops.generic.select_circle"
    bl_icon = "batoms_molecule_edit_atom"
    bl_widget = None

    def draw_settings(context, layout, tool):
        props = tool.operator_properties("batoms.molecule_edit_bond")
        layout.prop(props, "order")
        layout.prop(props, "rotate")
        layout.prop(props, "hydrogen")
