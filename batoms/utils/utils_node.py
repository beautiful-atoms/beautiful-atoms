import bpy
import numpy as np
# TODO: 4.2+ support
from .. import __package__ as batoms


def get_node_by_name(nodes, name, type=None):
    node = nodes.get(name)
    if node is None:
        node = nodes.new(type)
        node.name = name
    return node


def create_node_tree(name, node_group_type="GeometryNodeTree", interface=[]):
    """Create a node_tree by name and type.
    Add default interface"""
    # create a new node_tree if node_tree is None
    node_tree = bpy.data.node_groups.get(name)
    if node_tree is None:
        node_tree = bpy.data.node_groups.new(name, type=node_group_type)
        # add default interface if not exist
        for i in interface:
            if bpy.app.version_string >= "4.0.0":
                node_tree.interface.new_socket(name=i[0], socket_type=i[1], in_out=i[2])
            else:
                if i[2] == "INPUT":
                    node_tree.inputs.new(i[1], i[0])
                elif i[2] == "OUTPUT":
                    node_tree.outputs.new(i[1], i[0])
        # create input and output node
        node_tree.nodes.new("NodeGroupInput")
        node_tree.nodes.new("NodeGroupOutput")
    return node_tree


def get_socket_by_identifier(node, identifier, type="inputs"):
    """Get sockets by identifier"""
    for inp in getattr(node, type):
        if inp.identifier == identifier:
            return inp
    return None


def get_node_with_node_tree_by_name(
    nodes,
    name,
    node_type="GeometryNodeGroup",
    node_group_type="GeometryNodeTree",
    interface=[],
):
    """Get the node with a node_tree by name.
    If None, Create a new one based on type.
    Return:
        node (bpy.types.Node): A node with a node_tree
    """
    node = nodes.get(name)
    # create a new node if node is None
    if node is None:
        node = nodes.new(type=node_type)
        node.name = name
    node_tree = create_node_tree(name, node_group_type, interface)
    node.node_tree = node_tree
    return node


def get_cell_node(parent_tree, label):
    """Get the position of the cell."""
    from batoms.utils.utils_node import (
        get_node_by_name,
        create_node_tree,
    )

    default_interface = [
        ["Cell", "NodeSocketObject", "INPUT"],
        ["Matrix", "NodeSocketMatrix", "OUTPUT"],
        ["Transpose", "NodeSocketMatrix", "OUTPUT"],
        ["Invert", "NodeSocketMatrix", "OUTPUT"],
    ]
    name = "Cell_Array_%s" % (label)
    node = get_node_by_name(parent_tree.nodes, name=name, type="GeometryNodeGroup")
    node_tree = create_node_tree(name=name, interface=default_interface)
    node.node_tree = node_tree
    nodes = node_tree.nodes
    links = node_tree.links
    GroupInput = nodes[0]
    GroupOutput = nodes[1]
    # link the input to parent node
    # ------------------------------------------------------------------
    CellObject = get_node_by_name(
        nodes, f"{label}_CellObject", "GeometryNodeObjectInfo"
    )
    links.new(GroupInput.outputs["Cell"], CellObject.inputs["Object"])
    Position = get_node_by_name(
        nodes, "%s_Position" % (label), "GeometryNodeInputPosition"
    )
    combined_matrix = get_node_by_name(
        nodes,
        "%s_CombineMatrix" % (label),
        "FunctionNodeCombineMatrix",
    )
    transpose_matrix = get_node_by_name(
        nodes,
        "%s_TransposeMatrix" % (label),
        "FunctionNodeTransposeMatrix",
    )
    invert_matrix = get_node_by_name(
        nodes,
        "%s_InvertMatrix" % (label),
        "FunctionNodeInvertMatrix",
    )
    links.new(combined_matrix.outputs[0], transpose_matrix.inputs[0])
    links.new(transpose_matrix.outputs[0], invert_matrix.inputs[0])
    links.new(transpose_matrix.outputs[0], GroupOutput.inputs["Transpose"])
    links.new(invert_matrix.outputs[0], GroupOutput.inputs["Invert"])

    for i in range(3):
        PositionAtIndex = get_node_by_name(
            nodes, f"{label}_PositionAtIndex_{i}", "GeometryNodeSampleIndex"
        )
        PositionAtIndex.data_type = "FLOAT_VECTOR"
        PositionAtIndex.inputs["Index"].default_value = i + 1
        SeparateXYZ = get_node_by_name(
            nodes, f"{label}_SeparateXYZ_{i}", "ShaderNodeSeparateXYZ"
        )
        links.new(CellObject.outputs["Geometry"], PositionAtIndex.inputs[0])
        links.new(Position.outputs["Position"], PositionAtIndex.inputs["Value"])
        links.new(PositionAtIndex.outputs["Value"], SeparateXYZ.inputs[0])
        for j in range(3):
            links.new(
                SeparateXYZ.outputs[j],
                combined_matrix.inputs[i * 4 + j],
            )
        links.new(combined_matrix.outputs["Matrix"], GroupOutput.inputs["Matrix"])

    return node


def get_projected_position(
    node_tree: bpy.types.NodeTree,
    position_socket: bpy.types.NodeSocket,
    cell_obj: bpy.types.Object,
    label: str,
    scaled: bool = True,
) -> bpy.types.Node:
    """Get the projected position of the atom.
    If scaled is True, the position is scaled by the cell.
    If scaled is False, the scaled position is transformed to the cartesian position.
    """
    cell_node = get_cell_node(node_tree, label)
    cell_node.inputs["Cell"].default_value = cell_obj
    name = "%s_ProjectPoint_%s" % (label, "scaled" if scaled else "cartesian")
    project_point = get_node_by_name(
        node_tree.nodes,
        name,
        "FunctionNodeProjectPoint",
    )
    node_tree.links.new(position_socket, project_point.inputs[0])
    key = "Invert" if scaled else "Matrix"
    node_tree.links.new(cell_node.outputs[key], project_point.inputs[1])
    return project_point


def vectorDotMatrix(
    node_tree: bpy.types.NodeTree,
    vector_output: bpy.types.NodeSocket,
    matrix: np.ndarray,
    label: str,
    name: str,
) -> bpy.types.Node:
    """Create a vector dot matrix operation."""
    project_point = get_node_by_name(
        node_tree.nodes,
        "%s_ProjectPoint_%s" % (label, name),
        "FunctionNodeProjectPoint",
    )
    combined_matrix = get_node_by_name(
        node_tree.nodes,
        "%s_CombineMatrix_%s" % (label, name),
        "FunctionNodeCombineMatrix",
    )
    for i in range(3):
        for j in range(3):
            index = i * 4 + j
            combined_matrix.inputs[index].default_value = matrix[i, j]
    node_tree.links.new(vector_output, project_point.inputs[0])
    node_tree.links.new(combined_matrix.outputs[0], project_point.inputs[1])
    return project_point
