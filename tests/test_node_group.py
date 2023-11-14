import bpy


def test_create_node_group():
    from batoms.utils.butils import get_node_with_node_tree_by_name

    parent = bpy.data.node_groups.new("node_group1", "GeometryNodeTree")
    node = get_node_with_node_tree_by_name(parent.nodes, "node_group2")
    node.node_tree.nodes.new("GeometryNodeInputNamedAttribute")
    assert len(parent.nodes) == 1
    # default input and output nodes
    assert len(node.node_tree.nodes) == 3
