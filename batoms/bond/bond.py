"""Definition of the Bond class.

This module defines the Bond object in the Batoms package.

"""

import bpy
from time import time
from batoms.utils.butils import object_mode, get_node_by_name
from batoms.utils import string2Number, number2String
import numpy as np
from batoms.base.object import ObjectGN
from batoms.base.collection import BaseCollection
from .setting import BondSettings
from batoms.bond.search_bond import SearchBond, default_search_bond_datas

import logging

logger = logging.getLogger(__name__)


default_bond_attributes = [
    {"name": "atoms_index0", "data_type": "INT", "domain": "EDGE"},
    {"name": "atoms_index1", "data_type": "INT", "domain": "EDGE"},
    {"name": "atoms_index2", "data_type": "INT", "domain": "EDGE"},
    {"name": "atoms_index3", "data_type": "INT", "domain": "EDGE"},
    {"name": "species_index0", "data_type": "INT", "domain": "EDGE"},
    {"name": "species_index1", "data_type": "INT", "domain": "EDGE"},
    {"name": "bond_offset0", "data_type": "FLOAT_VECTOR", "domain": "EDGE"},
    {"name": "bond_offset1", "data_type": "FLOAT_VECTOR", "domain": "EDGE"},
    {"name": "bond_offset2", "data_type": "FLOAT_VECTOR", "domain": "EDGE"},
    {"name": "bond_offset3", "data_type": "FLOAT_VECTOR", "domain": "EDGE"},
    {"name": "bond_order", "data_type": "INT", "domain": "EDGE"},
    {"name": "bond_style", "data_type": "INT", "domain": "EDGE"},
    {"name": "bond_show", "data_type": "INT", "domain": "EDGE"},
    {"name": "bond_model_style", "data_type": "INT", "domain": "EDGE"},
    {"name": "polyhedra", "data_type": "INT", "domain": "EDGE"},
    {"name": "second_bond", "data_type": "INT", "domain": "EDGE"},
]

default_bond_datas = {
    "atoms_index0": np.ones(0, dtype=int),
    "atoms_index1": np.ones(0, dtype=int),
    "atoms_index2": np.ones(0, dtype=int),
    "atoms_index3": np.ones(0, dtype=int),
    "species_index0": np.ones(0, dtype=int),
    "species_index1": np.ones(0, dtype=int),
    "centers": np.zeros((1, 0, 3)),
    # 'vectors':np.zeros((0, 3)),
    "bond_offset0": np.zeros((0, 3)),
    "bond_offset1": np.zeros((0, 3)),
    "bond_offset2": np.zeros((0, 3)),
    "bond_offset3": np.zeros((0, 3)),
    # 'eulers':np.eye(3),
    # 'lengths':np.zeros((0, 3)),
    # 'widths': np.ones(0, dtype=float),
    "bond_show": np.zeros(0, dtype=int),
    "bond_order": np.zeros(0, dtype=int),
    "bond_style": np.zeros(0, dtype=int),
    "bond_model_style": np.ones(0, dtype=int),
    "polyhedra": np.ones(0, dtype=float),
    "second_bond": np.ones(0, dtype=int),
}


class Bond(BaseCollection, ObjectGN):
    """Bbond Class

    A Bbond object is linked to this main collection in Blender.

    Parameters:


    """

    def __init__(
        self,
        label=None,
        batoms=None,
        bond_datas=None,
        location=np.array([0, 0, 0]),
    ):
        #
        self.batoms = batoms
        self.label = label
        ObjectGN.__init__(self, label)
        BaseCollection.__init__(self, coll_name=label)
        self.settings = BondSettings(self.label, batoms=batoms, bonds=self)
        self.build_geometry_node()
        self._search_bond = SearchBond(self.label, batoms=batoms, load=True)

    @property
    def bond_node(self):
        """Get the top level node of bond node group."""
        from batoms.utils.butils import get_node_with_node_tree_by_name

        default_interface = [
            ["Geometry", "NodeSocketGeometry", "INPUT"],
            ["Geometry", "NodeSocketGeometry", "OUTPUT"],
        ]
        name = f"Bond_{self.label}"
        node = get_node_with_node_tree_by_name(
            self.gn_node_group.nodes, name=name, interface=default_interface
        )
        return node

    def get_bond_pair_node(self, name):
        """Get the node of bond pair node group.
        Create the node if not exist."""
        from batoms.utils.butils import get_node_with_node_tree_by_name

        # create group input and output sockets
        default_interface = [
            ["Geometry", "NodeSocketGeometry", "INPUT"],
            ["Species_index0", "NodeSocketInt", "INPUT"],
            ["Species0", "NodeSocketInt", "INPUT"],
            ["Species_index1", "NodeSocketInt", "INPUT"],
            ["Species1", "NodeSocketInt", "INPUT"],
            ["Model_style0", "NodeSocketInt", "INPUT"],
            ["Model_style1", "NodeSocketInt", "INPUT"],
            ["Order", "NodeSocketInt", "INPUT"],
            ["Style", "NodeSocketInt", "INPUT"],
            ["Length", "NodeSocketFloat", "INPUT"],
            ["Rotation", "NodeSocketVector", "INPUT"],
            ["Geometry", "NodeSocketGeometry", "OUTPUT"],
        ]
        # Create Node Group if not exit
        node = get_node_with_node_tree_by_name(
            self.bond_node.node_tree.nodes, name=name, interface=default_interface
        )
        return node

    def build_geometry_node(self):
        """
        v = (p1 + o1) - (p2 + o2)
        len(v) < max_length
        v align euler
        """
        from batoms.utils.butils import get_node_by_name, get_socket_by_identifier

        tstart = time()
        parent = self.gn_node_group
        node = self.bond_node
        nodes = node.node_tree.nodes
        links = node.node_tree.links
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        # link the input to parent node
        JoinGeometry = get_node_by_name(
            parent.nodes,
            "%s_JoinGeometry" % self.batoms.label,
            "GeometryNodeJoinGeometry",
        )
        parent.links.new(
            parent.nodes["Group Input"].outputs["Geometry"], node.inputs["Geometry"]
        )
        # link the outputs to parent node
        parent.links.new(node.outputs["Geometry"], JoinGeometry.inputs["Geometry"])
        # ------------------------------------------------------------------
        # this join geometry is belong to this node group
        JoinGeometry = get_node_by_name(
            nodes, "%s_JoinGeometry" % self.label, "GeometryNodeJoinGeometry"
        )
        links.new(JoinGeometry.outputs["Geometry"], GroupOutput.inputs["Geometry"])
        # ------------------------------------------------------------------
        self.build_bond_mesh()
        # get the species by atoms_index
        # order
        SpeciesIndexAttribute = get_node_by_name(
            nodes,
            "%s_NamedAttribute_species_index" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        SpeciesIndexAttribute.inputs["Name"].default_value = "species_index"
        SpeciesIndexAttribute.data_type = "INT"
        SpeciesAtIndexs = []
        for i in range(2):
            AtomsIndexAttribute = get_node_by_name(
                nodes,
                "%s_NamedAttribute_atoms_index%s" % (self.label, i),
                "GeometryNodeInputNamedAttribute",
            )
            AtomsIndexAttribute.inputs["Name"].default_value = f"atoms_index{i}"
            AtomsIndexAttribute.data_type = "INT"
            SpeciesAtIndex = get_node_by_name(
                nodes,
                "%s_SpeciesAtIndex%s" % (self.label, i),
                "GeometryNodeSampleIndex",
            )
            SpeciesAtIndex.data_type = "INT"
            SpeciesAtIndexs.append(SpeciesAtIndex)
            links.new(GroupInput.outputs["Geometry"], SpeciesAtIndex.inputs["Geometry"])
            output_socket = get_socket_by_identifier(
                AtomsIndexAttribute, "Attribute_Int", type="outputs"
            )
            links.new(output_socket, SpeciesAtIndex.inputs["Index"])
            input_socket = get_socket_by_identifier(SpeciesAtIndex, "Value_Int")
            output_socket = get_socket_by_identifier(
                SpeciesIndexAttribute, "Attribute_Int", type="outputs"
            )
            links.new(output_socket, input_socket)
        logger.debug("Build geometry nodes for bonds: %s" % (time() - tstart))
        #

    def build_bond_mesh(self):
        """Create geometry node for Bond Mesh.
        1) Transform the edge to mesh points.
        2) Set the position of the points to the center of the bond.
        3) Calculate the bond vector and length.
        """
        from batoms.utils.butils import (
            get_socket_by_identifier,
            get_node_by_name,
            create_node_tree,
        )

        default_interface = [
            ["Geometry", "NodeSocketGeometry", "INPUT"],
            ["Geometry", "NodeSocketGeometry", "OUTPUT"],
            ["Length", "NodeSocketFloat", "OUTPUT"],
            ["Rotation", "NodeSocketVector", "OUTPUT"],
        ]
        parent_tree = self.bond_node.node_tree
        name = "Bond_Mesh_%s" % (self.label)
        node = get_node_by_name(parent_tree.nodes, name=name, type="GeometryNodeGroup")
        node_tree = create_node_tree(name=name, interface=default_interface)
        node.node_tree = node_tree
        nodes = node_tree.nodes
        links = node_tree.links
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        # link the input to parent node
        parent_tree.links.new(
            parent_tree.nodes[0].outputs["Geometry"], node.inputs["Geometry"]
        )
        # -------------------------------------------------------------
        # calculate bond vector, length, rotation based on the index
        # Get four positions from batoms
        # atoms_index0, atoms_index1 of bond
        # atoms_index2, atoms_index3 of the second bond, for high order bond plane
        PositionBatoms = get_node_by_name(
            nodes, "%s_PositionBatoms" % (self.label), "GeometryNodeInputPosition"
        )
        PositionsAtIndexs = []
        AtomsIndexAttributes = []
        for i in range(4):
            tmp = get_node_by_name(
                nodes,
                "%s_PositionsAtIndex%s" % (self.label, i),
                "GeometryNodeSampleIndex",
            )
            tmp.data_type = "FLOAT_VECTOR"
            links.new(GroupInput.outputs["Geometry"], tmp.inputs["Geometry"])
            PositionsAtIndexs.append(tmp)
            tmp = get_node_by_name(
                nodes,
                "%s_NamedAttribute_atoms_index%s" % (self.label, i),
                "GeometryNodeInputNamedAttribute",
            )
            tmp.inputs["Name"].default_value = f"atoms_index{i}"
            tmp.data_type = "INT"
            AtomsIndexAttributes.append(tmp)
        for i in range(4):
            socket = get_socket_by_identifier(PositionsAtIndexs[i], "Value_Vector")
            links.new(PositionBatoms.outputs[0], socket)
            socket = get_socket_by_identifier(
                AtomsIndexAttributes[i], "Attribute_Int", type="outputs"
            )
            links.new(socket, PositionsAtIndexs[i].inputs["Index"])
        # ------------------------------------------------------------------
        # add positions with offsets
        # first step: get offsets from named attribute
        OffsetsAttributes = []
        for i in range(4):
            tmp = get_node_by_name(
                nodes,
                "%s_NamedAttribute_offsets%s" % (self.label, i),
                "GeometryNodeInputNamedAttribute",
            )
            tmp.inputs["Name"].default_value = f"bond_offset{i}"
            tmp.data_type = "FLOAT_VECTOR"
            OffsetsAttributes.append(tmp)
        # second step: we need five add operations
        # four: Get the positions with offset for four atoms
        # one: Get center = (positions1 + positions2)/2
        VectorAdd = []
        for i in range(5):
            tmp = get_node_by_name(
                nodes, "%s_VectorAdd%s" % (self.label, i), "ShaderNodeVectorMath"
            )
            tmp.operation = "ADD"
            VectorAdd.append(tmp)
        for i in range(4):
            socket = get_socket_by_identifier(
                PositionsAtIndexs[i], "Value_Vector", type="outputs"
            )
            links.new(socket, VectorAdd[i].inputs[0])
            links.new(OffsetsAttributes[i].outputs[0], VectorAdd[i].inputs[1])
        # divide by 2 to get the center
        VectorDivide = get_node_by_name(
            nodes, "VectorDivide_%s" % self.label, "ShaderNodeVectorMath"
        )
        VectorDivide.operation = "DIVIDE"
        VectorDivide.inputs[1].default_value = (2, 2, 2)
        links.new(VectorAdd[0].outputs[0], VectorAdd[4].inputs[0])
        links.new(VectorAdd[1].outputs[0], VectorAdd[4].inputs[1])
        links.new(VectorAdd[4].outputs[0], VectorDivide.inputs[0])
        # set center of the bond
        SetPosition = get_node_by_name(
            nodes, "%s_SetPosition" % self.label, "GeometryNodeSetPosition"
        )
        MeshToPoints = get_node_by_name(
            nodes, "%s_MeshToPoints_edge" % (self.label), "GeometryNodeMeshToPoints"
        )
        MeshToPoints.mode = "EDGES"
        links.new(GroupInput.outputs["Geometry"], MeshToPoints.inputs["Mesh"])
        links.new(MeshToPoints.outputs["Points"], SetPosition.inputs["Geometry"])
        links.new(VectorDivide.outputs[0], SetPosition.inputs["Position"])
        links.new(SetPosition.outputs["Geometry"], GroupOutput.inputs["Geometry"])
        # get the vector for the bond and the length
        # also the vector for the second bond
        VectorSubtract = []
        for i in range(2):
            tmp = get_node_by_name(
                nodes, "%s_VectorSubtract%s" % (self.label, i), "ShaderNodeVectorMath"
            )
            tmp.operation = "SUBTRACT"
            VectorSubtract.append(tmp)
        VectorLength = get_node_by_name(
            nodes, "%s_VectorLength" % self.label, "ShaderNodeVectorMath"
        )
        VectorLength.operation = "LENGTH"
        VectorCross0 = get_node_by_name(
            nodes, "%s_VectorCross0" % self.label, "ShaderNodeVectorMath"
        )
        VectorCross0.operation = "CROSS_PRODUCT"
        links.new(VectorAdd[0].outputs[0], VectorSubtract[0].inputs[0])
        links.new(VectorAdd[1].outputs[0], VectorSubtract[0].inputs[1])
        links.new(VectorAdd[2].outputs[0], VectorSubtract[1].inputs[0])
        links.new(VectorAdd[3].outputs[0], VectorSubtract[1].inputs[1])
        # calc the bond length, use it to build scale
        links.new(VectorSubtract[0].outputs[0], VectorLength.inputs[0])
        # cross for rotation, for high order bond
        links.new(VectorSubtract[0].outputs[0], VectorCross0.inputs[0])
        links.new(VectorSubtract[1].outputs[0], VectorCross0.inputs[1])
        # get Euler for rotation
        # we need align two vectors to fix a plane
        # we build the instancer by fix bond diection to Z, and
        # high order bond shift to X
        # thus the the normal of high order bond plane is Y
        AlignEuler = []
        for i in range(2):
            tmp = get_node_by_name(
                nodes,
                "%s_AlignEuler%s" % (self.label, i),
                "FunctionNodeAlignEulerToVector",
            )
            AlignEuler.append(tmp)
        AlignEuler[0].axis = "Z"
        AlignEuler[1].axis = "Y"
        # We should fix Z when align Y
        AlignEuler[1].pivot_axis = "Z"
        links.new(VectorSubtract[0].outputs[0], AlignEuler[0].inputs["Vector"])
        links.new(AlignEuler[0].outputs[0], AlignEuler[1].inputs["Rotation"])
        links.new(VectorCross0.outputs[0], AlignEuler[1].inputs["Vector"])
        links.new(VectorLength.outputs["Value"], GroupOutput.inputs["Length"])
        links.new(AlignEuler[1].outputs[0], GroupOutput.inputs["Rotation"])

    def update_geometry_nodes(self):
        tstart = time()
        # for sp in self.settings:
        # self.add_bond_pair_node(sp.as_dict())
        # only build pair in bondlists
        bondlists = self.arrays
        if len(bondlists["atoms_index0"]) == 0:
            return
        pairs = np.concatenate(
            (
                bondlists["species_index0"].reshape(-1, 1),
                bondlists["species_index1"].reshape(-1, 1),
                bondlists["bond_order"].reshape(-1, 1),
                bondlists["bond_style"].reshape(-1, 1),
            ),
            axis=1,
        )
        pairs = np.unique(pairs, axis=0)
        pairs = pairs.reshape(-1, 4)
        for pair in pairs:
            sp1 = number2String(pair[0])
            sp2 = number2String(pair[1])
            order = pair[2]
            style = pair[3]
            sp = self.settings["%s-%s" % (sp1, sp2)]
            self.add_bond_pair_node(sp.as_dict(), order, style)

        logger.debug("Update geometry nodes for bonds: %s" % (time() - tstart))

    def add_bond_pair_node(self, sp, order=None, style=None):
        """Add geometry node for the bond pair."""
        from batoms.utils.butils import get_node_by_name, get_socket_by_identifier

        parent_tree = self.bond_node.node_tree
        # tstart = time()
        if not order:
            order = sp["order"]
        if not style:
            style = int(sp["style"])
        name = "Bond_%s_%s_%s_%s" % (self.label, sp["name"], order, style)
        node = self.get_bond_pair_node(name)
        nodes = node.node_tree.nodes
        links = node.node_tree.links
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        # link the input to parent node
        BondMesh = get_node_by_name(
            parent_tree.nodes, "Bond_Mesh_%s" % (self.label), "GeometryNodeSetPosition"
        )
        CompareSpecies = []
        for i in range(2):
            SpeciesAtIndex = get_node_by_name(
                parent_tree.nodes,
                "%s_SpeciesAtIndex%s" % (self.label, i),
                "GeometryNodeSampleIndex",
            )
            output_socket = get_socket_by_identifier(
                SpeciesAtIndex, "Value_Int", type="outputs"
            )
            parent_tree.links.new(output_socket, node.inputs[f"Species_index{i}"])
            # e.g C-H bond, species0 is C, species1 is H
            tmp = get_node_by_name(
                nodes,
                "%s_CompareSpecies_%s_%s" % (self.label, sp["species1"], i),
                "FunctionNodeCompare",
            )
            tmp.operation = "EQUAL"
            tmp.data_type = "INT"
            node.inputs[f"Species{i}"].default_value = string2Number(
                sp[f"species{i + 1}"]
            )
            CompareSpecies.append(tmp)
            A_socket = get_socket_by_identifier(tmp, "A_INT")
            B_socket = get_socket_by_identifier(tmp, "B_INT")
            links.new(GroupInput.outputs[f"Species_index{i}"], A_socket)
            links.new(GroupInput.outputs[f"Species{i}"], B_socket)
        node.inputs["Order"].default_value = order
        node.inputs["Style"].default_value = style
        parent_tree.links.new(BondMesh.outputs["Geometry"], node.inputs["Geometry"])
        parent_tree.links.new(BondMesh.outputs["Length"], node.inputs["Length"])
        parent_tree.links.new(BondMesh.outputs["Rotation"], node.inputs["Rotation"])
        # link the outputs to parent node
        JoinGeometry = get_node_by_name(
            parent_tree.nodes,
            "%s_JoinGeometry" % self.label,
            "GeometryNodeJoinGeometry",
        )
        parent_tree.links.new(node.outputs["Geometry"], JoinGeometry.inputs["Geometry"])
        #
        InstanceOnPoint = get_node_by_name(
            nodes, "InstanceOnPoint_%s" % name, "GeometryNodeInstanceOnPoints"
        )
        ObjectInstancer = get_node_by_name(
            nodes, "ObjectInfo_%s" % name, "GeometryNodeObjectInfo"
        )
        ObjectInstancer.inputs["Object"].default_value = self.settings.instancers[
            sp["name"]
        ]["%s_%s" % (order, style)]
        #
        # order
        BondOrderAttribute = get_node_by_name(
            nodes,
            "%s_NamedAttribute_bond_order" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        BondOrderAttribute.inputs["Name"].default_value = "bond_order"
        BondOrderAttribute.data_type = "INT"
        CompareOrder = get_node_by_name(
            nodes, "CompareOrder_%s_%s" % (self.label, order), "FunctionNodeCompare"
        )
        CompareOrder.operation = "EQUAL"
        CompareOrder.data_type = "INT"
        socket = get_socket_by_identifier(CompareOrder, "B_INT")
        links.new(GroupInput.outputs["Order"], socket)
        CompareOrder.inputs[1].default_value = order
        output_socket = get_socket_by_identifier(
            BondOrderAttribute, "Attribute_Int", type="outputs"
        )
        input_socket = get_socket_by_identifier(CompareOrder, "A_INT")
        links.new(output_socket, input_socket)
        # style
        BondStyleAttribute = get_node_by_name(
            nodes,
            "%s_NamedAttribute_bond_style" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        BondStyleAttribute.inputs["Name"].default_value = "bond_style"
        BondStyleAttribute.data_type = "INT"
        CompareStyle = get_node_by_name(
            nodes, "CompareStyle_%s_%s" % (self.label, style), "FunctionNodeCompare"
        )
        CompareStyle.operation = "EQUAL"
        CompareStyle.data_type = "INT"
        socket = get_socket_by_identifier(CompareStyle, "B_INT")
        links.new(GroupInput.outputs["Style"], socket)
        CompareStyle.inputs[1].default_value = style
        output_socket = get_socket_by_identifier(
            BondStyleAttribute, "Attribute_Int", type="outputs"
        )
        input_socket = get_socket_by_identifier(CompareStyle, "A_INT")
        links.new(output_socket, input_socket)
        BoolSpecies = get_node_by_name(
            nodes, "%s_BooleanMath_species" % name, "FunctionNodeBooleanMath"
        )
        BoolOrder = get_node_by_name(
            nodes, "%s_BooleanMath_order" % name, "FunctionNodeBooleanMath"
        )
        BoolStyle = get_node_by_name(
            nodes, "%s_BooleanMath_style" % name, "FunctionNodeBooleanMath"
        )
        BoolModelStyle = get_node_by_name(
            nodes, "%s_BooleanMath_modelstyle" % name, "FunctionNodeBooleanMath"
        )
        BoolShow = get_node_by_name(
            nodes, "%s_BooleanMath_show" % name, "FunctionNodeBooleanMath"
        )
        BoolBondLength = get_node_by_name(
            nodes, "%s_BoolBondLength" % name, "FunctionNodeBooleanMath"
        )
        # bondlength larger than max will not show
        LessBondLength = get_node_by_name(
            nodes, "%s_LessBondLength" % name, "ShaderNodeMath"
        )
        LessBondLength.operation = "LESS_THAN"
        LessBondLength.inputs[1].default_value = sp["max"]
        #
        links.new(GroupInput.outputs["Geometry"], InstanceOnPoint.inputs["Points"])
        # scale
        CombineXYZ = get_node_by_name(
            nodes, "%s_CombineXYZ" % self.label, "ShaderNodeCombineXYZ"
        )
        CombineXYZ.inputs[0].default_value = 1
        CombineXYZ.inputs[1].default_value = 1
        links.new(GroupInput.outputs["Length"], CombineXYZ.inputs["Z"])
        # show
        BondShowAttribute = get_node_by_name(
            nodes,
            "%s_NamedAttribute_bond_show" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        BondShowAttribute.inputs["Name"].default_value = "bond_show"
        BondShowAttribute.data_type = "INT"
        socket = get_socket_by_identifier(
            BondShowAttribute, "Attribute_Int", type="outputs"
        )
        links.new(socket, BoolShow.inputs[0])
        # model style
        BondModelStyleAttribute = get_node_by_name(
            nodes,
            "%s_NamedAttribute_bond_model_style" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        BondModelStyleAttribute.inputs["Name"].default_value = "bond_model_style"
        BondModelStyleAttribute.data_type = "INT"
        socket = get_socket_by_identifier(
            BondModelStyleAttribute, "Attribute_Int", type="outputs"
        )
        links.new(socket, BoolModelStyle.inputs[0])
        for i in range(2):
            links.new(CompareSpecies[i].outputs[0], BoolSpecies.inputs[i])
        links.new(BoolSpecies.outputs[0], BoolOrder.inputs[0])
        links.new(CompareOrder.outputs[0], BoolOrder.inputs[1])
        links.new(BoolOrder.outputs[0], BoolStyle.inputs[0])
        links.new(CompareStyle.outputs[0], BoolStyle.inputs[1])
        links.new(BoolStyle.outputs[0], BoolModelStyle.inputs[1])
        links.new(BoolModelStyle.outputs[0], BoolShow.inputs[1])
        links.new(GroupInput.outputs["Length"], LessBondLength.inputs[0])
        links.new(BoolShow.outputs["Boolean"], BoolBondLength.inputs[0])
        links.new(LessBondLength.outputs[0], BoolBondLength.inputs[1])
        links.new(
            BoolBondLength.outputs["Boolean"], InstanceOnPoint.inputs["Selection"]
        )
        links.new(CombineXYZ.outputs[0], InstanceOnPoint.inputs["Scale"])
        #
        links.new(GroupInput.outputs["Rotation"], InstanceOnPoint.inputs["Rotation"])
        #
        links.new(
            ObjectInstancer.outputs["Geometry"], InstanceOnPoint.inputs["Instance"]
        )
        links.new(InstanceOnPoint.outputs["Instances"], GroupOutput.inputs["Geometry"])
        # print('Add geometry nodes for bonds: %s'%(time() - tstart))

    def update_geometry_node_instancer(self):
        """
        Make sure all pair has a geometry node flow
        and the instancer and material are updated.
        """
        tstart = time()
        #
        bondlists = self.arrays
        if len(bondlists["atoms_index0"]) == 0:
            return
        pairs = np.concatenate(
            (
                bondlists["species_index0"].reshape(-1, 1),
                bondlists["species_index1"].reshape(-1, 1),
                bondlists["bond_order"].reshape(-1, 1),
                bondlists["bond_style"].reshape(-1, 1),
            ),
            axis=1,
        )
        pairs = np.unique(pairs, axis=0)
        pairs = pairs.reshape(-1, 4)
        for pair in pairs:
            sp1 = number2String(pair[0])
            sp2 = number2String(pair[1])
            order = pair[2]
            style = pair[3]
            order_style = "%s_%s" % (order, style)
            sp = self.settings["%s-%s" % (sp1, sp2)]
            sp = sp.as_dict()
            # update geometry node
            if self.settings.instancers[sp["name"]][order_style] is None:
                self.settings.build_instancer(sp, order, style)
                self.add_bond_pair_node(sp)
            # always update materials
            self.settings.build_materials(sp, material_style=sp["material_style"])
            self.settings.assign_materials(sp, sp["order"], sp["style"])
            # compare radius
            if not np.isclose(
                sp["width"],
                self.settings.instancers[sp["name"]][order_style].Bbond.width,
            ):
                self.settings.build_instancer(sp)
                self.settings.assign_materials(sp, sp["order"], sp["style"])
            # update  instancers
            name = "Bond_%s_%s_%s_%s" % (
                self.label,
                sp["name"],
                sp["order"],
                sp["style"],
            )
            node = self.get_bond_pair_node(name)
            ObjectInstancer = get_node_by_name(
                node.node_tree.nodes, "ObjectInfo_%s" % name, "GeometryNodeObjectInfo"
            )
            ObjectInstancer.inputs["Object"].default_value = self.settings.instancers[
                sp["name"]
            ][order_style]
        logger.debug("update bond instancer: %s" % (time() - tstart))

    def update(self, bondlists=None, orders=None):
        """
        Draw bonds.
        calculate bond in all farmes, and merge all bondlists.
        Draw bonds in the bondlists, only show the bond if it
        is smaller than the max bond length.
        """

        object_mode()
        # clean_coll_objects(self.coll, 'bond')
        positions = self.batoms.positions
        trajectory = self.batoms.get_trajectory()["positions"]
        arrays = self.batoms.arrays
        boundary_data = self.batoms.boundary.boundary_data
        show = arrays["show"].astype(bool)
        species = arrays["species"][show]
        # trajectory_boundary = self.batoms.get_trajectory(self.batoms.batoms_boundary)
        # trajectory_search = self.batoms.get_trajectory(self.batoms.batoms_search)
        if self.batoms.nframe == 0:
            trajectory = np.array([positions])
        nframe = len(trajectory)
        bond_datas = {}
        tstart = time()
        setting = self.settings.as_dict()
        if bondlists is None:
            for f in range(nframe):
                # print('update bond: ', f)
                positions = trajectory[f, show, :]
                # build bondlist for unit cell
                (
                    bondlist,
                    bonddatas,
                    peciesBondDatas,
                    molPeciesDatas,
                ) = self.build_bondlists(
                    species, positions, self.batoms.cell, self.batoms.pbc, setting
                )
                # build bondlist for boundary atoms
                # for molecule with cell == [0, 0, 0], skip
                if boundary_data is not None:
                    bondlist = self.build_bondlists_with_boundary(
                        boundary_data,
                        bondlist,
                        bonddatas,
                        peciesBondDatas,
                        molPeciesDatas,
                    )
                    bondlist = self.check_boundary(bondlist)
                # search molecule
                # search bond
                if self.show_search:
                    self.search_bond.hide = False
                    self.search_bond.update(
                        bondlist,
                        peciesBondDatas,
                        molPeciesDatas,
                        arrays,
                        self.batoms.cell,
                    )
                if f == 0:
                    bondlists = bondlist
                else:
                    bondlists = np.append(bondlists, bondlist, axis=0)
            bondlists = np.unique(bondlists, axis=0)
        bond_datas = self.calc_bond_data(
            species,
            trajectory[:, show, :],
            self.batoms.cell,
            bondlists,
            self.settings,
            arrays["model_style"][show],
        )
        if orders:
            bond_datas.update({"bond_order": orders})
        if len(bond_datas) == 0:
            return
        self.set_arrays(bond_datas)
        logger.debug("draw bond: {0:10.2f} s".format(time() - tstart))

    def get_arrays(self):
        """ """
        object_mode()
        # tstart = time()
        arrays = self.get_attributes(domains=["EDGE"])
        arrays.update(
            {
                "positions": self.positions,
            }
        )
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    def set_arrays(self, arrays):
        """Add edge for bond.
        Set the bond properties."""
        # if not attributes, set them
        # in the case of mergeing objects, e.g. au111 += h2o,
        # the attributes of edges are removed.
        if "atoms_index0" not in self.obj.data.attributes:
            for att in default_bond_attributes:
                self.add_attribute(**att)
        # blender delete all edges
        dnvert = len(self.obj.data.edges)
        self.delete_edges_bmesh(range(dnvert))
        #
        edges = np.vstack((arrays["atoms_index0"], arrays["atoms_index1"])).T
        self.add_edges(edges, self.batoms.obj)
        self.set_attributes(arrays)
        # self.set_trajectory(arrays)
        self.update_geometry_node_instancer()
        self.update_geometry_nodes()

    def get_trajectory(self):
        """ """
        trajectory = {}
        trajectory["centers"] = self.get_obj_trajectory(self.obj)
        # trajectory['offsets'] = self.get_obj_trajectory(self.obj_o)
        return trajectory

    def set_trajectory(self, trajectory=None, frame_start=0):
        if trajectory is None:
            trajectory = self._trajectory
        nframe = len(trajectory["centers"])
        if nframe == 0:
            return
        name = "%s_bond" % (self.label)
        obj = self.obj
        self.set_shape_key(name, obj, trajectory["centers"], frame_start=frame_start)
        #
        # name = '%s_bond_offset'%(self.label)
        # obj = self.obj_o
        # self.set_shape_key(name, obj, trajectory['offsets'], frame_start=frame_start)

    @property
    def bondlists(self):
        return self.get_bondlists()

    def get_bondlists(self):
        """ """
        object_mode()
        # tstart = time()
        arrays = self.arrays
        i = arrays["atoms_index0"].reshape(-1, 1)
        j = arrays["atoms_index1"].reshape(-1, 1)
        p = arrays["polyhedra"].reshape(-1, 1)
        bond_offset0 = arrays["bond_offset0"]
        bond_offset1 = arrays["bond_offset1"]
        bondlists = np.concatenate((i, j, bond_offset0, bond_offset1, p), axis=1)
        # bondlists = bondlists.astype(int)
        # print('get_arrays: %s'%(time() - tstart))
        return bondlists

    @property
    def search_bond(self):
        """search_bond object."""
        if self._search_bond is not None:
            return self._search_bond
        search_bond = SearchBond(self.label, batoms=self.batoms)
        self._search_bond = search_bond
        return search_bond

    @property
    def show_search(self):
        return self.settings.coll.Bbond.show_search

    @show_search.setter
    def show_search(self, show_search):
        self.settings.coll.Bbond.show_search = show_search
        if not show_search:
            self.search_bond.set_arrays(default_search_bond_datas.copy())
            self._search_bond = None
        else:
            self.update()

    @property
    def show_hydrogen_bond(self):
        return self.settings.coll.Bbond.show_hydrogen_bond

    @show_hydrogen_bond.setter
    def show_hydrogen_bond(self, show_hydrogen_bond):
        self.settings.coll.Bbond.show_hydrogen_bond = show_hydrogen_bond
        self.update()

    def __getitem__(self, indices):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        from batoms.bond.slicebonds import SliceBonds

        slicebonds = SliceBonds(self.label, indices, bonds=self)
        return slicebonds

    def __setitem__(self, indices, value):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[indices] = value
        self.set_positions(positions)

    def repeat(self, m, cell):
        """
        In-place repeat of atoms.

        >>> from batoms.bond import Bbond
        >>> c = Bbond('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
        >>> c.repeat([3, 3, 3], np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]))
        """
        if isinstance(m, int):
            m = (m, m, m)
        for x, vec in zip(m, cell):
            if x != 1 and not vec.any():
                raise ValueError("Cannot repeat along undefined lattice " "vector")
        M = np.product(m)
        n = len(self)
        positions = np.tile(
            self.positions, (M,) + (1,) * (len(self.positions.shape) - 1)
        )
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
        self.add_vertices(positions[n:])

    def copy(self, label, species):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:

        >>> h_new = h.copy(label = 'h_new', species = 'H')

        """
        object_mode()
        bbond = self.__class__(
            label,
            species,
            self.local_positions,
            location=self.obj.location,
            scale=self.scale,
            material=self.material,
        )
        return bbond

    def extend(self, other):
        """
        Extend bbond object by appending bbond from *other*.

        >>> from batoms.bonds import Bbond
        >>> h1 = Bbond('h2o', 'H_1', [[0, 0, 0], [2, 0, 0]])
        >>> h2 = Bbond('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
        >>> h = h1 + h2
        """
        # could also use self.add_vertices(other.positions)
        object_mode()
        bpy.ops.object.select_all(action="DESELECT")
        self.obj.select_set(True)
        other.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.join()

    def __iadd__(self, other):
        """
        >>> h1 += h2
        """
        self.extend(other)
        return self

    def __add__(self, other):
        """
        >>> h1 = h1 + h2
        """
        self += other
        return self

    def __iter__(self):
        bbond = self.obj
        for i in range(len(self)):
            yield bbond.matrix_world @ bbond.data.vertices[i].co

    def __repr__(self):
        s = "Bonds(Total: {:6d}, {}".format(len(self), self.arrays)
        return s

    def build_bondlists(self, species, positions, cell, pbc, setting):
        """
        build bondlist for atoms
        steps:
        1 build bondlist for atoms with pbc
        2 search connected_components (molecule),
          return pecies data and its neighbour
        3 add bondlist related with molecule

        """
        from batoms.neighborlist import bondlist_kdtree

        bondlists = np.zeros((0, 11), dtype=int)
        bonddatas = {}
        if len(setting) == 0:
            return bondlists, bonddatas, {}, {}
        #
        # tstart = time()
        # ==========================================================
        # step 1 build bondlist for atoms with pbc
        # nli: index1
        # nlj: index2
        # nlk: search bond style
        # nlp: polyhedra
        # nlt: bond type: hydrogen bond
        # nlSj: offset of atoms in nlj
        nli, nlj, nlk, nlp, nlt, nlSj = bondlist_kdtree(
            "ijkptS", species, positions, cell, pbc, setting
        )
        nb = len(nli)
        nlSi = np.zeros((nb, 3))
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        # TODO search0 is not used, please check why
        # search type 0,
        # search0 = np.where(
        # (nlk == 0) & (nlSj != np.array([0, 0, 0])).any(axis=1), False, True
        # )
        # 0  1  2:5       5:8          8     9
        # i, j, offset_i, offset_j, search, search_style
        bondlists = np.concatenate(
            (
                np.array([nli, nlj]).T,
                np.array(nlSi, dtype=int),
                np.array(nlSj),
                np.array(nlk).reshape(-1, 1),
                np.array(nlp).reshape(-1, 1),
                np.array(nlt).reshape(-1, 1),
            ),
            axis=1,
        )
        bondlists = bondlists.astype(int)
        # remove bond outside box for search0
        # not now, we need all bonds here, and then add a final check.
        # bondlists = bondlists[search0]
        # build bondatas, for each atom, save the bonds connect to it.
        argsort = bondlists[:, 0].argsort()
        bondlists = bondlists[argsort]
        u, indices = np.unique(bondlists[:, 0], return_index=True)
        indices = np.append(indices, len(bondlists))
        m = len(u)
        for i in range(m):
            bonddatas[u[i]] = bondlists[indices[i] : indices[i + 1]]
        #
        # ===================================================
        # 2 search connected_components (molecule),
        #   return pecies data and its neighbour
        peciesBondLists, molPeciesDatas = self.build_peciesBondLists(
            len(positions), bondlists
        )
        # build peciesBondDatas
        peciesBondDatas = {}
        for p in molPeciesDatas:
            peciesBondDatas[p] = []
        # since this is search bond type 2, bothway, first i
        argsort = peciesBondLists[:, 0].argsort()
        peciesBondLists = peciesBondLists[argsort]
        u, indices = np.unique(peciesBondLists[:, 0], return_index=True)
        indices = np.append(indices, len(peciesBondLists))
        m = len(u)
        # for each pecies, save its neighbour and offset
        for i in range(m):
            peciesBondDatas[u[i]] = peciesBondLists[indices[i] : indices[i + 1]]
        # bothway, then j
        argsort = peciesBondLists[:, 1].argsort()
        peciesBondLists = peciesBondLists[argsort]
        u, indices = np.unique(peciesBondLists[:, 1], return_index=True)
        indices = np.append(indices, len(peciesBondLists))
        m = len(u)
        # for each pecies, save its neighbour and offset
        for i in range(m):
            data = peciesBondLists[indices[i] : indices[i + 1]]
            data = data.reshape(-1, 11)
            data[:, 5:8] *= -1
            data[:, [0, 1]] = data[:, [1, 0]]
            peciesBondDatas[u[i]] = data
        # ==========================================================
        # 3 add bondlist related with molecule
        # add bondlist for molecule
        for mol in peciesBondLists:
            indices = molPeciesDatas[mol[1]]
            for i in indices:
                if i not in bonddatas:
                    continue
                data = bonddatas[i]
                n = len(data)
                bondlists = np.append(bondlists, data, axis=0)
                bondlists[-n:, 2:5] += mol[5:8]
                bondlists[-n:, 5:8] += mol[5:8]
            indices = molPeciesDatas[mol[0]]
            for i in indices:
                if i not in bonddatas:
                    continue
                data = bonddatas[i]
                n = len(data)
                bondlists = np.append(bondlists, data, axis=0)
                bondlists[-n:, 2:5] -= mol[5:8]
                bondlists[-n:, 5:8] -= mol[5:8]
        # print(bondlists)
        bondlists = bondlists.astype(int)
        bondlists = np.unique(bondlists, axis=0)
        self.peciesBondLists = peciesBondLists
        self.molPeciesDatas = molPeciesDatas
        self.peciesBondDatas = peciesBondDatas
        return bondlists, bonddatas, peciesBondDatas, molPeciesDatas

    def build_peciesBondLists(self, natom, bondlists):
        """
        search type 2: build molecule
        steps:
        1 search connected_components (molecules) inside atoms
        2 construct the molecules by its pecies and the coresponding offsets
        3 for
        3 return the molecules
        """
        from scipy.sparse import csgraph, csr_matrix

        molPeciesDatas = {}
        peciesBondLists = np.zeros((0, 11), dtype=int)
        # search type 2
        k = bondlists[:, 8]
        indices = np.where(k == 2)[0]
        ns2 = len(indices)
        if ns2 == 0:
            return peciesBondLists, molPeciesDatas
        bondlists1 = bondlists[indices, :]
        # ========================================================
        # 1 search connected_components (molecules) inside atoms
        # with crossed bond, entire molecule
        molDatas = {}
        ai = bondlists1[:, 0]
        aj = bondlists1[:, 1]
        data = np.ones(ns2, dtype=int)
        matrix = csr_matrix((data, (ai, aj)), shape=(natom, natom))
        n_components1, component_list1 = csgraph.connected_components(matrix)
        for i in range(n_components1):
            indices = np.where(component_list1 == i)[0]
            n = len(indices)
            if n < 2:
                continue
            molDatas[i] = {"sub": []}
            molDatas[i]["indices"] = indices
            # TODO: check here
            molDatas[i]["offsets"] = indices
        # without crossed bond, small pecies of molDatas
        mask = np.where(
            (bondlists1[:, 2:5] != bondlists1[:, 5:8]).any(axis=1), False, True
        )
        data = data[mask]
        ai = ai[mask]
        aj = aj[mask]
        matrix = csr_matrix((data, (ai, aj)), shape=(natom, natom))
        n_components2, component_list2 = csgraph.connected_components(matrix)
        #
        for i in range(n_components2):
            indices = np.where(component_list2 == i)[0]
            if component_list1[indices[0]] in molDatas:
                # this pecies belong to molDatas
                molDatas[component_list1[indices[0]]]["sub"].append(i)
                molPeciesDatas[i] = indices
        # pprint(molPeciesDatas)
        # cross box bond, find direct neighbour pecies
        bondlists2 = bondlists1[~mask]
        n = len(bondlists2)
        ai = bondlists2[:, 0].astype(int)
        aj = bondlists2[:, 1].astype(int)
        peciesBondLists = np.zeros((n, 11), dtype=int)
        molBondDicts = {}
        for i in range(n):
            # molDatasId = component_list1[bondlists2[i, 0]]
            i1 = component_list2[ai[i]]
            j1 = component_list2[aj[i]]
            peciesBondLists[i, 0] = i1
            peciesBondLists[i, 1] = j1
            peciesBondLists[i, 2:] = bondlists2[i, 2:]
            molBondDicts[(i1, j1)] = bondlists2[i, 2:]
        peciesBondLists = np.unique(peciesBondLists, axis=0)
        # ========================================================
        # find indirect neighbour pecies inside mol1
        # find connected path between pecies
        ai = peciesBondLists[:, 0].astype(int)
        aj = peciesBondLists[:, 1].astype(int)
        nml = len(ai)
        data = np.ones(nml, dtype=int)
        graph = csr_matrix((data, (ai, aj)), shape=(n_components2, n_components2))
        dist_matrix, predecessors = csgraph.shortest_path(
            csgraph=graph, return_predecessors=True
        )
        # print('dist_matrix: ', dist_matrix)
        # print('predecessors: ', predecessors)
        dist_matrix = dist_matrix.astype(int)
        mollist = np.zeros(11, dtype=int)
        for i, data in molDatas.items():
            indices = data["sub"]
            n = len(indices)
            if n < 2:
                continue
            for j in range(n - 1):
                for k in range(j + 1, n):
                    if dist_matrix[indices[j], indices[k]] == 1:
                        continue
                    path = [indices[k]]
                    end = indices[k]
                    last = indices[k]
                    for i1 in range(dist_matrix[indices[j], end]):
                        last = predecessors[indices[j], last]
                        path = [last] + path
                    offsets = np.zeros(3)
                    for i2 in range(0, len(path) - 1):
                        offsets += molBondDicts[(path[i2], path[i2 + 1])][3:6]
                    mollist[0] = path[0]
                    mollist[1] = path[-1]
                    mollist[5:8] = offsets
                    peciesBondLists = np.append(
                        peciesBondLists, np.array([mollist]), axis=0
                    )
        self.molDatas = molDatas
        return peciesBondLists, molPeciesDatas

    def build_bondlists_with_boundary(
        self, arrays, bondlists, bonddatas, peciesBondDatas, molPeciesDatas
    ):
        """
        build extra bondlists based on boundary atoms
        """
        n = len(arrays["positions"])
        if n == 0:
            return bondlists
        # search bond type 0 and 1
        for i in range(n):
            if arrays["indices"][i] not in bonddatas:
                # print('no bonds: ', arrays['indices'][i])
                continue
            data = bonddatas[arrays["indices"][i]]
            n = len(data)
            bondlists = np.append(bondlists, data, axis=0)
            bondlists[-n:, 2:5] = bondlists[-n:, 2:5] + arrays["boundary_offset"][i]
            bondlists[-n:, 5:8] = bondlists[-n:, 5:8] + arrays["boundary_offset"][i]
            # todo: in this case, some of the atoms overlap
            # with original atoms.
            # since it doesn't influence the 3d view, we
            # just leave it like this
            # bondlists[-n:, 8] = 1
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        bondlists = np.unique(bondlists, axis=0)
        # search bond type 2
        # divide boudanry atoms by boundary_offset
        boundary_offset = np.unique(arrays["boundary_offset"], axis=0)
        n = len(boundary_offset) * len(molPeciesDatas)
        peciesArrays = np.zeros((n, 4), dtype=int)
        n = 0
        for offset in boundary_offset:
            indices = np.where((arrays["boundary_offset"] == offset).all(axis=1))[0]
            indices = arrays["indices"][indices]
            for p, indices0 in molPeciesDatas.items():
                indices1 = np.intersect1d(indices0, indices)
                if len(indices1) > 0:
                    peciesArrays[n][0] = p
                    peciesArrays[n][1:4] = offset
                    n += 1
        peciesArrays = peciesArrays[:n]
        npa = len(peciesArrays)
        if npa == 0:
            return bondlists
        # search bond type 2
        for i in range(npa):
            peciesId = peciesArrays[i, 0]
            p = peciesBondDatas[peciesId]
            # repeat itself
            indices = molPeciesDatas[peciesId]
            for j in indices:
                if j not in bonddatas:
                    continue
                data = bonddatas[j]
                n = len(data)
                bondlists = np.append(bondlists, data, axis=0)
                bondlists[-n:, 2:5] += peciesArrays[i, 1:4]
                bondlists[-n:, 5:8] += peciesArrays[i, 1:4]
            # repeat its neighbour pecies
            for pb in p:
                peciesId = pb[1]
                indices = molPeciesDatas[peciesId]
                for j in indices:
                    if j not in bonddatas:
                        continue
                    data = bonddatas[j]
                    n = len(data)
                    bondlists = np.append(bondlists, data, axis=0)
                    bondlists[-n:, 2:5] += peciesArrays[i, 1:4] + pb[5:8]
                    bondlists[-n:, 5:8] += peciesArrays[i, 1:4] + pb[5:8]
            # todo: in this case, some of the atoms overlap with
            # original atoms.
            # since it doesn't influence the 3d view,
            # we just leave it like this
            # bondlists[-n:, 8] = 1
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        bondlists = np.unique(bondlists, axis=0)
        # search bond type 2

        return bondlists

    def check_boundary(self, bondlists, eps=1e-6):
        """check boundary for bond search 0

        Args:
            eps (float): default 1e-6
        """
        from ase.geometry import complete_cell

        arrays = self.batoms.arrays
        positions = arrays["positions"]
        cell = self.batoms.cell
        # get scaled positions
        positions = np.linalg.solve(complete_cell(cell).T, positions.T).T
        npositions = positions[bondlists[:, 1]] + bondlists[:, 5:8]
        boundary = self.batoms.boundary.boundary
        # boundary condition
        mask1 = np.where(
            (npositions[:, 0] > boundary[0][0] - eps)
            & (npositions[:, 0] < boundary[0][1] + eps)  # noqa W503
            & (npositions[:, 1] > boundary[1][0] - eps)  # noqa W503
            & (npositions[:, 1] < boundary[1][1] + eps)  # noqa W503
            & (npositions[:, 2] > boundary[2][0] - eps)  # noqa W503
            & (npositions[:, 2] < boundary[2][1] + eps),  # noqa W503
            False,
            True,
        )
        mask2 = np.where(bondlists[:, 8] == 0, True, False)
        mask = ~(mask1 & mask2)
        bondlists = bondlists[mask]
        return bondlists

    def search_molecule(self, natom, bondlists):
        """ """
        from scipy.sparse import csgraph, csr_matrix

        mols = {}
        # search type
        k = bondlists[:, -1]
        indices = np.where(k == 2)[0]
        ns2 = len(indices)
        if ns2 == 0:
            return mols
        i = bondlists[indices, 0]
        j = bondlists[indices, 1]
        data = np.ones(ns2, dtype=int)
        matrix = csr_matrix((data, (i, j)), shape=(natom, natom))
        # print(matrix)
        n_components, component_list = csgraph.connected_components(matrix)
        # print(n_components)
        # print(component_list)
        for i in range(n_components):
            indices = np.where(component_list == i)[0]
            if len(indices) < 2:
                continue
            mols[i] = indices
        return mols

    def calc_bond_data(
        self, speciesarray, positions, cell, bondlists, bondsetting, model_styles
    ):
        """
        positions: 3-demision
        todo: support frame for positions and offsets.
        """
        tstart = time()
        if not self.show_hydrogen_bond:
            bondlists = bondlists[bondlists[:, 10] != 1]
        # properties
        atoms_index0 = np.array(bondlists[:, 0], dtype=int)
        atoms_index1 = np.array(bondlists[:, 1], dtype=int)
        # second bond for high order bond
        nb = len(atoms_index0)
        second_bond = np.roll(np.arange(nb), 1)
        # atoms_index2 and atoms_index3 will be removed for Blender 3.1
        atoms_index2 = np.roll(atoms_index0, 1)
        atoms_index3 = np.roll(atoms_index1, 1)
        shows = np.ones(nb, dtype=int)
        orders = np.zeros(nb, dtype=int)  # 1, 2, 3
        styles = np.zeros(nb, dtype=int)  # 0, 1, 2
        widths = np.ones(nb, dtype=float)
        species_index0 = np.zeros(nb, dtype=int)
        species_index1 = np.zeros(nb, dtype=int)
        bond_model_styles = model_styles[atoms_index0]
        if nb == 0:
            return default_bond_datas
        # ------------------------------------
        # offsets
        bond_offset0 = np.dot(bondlists[:, 2:5], cell)
        bond_offset1 = np.dot(bondlists[:, 5:8], cell)
        # polyhedra
        polyhedras = np.array(bondlists[:, 9], dtype=int)
        # bond center
        if len(positions.shape) == 2:
            positions = np.array([positions])
        # ---------------------------------------------
        for b in bondsetting:
            spi = b.species1
            spj = b.species2
            indi = speciesarray[atoms_index0] == spi
            indj = speciesarray[atoms_index1] == spj
            ind = indi & indj
            orders[ind] = b.order
            styles[ind] = int(b.style)
            widths[ind] = b.width
            polyhedras[ind] = b.polyhedra
            species_index0[ind] = string2Number(spi)
            species_index1[ind] = string2Number(spj)
            if b.order > 1:
                second_bond, atoms_index2, atoms_index3 = self.secondBond(
                    b, speciesarray, second_bond, atoms_index2, atoms_index3, bondlists
                )
        # offsets for second bond
        # this will be remove for Blender 3.1
        bond_offset2 = np.dot(bondlists[:, 2:5][second_bond], cell)
        bond_offset3 = np.dot(bondlists[:, 5:8][second_bond], cell)
        datas = {
            "atoms_index0": atoms_index0,
            "atoms_index1": atoms_index1,
            "atoms_index2": atoms_index2,
            "atoms_index3": atoms_index3,
            "second_bond": second_bond,
            "species_index0": species_index0,
            "species_index1": species_index1,
            "bond_offset0": bond_offset0,
            "bond_offset1": bond_offset1,
            "bond_offset2": bond_offset2,
            "bond_offset3": bond_offset3,
            # 'widths': widths,
            "bond_show": shows,
            "bond_order": orders,
            "bond_style": styles,
            "bond_model_style": bond_model_styles,
            "polyhedra": polyhedras,
        }
        # print('datas: ', datas)
        logger.debug("calc_bond_data: {0:10.2f} s".format(time() - tstart))
        return datas

    def secondBond(
        self, b, speciesarray, second_bond, atoms_index2, atoms_index3, bondlists
    ):
        """
        determine the plane of high order bond.
        v1 = atoms2 - atoms1
        v2 = atoms4 - atoms3
        normal = np.cross(v1, v2)
        """
        spi = b.species1
        spj = b.species2
        # find all bonds belong to this pair spi--spj
        # here use &
        indi = speciesarray[bondlists[:, 0]] == spi
        indj = speciesarray[bondlists[:, 1]] == spj
        indices1 = np.where(indi & indj)[0]
        # find all bonds connect to this bond pair
        # here use |
        indi1 = speciesarray[bondlists[:, 1]] == spi
        indj2 = speciesarray[bondlists[:, 0]] == spj
        indices2 = np.where(indi | indj | indi1 | indj2)[0]
        bondlists2 = bondlists[indices2]
        for i in indices1:
            # find another bond for this bond: bondlists[i]
            indices3 = np.where(
                (bondlists2[:, 0] == bondlists[i, 0])  # noqa W503
                | (bondlists2[:, 1] == bondlists[i, 0])  # noqa W503
                | (bondlists2[:, 0] == bondlists[i, 1])  # noqa W503
                | (bondlists2[:, 1] == bondlists[i, 1])  # noqa W503
            )[0]
            indices3 = indices2[indices3]
            if len(indices3) > 1:
                indices4 = np.where(indices3 != i)[0]
                second_bond[i] = indices3[indices4[0]]
                atoms_index2[i] = bondlists[second_bond[i], 0]
                atoms_index3[i] = bondlists[second_bond[i], 1]
        return second_bond, atoms_index2, atoms_index3

    def high_order_bond_plane(
        self, b, speciesarray, positions, nvec, offsets, bondlists
    ):
        """
        determine the plane of high order bond
        """
        spi = b.species1
        spj = b.species2
        indi = speciesarray[bondlists[:, 0]] == spi
        indj = speciesarray[bondlists[:, 1]] == spj
        bondlists1 = bondlists[indi & indj]
        nbond = len(bondlists1)
        indi1 = speciesarray[bondlists[:, 1]] == spi
        indj2 = speciesarray[bondlists[:, 0]] == spj
        bondlists2 = bondlists[indi | indj | indi1 | indj2]
        for i in range(nbond):
            # find another bond
            mask = np.logical_not(
                (
                    (bondlists2[:, 0] == bondlists1[i, 0])
                    & (bondlists2[:, 1] == bondlists1[i, 1])  # noqa W503
                )
                | (  # noqa W503
                    (bondlists2[:, 0] != bondlists1[i, 0])
                    & (bondlists2[:, 0] != bondlists1[i, 1])  # noqa W503
                    & (bondlists2[:, 1] != bondlists1[i, 0])  # noqa W503
                    & (bondlists2[:, 1] != bondlists1[i, 1])  # noqa W503
                )
            )
            localbondlist = bondlists2[mask]
            if len(localbondlist) == 0:
                second_bond = np.array([0.0, 0.0, 1])
            else:
                second_bond = (
                    positions[localbondlist[0, 0]] - positions[localbondlist[0, 1]]
                )
            norml = np.cross(second_bond, nvec[i]) + np.array([1e-8, 0, 0])
            offset = np.cross(norml, nvec[i]) + np.array([1e-8, 0, 0])
            offsets[i] = offset / np.linalg.norm(offset)
            # print(offsets[i], offset)

        return offsets

    def bond_order_auto_set(self):
        """Set bond order by pybel
        #TODO support show, select
        """
        from openbabel import pybel

        species = self.batoms.arrays["species"]

        mol = self.batoms.as_pybel(export_bond=False)
        bondlists = [
            [
                b.GetBeginAtom().GetIndex(),
                b.GetEndAtom().GetIndex(),
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
            for b in pybel.ob.OBMolBondIter(mol.OBMol)
        ]
        # TODO detect hydrogen bond
        cutoff_dict = self.settings.cutoff_dict
        nbond = len(bondlists)
        remove = []
        bondlists = np.array(bondlists)
        for i in range(nbond):
            b = bondlists[i]
            if (species[b[0]], species[b[1]]) not in cutoff_dict:
                # swtich bond
                if (species[b[1]], species[b[0]]) in cutoff_dict:
                    bondlists[i][0], bondlists[i][1] = bondlists[i][1], bondlists[i][0]
                else:
                    # remove bond
                    remove.append(i)
        np.delete(bondlists, remove)
        orders = [b.GetBondOrder() for b in pybel.ob.OBMolBondIter(mol.OBMol)]
        self.update(bondlists, orders=orders)

    @property
    def setting(self):
        from batoms.utils import deprecated

        """setting object."""
        deprecated(
            '"setting" will be deprecated in the furture, please use "settings".'
        )
        return self.settings

    def as_dict(self):
        """ """
        data = {"array": None}
        if len(self) > 1:
            data["array"] = dict(self.arrays)
        data["show_search"] = self.show_search
        data["show_hydrogen_bond"] = self.show_hydrogen_bond
        data["settings"] = self.settings.as_dict()
        return data
