"""
# TODO add locaiton in geometry node
# TODO rearrange code to handle offsets for all bonds, boundary
       search_bond and so on
"""
from time import time

import bpy
import numpy as np
from ase.geometry import complete_cell
from batoms.base.object import ObjectGN
from batoms.utils import number2String, string2Number
from batoms.utils.butils import compareNodeType
import logging

# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

default_attributes = [
    {"name": "atoms_index", "data_type": "INT"},
    {"name": "species_index", "data_type": "INT"},
    {"name": "show", "data_type": "INT"},
    {"name": "select", "data_type": "INT"},
    {"name": "model_style", "data_type": "INT"},
    {"name": "scale", "data_type": "FLOAT"},
    {"name": "radius_style", "data_type": "INT"},
    {"name": "boundary_offset", "data_type": "FLOAT_VECTOR"},
]

default_GroupInput = [
    ["atoms_index", "NodeSocketInt"],
    ["species_index", "NodeSocketInt"],
    ["show", "NodeSocketBool"],
    ["select", "NodeSocketInt"],
    ["model_style", "NodeSocketInt"],
    ["scale", "NodeSocketFloat"],
    ["radius_style", "NodeSocketInt"],
    ["boundary_offset", "NodeSocketVector"],
]

default_boundary_datas = {
    "atoms_index": np.ones(0, dtype=int),
    "species_index": np.ones(0, dtype=int),
    "species": np.ones(0, dtype="U4"),
    "positions": np.zeros((0, 3)),
    "scales": np.zeros(0),
    "boundary_offset": np.zeros((0, 3)),
    "model_styles": np.ones(0, dtype=int),
    "shows": np.ones(0, dtype=int),
    "selects": np.ones(0, dtype=int),
}


class Boundary(ObjectGN):
    def __init__(
        self,
        label=None,
        batoms=None,
        boundary=None,
        boundary_datas=None,
        location=(0, 0, 0),
        load=False,
    ):
        """_summary_

        Args:
            label (str, optional):
                _description_. Defaults to None.
            boundary (list, optional):
                search atoms at the boundary.
                Defaults to np.array([[-0.0, 1.0], [-0.0, 1.0], [-0.0, 1.0]]).
            boundary_datas (_type_, optional):
                _description_. Defaults to None.
            batoms (_type_, optional):
                _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = "boundary"
        ObjectGN.__init__(self, label, name)
        if not load or not self.loadable():
            if boundary_datas is not None:
                self.build_object(boundary_datas)  # , location=location)
            else:
                self.build_object(default_boundary_datas)  # , location=location)
            if boundary is not None:
                self.boundary = boundary

    def loadable(self):
        """Check loadable or not"""
        # object exist
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            return False
        # batoms exist, and flag is True
        coll = bpy.data.collections.get(self.label)
        if coll is None:
            return False
        return coll.batoms.boundary.flag

    def build_object(self, datas, location=[0, 0, 0], attributes={}):
        """
        build child object and add it to main objects.
        """
        # tstart = time()
        if len(datas["positions"].shape) == 2:
            self._trajectory = {
                "positions": np.array([datas["positions"]]),
            }
            positions = datas["positions"]
        elif len(datas["positions"].shape) == 3:
            self._trajectory = {
                "positions": datas["positions"],
            }
            positions = datas["positions"][0]
        else:
            raise Exception("Shape of positions is wrong!")
        #
        attributes.update(
            {
                "atoms_index": datas["atoms_index"],
                "species_index": datas["species_index"],
                "show": datas["shows"],
                "model_style": datas["model_styles"],
                "select": datas["selects"],
                "scale": datas["scales"],
            }
        )
        name = "%s_boundary" % self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(positions, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        # Add attributes
        for att in default_attributes:
            self.add_attribute(**att)
        self.batoms.coll.objects.link(obj)
        obj.location = location
        obj.batoms.type = "BOUNDARY"
        obj.batoms.label = self.batoms.label
        obj.parent = self.batoms.obj
        self.batoms.coll.batoms.boundary.flag = True
        #
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.init_geometry_node_modifier(default_GroupInput)
        self.build_geometry_node()
        self.set_trajectory()
        # print('boundary: build_object: {0:10.2f} s'.format(time() - tstart))

    def build_geometry_node(self):
        """ """
        from batoms.utils.butils import get_node_by_name

        links = self.gn_node_group.links
        nodes = self.gn_node_group.nodes
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        # ------------------------------------------------------------------
        JoinGeometry = get_node_by_name(
            nodes, "%s_JoinGeometry" % self.label, "GeometryNodeJoinGeometry"
        )
        links.new(GroupInput.outputs["Geometry"], JoinGeometry.inputs["Geometry"])
        links.new(JoinGeometry.outputs["Geometry"], GroupOutput.inputs["Geometry"])
        # ------------------------------------------------------------------
        # transform postions of batoms to boundary
        ObjectBatoms = get_node_by_name(
            nodes, "%s_ObjectBatoms" % self.label, "GeometryNodeObjectInfo"
        )
        ObjectBatoms.inputs["Object"].default_value = self.batoms.obj
        PositionBatoms = get_node_by_name(
            nodes, "%s_PositionBatoms" % (self.label), "GeometryNodeInputPosition"
        )
        TransferBatoms = get_node_by_name(
            nodes, "%s_TransferBatoms" % (self.label), "GeometryNodeSampleIndex"
        )
        TransferBatoms.data_type = "FLOAT_VECTOR"
        links.new(ObjectBatoms.outputs["Geometry"], TransferBatoms.inputs[0])
        links.new(PositionBatoms.outputs["Position"], TransferBatoms.inputs[3])
        links.new(GroupInput.outputs[1], TransferBatoms.inputs["Index"])
        # ------------------------------------------------------------------
        # calculate offset for boundary atoms
        OffsetAttribute = get_node_by_name(
            nodes,
            "%s_boundary_offset" % (self.label),
            "GeometryNodeInputNamedAttribute",
        )
        OffsetAttribute.inputs["Name"].default_value = "boundary_offset"
        OffsetAttribute.data_type = "FLOAT_VECTOR"
        # get arrays of cell
        cell_node = self.get_cell_node(self.gn_node_group)
        cell_node.inputs["Cell"].default_value = self.batoms.cell.obj
        dot_node = self.vector_dot_cell(self.gn_node_group)
        links.new(OffsetAttribute.outputs[0], dot_node.inputs["Vector"])
        for i in range(3):
            links.new(
                cell_node.outputs["A%d" % (i + 1)], dot_node.inputs["Array%d" % (i + 1)]
            )
        # we need one add operation to get the positions with offset
        VectorAdd = get_node_by_name(
            nodes, "%s_VectorAdd" % (self.label), "ShaderNodeVectorMath"
        )
        # ------------------------------------------------------------------
        # add positions with offsets
        VectorAdd.operation = "ADD"
        links.new(TransferBatoms.outputs[2], VectorAdd.inputs[0])
        links.new(dot_node.outputs[0], VectorAdd.inputs[1])
        # set positions
        SetPosition = get_node_by_name(
            nodes, "%s_SetPosition" % self.label, "GeometryNodeSetPosition"
        )
        links.new(GroupInput.outputs["Geometry"], SetPosition.inputs["Geometry"])
        links.new(VectorAdd.outputs[0], SetPosition.inputs["Position"])
        #
        # ------------------------------------------------------------------
        # transform scale of batoms to boundary
        if bpy.app.version_string >= "3.2.0":
            ScaleBatoms = get_node_by_name(
                nodes,
                "%s_ScaleBatoms" % (self.label),
                "GeometryNodeInputNamedAttribute",
            )
            # need to be "FLOAT_VECTOR", because scale is "FLOAT_VECTOR"
            ScaleBatoms.data_type = "FLOAT_VECTOR"
            ScaleBatoms.inputs[0].default_value = "scale"
            TransferScale = get_node_by_name(
                nodes, "%s_TransferScale" % (self.label), "GeometryNodeSampleIndex"
            )
            TransferScale.data_type = "FLOAT_VECTOR"
            links.new(ObjectBatoms.outputs["Geometry"], TransferScale.inputs[0])
            links.new(ScaleBatoms.outputs["Attribute"], TransferScale.inputs[3])
            links.new(GroupInput.outputs[1], TransferScale.inputs["Index"])

    def vector_dot_cell(self, parent_tree):
        """Dot product of a vector (1x3) and a cell (3x3).
        Note, the vector could be a vector field, (Nx3)"""
        from batoms.utils.butils import (
            get_node_by_name,
            create_node_tree,
        )

        default_interface = [
            ["Vector", "NodeSocketVector", "INPUT"],
            ["Array1", "NodeSocketVector", "INPUT"],
            ["Array2", "NodeSocketVector", "INPUT"],
            ["Array3", "NodeSocketVector", "INPUT"],
            ["Vector", "NodeSocketVector", "OUTPUT"],
        ]
        name = "Vector_dot_cell_%s" % (self.label)
        node = get_node_by_name(parent_tree.nodes, name=name, type="GeometryNodeGroup")
        node_tree = create_node_tree(name=name, interface=default_interface)
        node.node_tree = node_tree
        nodes = node_tree.nodes
        links = node_tree.links
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        #
        SeparateXYZs = []
        CombineXYZs = []
        DotProdcuts = []
        for i in range(3):
            # SeparateXYZ to get the transpose of the cell
            SeparateXYZ = get_node_by_name(
                nodes, f"{self.label}_SeparateXYZ_{i}", "ShaderNodeSeparateXYZ"
            )
            # Dot product of the vector and the array
            DotProdcut = get_node_by_name(
                nodes, f"{self.label}_DotProduct_{i}", "ShaderNodeVectorMath"
            )
            DotProdcut.operation = "DOT_PRODUCT"
            # then combine the vectors, and output to the group output
            CombineXYZ = get_node_by_name(
                nodes, f"{self.label}_CombineXYZ_{i}", "ShaderNodeCombineXYZ"
            )
            links.new(GroupInput.outputs[f"Array{i + 1}"], SeparateXYZ.inputs["Vector"])
            links.new(GroupInput.outputs["Vector"], DotProdcut.inputs[0])
            links.new(CombineXYZ.outputs["Vector"], DotProdcut.inputs[1])
            CombineXYZs.append(CombineXYZ)
            SeparateXYZs.append(SeparateXYZ)
            DotProdcuts.append(DotProdcut)
        for i in range(3):
            for j in range(3):
                links.new(SeparateXYZs[i].outputs[j], CombineXYZs[j].inputs[i])
        # output final vector
        CombineXYZ = get_node_by_name(
            nodes, f"{self.label}_CombineXYZ", "ShaderNodeCombineXYZ"
        )
        for i in range(3):
            links.new(DotProdcuts[i].outputs["Value"], CombineXYZ.inputs[i])
        links.new(CombineXYZ.outputs["Vector"], GroupOutput.inputs["Vector"])
        return node

    def get_cell_node(self, parent_tree):
        """Get the position of the cell."""
        from batoms.utils.butils import (
            get_socket_by_identifier,
            get_node_by_name,
            create_node_tree,
        )

        default_interface = [
            ["Cell", "NodeSocketObject", "INPUT"],
            ["A1", "NodeSocketVector", "OUTPUT"],
            ["A2", "NodeSocketVector", "OUTPUT"],
            ["A3", "NodeSocketVector", "OUTPUT"],
        ]
        name = "Cell_Array_%s" % (self.label)
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
            nodes, f"{self.label}_CellObject", "GeometryNodeObjectInfo"
        )
        links.new(GroupInput.outputs["Cell"], CellObject.inputs["Object"])
        Position = get_node_by_name(
            nodes, "%s_Position" % (self.label), "GeometryNodeInputPosition"
        )
        for i in range(3):
            PositionAtIndex = get_node_by_name(
                nodes, f"{self.label}_PositionAtIndex_{i}", "GeometryNodeSampleIndex"
            )
            PositionAtIndex.data_type = "FLOAT_VECTOR"
            PositionAtIndex.inputs["Index"].default_value = i + 1
            links.new(CellObject.outputs["Geometry"], PositionAtIndex.inputs[0])
            input_socket = get_socket_by_identifier(PositionAtIndex, "Value_Vector")
            links.new(Position.outputs["Position"], input_socket)
            output_socket = get_socket_by_identifier(
                PositionAtIndex, "Value_Vector", type="outputs"
            )
            links.new(output_socket, GroupOutput.inputs["A%d" % (i + 1)])

        return node

    def add_geometry_node(self, spname):
        """ """
        from batoms.utils.butils import get_node_by_name

        links = self.gn_node_group.links
        nodes = self.gn_node_group.nodes
        GroupInput = nodes[0]
        SetPosition = get_node_by_name(nodes, "%s_SetPosition" % self.label)
        JoinGeometry = get_node_by_name(
            nodes, "%s_JoinGeometry" % self.label, "GeometryNodeJoinGeometry"
        )
        CompareSpecies = get_node_by_name(
            nodes, "CompareFloats_%s_%s" % (self.label, spname), compareNodeType
        )
        CompareSpecies.operation = "EQUAL"
        # CompareSpecies.data_type = 'INT'
        CompareSpecies.inputs[1].default_value = string2Number(spname)
        InstanceOnPoint = get_node_by_name(
            nodes,
            "InstanceOnPoint_%s_%s" % (self.label, spname),
            "GeometryNodeInstanceOnPoints",
        )
        ObjectInfo = get_node_by_name(
            nodes, "ObjectInfo_%s_%s" % (self.label, spname), "GeometryNodeObjectInfo"
        )
        ObjectInfo.inputs["Object"].default_value = self.batoms.species.instancers[
            spname
        ]
        BoolShow = get_node_by_name(
            nodes,
            "BooleanMath_%s_%s_1" % (self.label, spname),
            "FunctionNodeBooleanMath",
        )
        #
        links.new(SetPosition.outputs["Geometry"], InstanceOnPoint.inputs["Points"])
        links.new(GroupInput.outputs[2], CompareSpecies.inputs[0])
        links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        # transfer scale
        TransferScale = get_node_by_name(
            nodes, "%s_TransferScale" % (self.label), "GeometryNodeSampleIndex"
        )
        links.new(TransferScale.outputs[2], InstanceOnPoint.inputs["Scale"])
        links.new(CompareSpecies.outputs[0], BoolShow.inputs[1])
        links.new(BoolShow.outputs["Boolean"], InstanceOnPoint.inputs["Selection"])
        links.new(ObjectInfo.outputs["Geometry"], InstanceOnPoint.inputs["Instance"])
        links.new(InstanceOnPoint.outputs["Instances"], JoinGeometry.inputs["Geometry"])

    def update(self):
        """Main function to update boundary
        1) for each frame, search the boundary list
        2) mearch all boundary list, and find the unique boundary list
        3) calculate all boundary data
        """
        # object_mode()
        # clean_coll_objects(self.coll, 'bond')
        self.hide = False
        if not self.active:
            self.set_arrays(default_boundary_datas)
            return
        positions = self.batoms.positions
        trajectory = self.batoms.get_trajectory()["positions"]
        images = self.batoms.as_ase()
        nframe = self.batoms.nframe
        if nframe == 0:
            images = [images]
            trajectory = np.array([positions])
            nframe = 1
        tstart = time()
        for f in range(nframe):
            # print('update boundary: ', f)
            # print(images[f].positions)
            # use local positions for boundary search
            boundary_list = search_boundary(images[f], self.boundary)
            if f == 0:
                boundary_lists = boundary_list
            else:
                boundary_lists = np.append(boundary_lists, boundary_list, axis=0)
        boundary_lists = np.unique(boundary_lists, axis=0)
        boundary_datas = self.calc_boundary_data(
            boundary_lists, images[0].arrays, trajectory, self.batoms.cell
        )
        # update unit cell

        #
        self.set_arrays(boundary_datas)
        # self.batoms.draw()
        logger.debug("update boundary: {0:10.2f} s".format(time() - tstart))

    def update_geometry_node_instancer(self, spname, instancer):
        """When instances are re-build, we need also update
        the geometry node.

        Args:
            spname (str): name of the species
        """
        from batoms.utils.butils import get_node_by_name

        # update  instancers
        ObjectInfo = get_node_by_name(
            self.gn_node_group.nodes,
            "ObjectInfo_%s_%s" % (self.label, spname),
            "GeometryNodeObjectInfo",
        )
        ObjectInfo.inputs["Object"].default_value = instancer
        logger.debug("update boundary instancer: {}".format(spname))

    def update_gn_cell(self):
        from batoms.utils.butils import get_node_by_name

        # update cell
        cell = self.batoms.cell.array
        # set positions
        nodes = self.gn_node_group.nodes
        for i in range(3):
            tmp = get_node_by_name(
                nodes, "%s_VectorDot%s_%s" % (self.label, i, ""), "ShaderNodeVectorMath"
            )
            tmp.operation = "DOT_PRODUCT"
            tmp.inputs[1].default_value = cell[:, i]

    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            # , self.batoms.location)
            self.build_object(default_boundary_datas)
        return obj

    @property
    def active(self):
        return self.batoms.coll.batoms.boundary.active

    @active.setter
    def active(self, value):
        self.batoms.coll.batoms.boundary.active = value

    @property
    def boundary(self):
        return self.get_boundary()

    @boundary.setter
    def boundary(self, boundary):
        if boundary is not None:
            if isinstance(boundary, (int, float)):
                boundary = np.array([[-boundary, 1 + boundary]] * 3)
            elif len(boundary) == 3:
                if isinstance(boundary[0], (int, float)):
                    boundary = np.array(
                        [
                            [-boundary[0], 1 + boundary[0]],
                            [-boundary[1], 1 + boundary[1]],
                            [-boundary[2], 1 + boundary[2]],
                        ]
                    )
                elif len(boundary[0]) == 2:
                    boundary = np.array(boundary)
            else:
                raise Exception("Wrong boundary setting!")
            if np.isclose(boundary[:].flatten(), np.array([0, 1, 0, 1, 0, 1])).all():
                self.active = False
            else:
                self.active = True
            self.batoms.coll.batoms.boundary.boundary = boundary[:].flatten()
        else:
            self.active = False
        self.update()
        # self.batoms.draw()

    def get_boundary(self):
        boundary = np.array(self.batoms.coll.batoms.boundary.boundary)
        return boundary.reshape(3, -1)

    def set_arrays(self, arrays):
        """ """
        attributes = self.attributes
        # same length
        dnvert = len(arrays["species_index"]) - len(attributes["species_index"])
        if dnvert > 0:
            self.add_vertices_bmesh(dnvert)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
        if len(arrays["positions"]) == 0:
            return
        self.positions = arrays["positions"][0]
        self.set_trajectory(arrays)
        self.update_mesh()
        species_index = [string2Number(sp) for sp in arrays["species"]]
        self.set_attributes(
            {
                "atoms_index": arrays["atoms_index"],
                "species_index": species_index,
                "scale": arrays["scales"],
                "show": arrays["shows"],
                "boundary_offset": arrays["boundary_offset"],
            }
        )
        species = np.unique(arrays["species"])
        for sp in species:
            self.add_geometry_node(sp)

    def get_arrays(self):
        """ """
        # object_mode()
        # tstart = time()
        arrays = self.attributes
        arrays.update({"positions": self.positions})
        # radius
        radius = self.batoms.radius
        arrays.update({"radius": np.zeros(len(self))})
        species = np.array(
            [number2String(i) for i in arrays["species_index"]], dtype="U20"
        )
        arrays["species"] = species
        for sp, value in radius.items():
            mask = np.where(arrays["species"] == sp)
            arrays["radius"][mask] = value
        # size
        arrays["size"] = arrays["radius"] * arrays["scale"]
        # main elements
        main_elements = self.batoms.species.main_elements
        elements = [main_elements[sp] for sp in arrays["species"]]
        arrays.update({"elements": np.array(elements, dtype="U20")})
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    @property
    def boundary_data(self):
        return self.get_boundary_data()

    def get_boundary_data(self, include_batoms=False):
        """ """
        # check cell voluem
        # if < 1e-6, invalid cell and boundary.
        if self.batoms.cell.volume < 1e-6:
            return None
        arrays = self.arrays
        boundary_data = {
            "positions": arrays["positions"],
            "species": arrays["species"],
            "indices": arrays["atoms_index"],
            "boundary_offset": arrays["boundary_offset"],
        }
        return boundary_data

    @property
    def bondlists(self):
        return self.get_bondlists()

    def get_bondlists(self):
        bondlists = self.batoms.bonds.arrays
        return bondlists

    def get_trajectory(self):
        """ """
        trajectory = {}
        trajectory["positions"] = self.get_obj_trajectory(self.obj)
        return trajectory

    def set_trajectory(self, trajectory=None, frame_start=0):
        if trajectory is None:
            trajectory = self._trajectory
        name = "%s_boundary" % (self.label)
        obj = self.obj
        self.set_shape_key(name, obj, trajectory["positions"], frame_start=frame_start)

    def calc_boundary_data(self, boundary_lists, arrays, positions, cell):
        """ """
        # tstart = time()
        # properties
        model_styles = arrays["model_style"][boundary_lists[:, 0]]
        shows = arrays["show"][boundary_lists[:, 0]]
        selects = arrays["select"][boundary_lists[:, 0]]
        scales = arrays["scale"][boundary_lists[:, 0]]
        species_indexs = arrays["species_index"][boundary_lists[:, 0]]
        species = arrays["species"][boundary_lists[:, 0]]
        # ------------------------------------
        offset_vectors = boundary_lists[:, 1:4]
        offsets = np.dot(offset_vectors, cell)
        # use global positions to build boundary positions
        if len(positions.shape) == 2:
            positions = np.array([positions])
        positions = positions[:, boundary_lists[:, 0]] + offsets
        datas = {
            "atoms_index": np.array(boundary_lists[:, 0]),
            "species_index": species_indexs,
            "species": species,
            "positions": positions,
            # 'offsets':offsets,
            "boundary_offset": offset_vectors,
            "model_styles": model_styles,
            "shows": shows,
            "selects": selects,
            "scales": scales,
        }
        # print('datas: ', datas)
        # print('calc_boundary_data: {0:10.2f} s'.format(time() - tstart))
        return datas

    def __getitem__(self, index):
        return self.boundary[index]

    def __setitem__(self, index, value):
        """Set boundary vectors."""
        boundary = self.boundary
        if isinstance(value, (int, float)):
            boundary = np.array([-boundary, 1 + boundary])
        boundary[index] = value
        self.boundary = boundary

    def __repr__(self) -> str:
        s = self.boundary.__repr__()
        return s

    def as_dict(self):
        """ """
        data = {"array": None}
        if len(self) > 1:
            data["array"] = dict(self.arrays)
        data["boundary"] = self.boundary
        return data


def search_boundary(
    atoms,
    boundary=[[0, 1], [0, 1], [0, 1]],
):
    """Search atoms in the boundary

    Args:
        atoms: _description_
        boundary (list, optional): _description_.
            Defaults to [[0, 1], [0, 1], [0, 1]].

    Returns:
        _type_: _description_
    """
    # tstart = time()
    cell = atoms.cell
    positions = atoms.positions
    species = atoms.arrays["species"]
    if isinstance(boundary, float):
        boundary = [
            [-boundary, 1 + boundary],
            [-boundary, 1 + boundary],
            [-boundary, 1 + boundary],
        ]
    boundary = np.array(boundary)
    # find supercell
    f = np.floor(boundary)
    c = np.ceil(boundary)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    M = np.product(ib[1] - ib[0] + 1)
    # get scaled positions
    positions = np.linalg.solve(complete_cell(cell).T, positions.T).T
    n = len(positions)
    npositions = np.tile(positions, (M - 1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    # index
    offsets = np.zeros((M * n, 4), dtype=int)
    ind0 = np.arange(n).reshape(-1, 1)
    species0 = species
    species = []
    # repeat the positions so that
    # it completely covers the boundary
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0:
                    continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                offsets[i0:i1] = np.append(ind0, np.array([[m0, m1, m2]] * n), axis=1)
                species.extend(species0)
                i0 = i1
    # boundary condition
    ind1 = np.where(
        (npositions[:, 0] > boundary[0][0])
        & (npositions[:, 0] < boundary[0][1])  # noqa: W503
        & (npositions[:, 1] > boundary[1][0])  # noqa: W503
        & (npositions[:, 1] < boundary[1][1])  # noqa: W503
        & (npositions[:, 2] > boundary[2][0])  # noqa: W503
        & (npositions[:, 2] < boundary[2][1])  # noqa: W503
    )[0]
    offsets_b = offsets[ind1]
    # print('search boundary: {0:10.2f} s'.format(time() - tstart))
    return offsets_b
