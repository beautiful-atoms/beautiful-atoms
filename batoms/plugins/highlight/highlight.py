"""Definition of the highlight class.

This module defines the highlight object in the Batoms package.

"""

import bpy
import numpy as np
from batoms.base.object import ObjectGN
from batoms.plugins.base import PluginObject
from .setting import HighlightSettings
from batoms.utils.butils import get_node_by_name
from batoms.utils import string2Number
import logging

logger = logging.getLogger(__name__)

default_attributes = [
    {"name": "atom_index", "data_type": "INT"},
    {"name": "select_index", "data_type": "INT"},
    {"name": "select", "data_type": "STRING"},
    {"name": "show", "data_type": "INT"},
    {"name": "scale", "data_type": "FLOAT"},
]

default_GroupInput = [
    ["atom_index", "NodeSocketInt"],
    ["select_index", "NodeSocketInt"],
    ["show", "NodeSocketBool"],
    ["scale", "NodeSocketFloat"],
]

default_highlight_datas = {
    "select_index": np.ones(0, dtype=int),
    "atom_index": np.ones(0, dtype=int),
    "centers": np.zeros((1, 0, 3)),
    "show": np.zeros(0, dtype=int),
    "scale": np.ones(0, dtype=float),
}


class Highlight(ObjectGN, PluginObject):
    def __init__(
        self,
        label=None,
        highlight_datas=None,
        location=np.array([0, 0, 0]),
        batoms=None,
    ):
        """Highlight Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = "highlight"
        ObjectGN.__init__(self, label, name)
        PluginObject.__init__(self, label, "highlight", batoms)
        flag = self.load()
        if not flag:
            if highlight_datas is None:
                highlight_datas = default_highlight_datas
            self.build_object(highlight_datas)
            self.settings = HighlightSettings(self.label, batoms=batoms, parent=self)
        else:
            self.settings = HighlightSettings(self.label, batoms=batoms, parent=self)
        self.settings.bpy_data.active = True

    def build_materials(self, name, color, node_inputs=None, material_style="default"):
        """ """
        from batoms.material import create_material

        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(
            name,
            color=color,
            node_inputs=node_inputs,
            material_style=material_style,
            backface_culling=False,
        )
        return mat

    def draw(self):
        """Draw highlight."""
        # delete old highlight
        highlight_datas = self.build_highlight()
        self.set_arrays(highlight_datas)

    def build_highlight(self):
        """Find cage base on radius range and boundary range
        Algorith:
            1) build a meshgrid
            2) calulated first neighbour distance by kdtree
            3) find the max distance (max radius sphere)
            4) remove meshgrid within this sphere,
            5) find next max radius, ..., untill radius < minCave
            6) remove all spheres contact with boundary
        Args:
            cell (_type_): _description_
            positions (_type_): _description_
            radius (_type_): _description_
            step (float, optional): _description_. Defaults to 1.0.
            boundary (list, optional): _description_. Defaults to [[0, 1], [0, 1], [0, 1]].

        Returns:
            _type_: _description_
        """
        arrays = self.batoms.arrays
        natom = 0
        for hl in self.settings.bpy_setting:
            indices = self.batoms.selects[hl.select].indices
            natom += len(indices)
        highlights = {
            "select_index": np.ones(natom, dtype=int),
            "atom_index": np.ones(natom, dtype=int),
            "centers": np.zeros((1, natom, 3)),
            "show": np.zeros(natom, dtype=int),
            "scale": np.ones(natom, dtype=float),
        }
        i = 0
        for hl in self.settings.bpy_setting:
            indices = self.batoms.selects[hl.select].indices
            n = len(indices)
            highlights["centers"][0][i : i + n] = arrays["positions"][indices]
            highlights["scale"][i : i + n] = arrays["scale"][indices] * hl.scale
            highlights["show"][i : i + n] = np.ones(len(indices))
            highlights["select_index"][i : i + n] = np.ones(
                len(indices)
            ) * string2Number(hl.name)
            highlights["atom_index"][i : i + n] = np.array(indices)
            i += n
        return highlights

    def build_object(self, highlight_datas, attributes={}):
        """Build the main Batoms object

        Args:
            label (str):
                Name of the object
            arrays (array):
                arrays of properties for each atoms
            location (list, optional):
                Location of the object. Defaults to [0, 0, 0].
        """
        if len(highlight_datas["centers"].shape) == 2:
            self._trajectory = {
                "centers": np.array([highlight_datas["centers"]]),
            }
            centers = highlight_datas["centers"]
        elif len(highlight_datas["centers"].shape) == 3:
            self._trajectory = {
                "centers": highlight_datas["centers"],
            }
            centers = highlight_datas["centers"][0]
        else:
            raise Exception("Shape of centers is wrong!")
        attributes.update(
            {
                "select_index": highlight_datas["select_index"],
                "show": highlight_datas["show"],
                "scale": highlight_datas["scale"],
            }
        )
        name = self.obj_name
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        obj = bpy.data.objects.new(name, mesh)
        obj.data.from_pydata(centers, [], [])
        obj.batoms.type = "HIGHLIGHT"
        obj.batoms.label = self.batoms.label
        self.batoms.coll.objects.link(obj)
        # Add attributes
        for att in default_attributes:
            self.add_attribute(**att)
        # add cell object as its child
        obj.parent = self.batoms.obj
        self.set_attributes(attributes)
        self.init_geometry_node_modifier(default_GroupInput)
        self.build_geometry_node()

    def build_geometry_node(self):
        """Geometry node for instancing sphere on vertices!"""
        nodes = self.gn_node_group.nodes
        links = self.gn_node_group.links
        GroupInput = nodes[0]
        GroupOutput = nodes[1]
        # print(gn.name)
        JoinGeometry = get_node_by_name(
            nodes, "%s_JoinGeometry" % self.label, "GeometryNodeJoinGeometry"
        )
        SeparateGeometry = get_node_by_name(
            nodes, "%s_SeparateGeometry" % self.label, "GeometryNodeSeparateGeometry"
        )
        links.new(GroupInput.outputs["Geometry"], SeparateGeometry.inputs["Geometry"])
        links.new(GroupInput.outputs[2], SeparateGeometry.inputs["Selection"])
        links.new(SeparateGeometry.outputs[0], JoinGeometry.inputs["Geometry"])
        links.new(JoinGeometry.outputs["Geometry"], GroupOutput.inputs["Geometry"])
        #
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
        #
        # set positions
        SetPosition = get_node_by_name(
            nodes, "%s_SetPosition" % self.label, "GeometryNodeSetPosition"
        )
        links.new(GroupInput.outputs["Geometry"], SetPosition.inputs["Geometry"])
        links.new(TransferBatoms.outputs[0], SetPosition.inputs["Position"])

    def add_geometry_node(self, slname, instancer):
        """Add geometry node for each select

        Args:
            slname (str):
                Name of the select
            instancer (bpy.data.object):
                Object to be instanced
        """
        from batoms.utils.butils import compareNodeType

        nodes = self.gn_node_group.nodes
        links = self.gn_node_group.links
        GroupInput = nodes[0]
        SetPosition = get_node_by_name(nodes, "%s_SetPosition" % self.label)
        JoinGeometry = get_node_by_name(
            nodes, "%s_JoinGeometry" % self.label, "GeometryNodeJoinGeometry"
        )
        CompareSelect = get_node_by_name(
            nodes, "CompareFloats_%s_%s" % (self.label, slname), compareNodeType
        )
        CompareSelect.operation = "EQUAL"
        # CompareSelect.data_type = 'INT'
        CompareSelect.inputs[1].default_value = string2Number(slname)
        InstanceOnPoint = get_node_by_name(
            nodes,
            "InstanceOnPoint_%s_%s" % (self.label, slname),
            "GeometryNodeInstanceOnPoints",
        )
        ObjectInfo = get_node_by_name(
            nodes, "ObjectInfo_%s_%s" % (self.label, slname), "GeometryNodeObjectInfo"
        )
        ObjectInfo.inputs["Object"].default_value = instancer
        BoolShow = get_node_by_name(
            nodes,
            "BooleanMath_%s_%s_1" % (self.label, slname),
            "FunctionNodeBooleanMath",
        )
        #
        links.new(SetPosition.outputs["Geometry"], InstanceOnPoint.inputs["Points"])
        links.new(GroupInput.outputs[2], CompareSelect.inputs[0])
        links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        links.new(GroupInput.outputs[4], InstanceOnPoint.inputs["Scale"])
        links.new(CompareSelect.outputs[0], BoolShow.inputs[1])
        links.new(BoolShow.outputs["Boolean"], InstanceOnPoint.inputs["Selection"])
        links.new(ObjectInfo.outputs["Geometry"], InstanceOnPoint.inputs["Instance"])
        links.new(InstanceOnPoint.outputs["Instances"], JoinGeometry.inputs["Geometry"])

    def update_geometry_node_instancer(self):
        """ """
        for sp in self.settings.bpy_setting:
            self.settings.build_instancer(sp.as_dict())

    def set_arrays(self, arrays):
        """ """
        # if len(arrays['positions']) == 0:
        #     return
        attributes = self.attributes
        print(arrays)
        # same length
        dnvert = len(arrays["select_index"]) - len(attributes["select_index"])
        if dnvert > 0:
            self.add_vertices_bmesh(dnvert)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
        self.set_trajectory(arrays)
        self.set_attributes({"select_index": arrays["select_index"]})
        self.set_attributes({"scale": arrays["scale"]})
        self.set_attributes({"show": arrays["show"]})
        self.set_attributes({"atom_index": arrays["atom_index"]})
        self.update_mesh()
        self.update_geometry_node_instancer()

    def set_trajectory(self, frames=None, frame_start=0):
        if frames is None:
            frames = self._trajectory
        nframe = len(frames["centers"])
        if nframe == 0:
            return
        name = "%s_highlight" % (self.label)
        obj = self.obj
        self.set_shape_key(name, obj, frames["centers"], frame_start=frame_start)

    @property
    def objs(self):
        objs = {}
        name = "{}_{}".format(self.label, self.name)
        obj = bpy.data.objects.get(name)
        objs[name] = obj
        return objs

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
        data = {}
        data["settings"] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data
