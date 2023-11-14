import bpy
import numpy as np
from batoms.utils.butils import (
    get_node_by_name,
    object_mode,
    set_look_at,
    update_object,
)

import bmesh
import logging

# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


default_object_attributes = []

default_object_datas = {}


class BaseObject:
    def __init__(self, obj_name, btype="batoms"):
        self.obj_name = obj_name
        self.btype = btype

    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError("%s object is not exist." % self.obj_name)
        return obj

    @property
    def bpy_data(self):
        return self.get_bpy_data()

    def get_bpy_data(self):
        bpy_data = getattr(self.obj, self.btype)
        return bpy_data

    @property
    def location(self):
        return self.get_location()

    def get_location(self):
        return np.array(self.obj.location)

    @location.setter
    def location(self, location):
        self.set_location(location)

    def set_location(self, location):
        self.obj.location = location

    @property
    def type(self):
        return self.get_type()

    # this is weak, because not work for mesh object

    @type.setter
    def type(self, type):
        self.set_type(type)

    def get_type(self):
        return self.obj.data.type

    def set_type(self, type):
        self.obj.data.type = type.upper()

    @property
    def hide(self):
        return self.get_hide()

    @hide.setter
    def hide(self, state):
        self.set_hide(state)

    def get_hide(self):
        return self.obj.hide_get()

    def set_hide(self, state):
        self.obj.hide_render = state
        self.obj.hide_set(state)

    @property
    def select(self):
        return self.get_select()

    @select.setter
    def select(self, state):
        self.set_select(state)

    def get_select(self):
        return self.obj.select_get()

    def set_select(self, state):
        self.obj.select_set(state)

    @property
    def scene(self):
        return self.get_scene()

    def get_scene(self):
        return bpy.context.scene

    @property
    def look_at(self):
        return self.get_look_at()

    @look_at.setter
    def look_at(self, look_at):
        self.set_look_at(look_at)

    def get_look_at(self):
        return np.array(self.bpy_data.look_at)

    def set_look_at(self, look_at):
        self.bpy_data.look_at = look_at
        set_look_at(self.obj, look_at, roll=0.0)

    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h.translate([0, 0, 5])
        """
        object_mode()
        bpy.ops.object.select_all(action="DESELECT")
        self.obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)

    def rotate(self, angle, axis="Z", orient_type="GLOBAL"):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:

        >>> h.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action="DESELECT")
        self.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.transform.rotate(
            value=angle, orient_axis=axis.upper(), orient_type=orient_type
        )

    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink=True)

    def delete_material(self, name):
        if name in bpy.data.materials:
            obj = bpy.data.materials.get(name)
            bpy.data.materials.remove(obj, do_unlink=True)

    def update_mesh(self, obj=None):
        """Update mesh of the object.
        By switching to edit mode and back to object mode."""

        object_mode()
        if obj is None:
            obj = self.obj
        # bm = bmesh.new()
        # bm.from_mesh(obj.data)
        # bm.verts.ensure_lookup_table()
        # bm.to_mesh(obj.data)
        # bm.clear()
        # bpy.context.view_layer.update()
        # I don't why the shape key is not udpate in background mode, so we need this.
        bpy.context.view_layer.objects.active = self.obj
        mode = self.obj.mode
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.object.mode_set(mode="OBJECT")
        bpy.ops.object.mode_set(mode=mode)

    def add_verts(self, count, obj=None):
        """Add vertices to the object."""

        object_mode()
        if obj is None:
            obj = self.obj
        obj.data.vertices.add(count)
        self.update_mesh(obj)

    def add_vertices_bmesh(self, count, obj=None):
        import bmesh
        import numpy as np

        object_mode()
        if obj is None:
            obj = self.obj
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        vertices = np.zeros((count, 3))
        for vert in vertices:
            bm.verts.new(vert)
        bm.to_mesh(obj.data)
        bm.clear()

    def delete_vertices_bmesh(
        self,
        index=[],
        obj=None,
    ):
        """
        delete verts
        """
        import bmesh

        object_mode()
        if obj is None:
            obj = self.obj
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        verts_select = [bm.verts[i] for i in index]
        bmesh.ops.delete(bm, geom=verts_select, context="VERTS")
        bm.to_mesh(obj.data)
        bm.clear()

    def add_edges(self, edges=None, obj=None):
        object_mode()
        if obj is None:
            obj = self.obj
        obj.data.edges.add(len(edges))
        edges = edges.reshape(len(edges) * 2)
        obj.data.edges.foreach_set("vertices", edges)
        self.update_mesh(obj)

    def delete_edges_bmesh(
        self,
        index=[],
        obj=None,
    ):
        """
        delete edges
        """
        import bmesh

        object_mode()
        if obj is None:
            obj = self.obj
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.edges.ensure_lookup_table()
        edges_select = [bm.edges[i] for i in index]
        bmesh.ops.delete(bm, geom=edges_select, context="EDGES_FACES")
        bm.to_mesh(obj.data)
        bm.clear()

    @property
    def shape_keys(self):
        base_name = "Basis_%s" % self.obj_name
        if self.obj.data.shape_keys is None and len(self) > 0:
            self.obj.shape_key_add(name=base_name)
        return self.obj.data.shape_keys

    def __len__(self):
        return len(self.obj.data.vertices)


class ObjectGN(BaseObject):
    """
    Object with Geometry Node.

    """

    def __init__(self, label, name=None):
        if name:
            self.name = name
            obj_name = "%s_%s" % (label, name)
        else:
            obj_name = label
        BaseObject.__init__(self, obj_name=obj_name)

    def build_object(self, arrays, attributes={}):
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_shape_key(self._trajectory)

    def load(self):
        flag = True
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            flag = False
        return flag

    @property
    def gn_modifier(self):
        return self.get_gn_modifier()

    @property
    def gn_node_group(self):
        """Get the geometry node group of the object."""
        return self.gn_modifier.node_group

    def get_gn_modifier(self):
        """Get the geometry node modifier of the object."""
        name = "GeometryNodes_%s" % self.obj_name
        modifier = self.obj.modifiers.get(name)
        if modifier is None:
            self.init_geometry_node_modifier()
        return modifier

    def init_geometry_node_modifier(self, inputs=[]):
        """Init geometry node modifier"""
        # blender 4.0 use a interface to add sockets, both input and output
        from batoms.utils.butils import build_gn_modifier

        name = "GeometryNodes_%s" % self.obj_name
        modifier = build_gn_modifier(self.obj, name)
        if bpy.app.version_string >= "4.0.0":
            interface = modifier.node_group.interface
            for input in inputs:
                socket = interface.new_socket(
                    name=input[0], socket_type=input[1], in_out="INPUT"
                )
                modifier["%s_use_attribute" % socket.identifier] = True
                modifier["%s_attribute_name" % socket.identifier] = input[0]
        else:
            for input in inputs:
                modifier.node_group.inputs.new(input[1], input[0])
                id = modifier.node_group.inputs[input[0]].identifier
                modifier["%s_use_attribute" % id] = True
                modifier["%s_attribute_name" % id] = input[0]
        return modifier

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        from batoms.utils.butils import build_gn_modifier

        name = "GeometryNodes_%s" % self.obj_name
        modifier = build_gn_modifier(self.obj, name)
        return modifier

    def vectorDotMatrix(self, gn_node_group, vector_output, matrix, name):
        """ """
        CombineXYZ = get_node_by_name(
            gn_node_group.nodes,
            "%s_CombineXYZ_%s" % (self.label, name),
            "ShaderNodeCombineXYZ",
        )
        #
        VectorDot = []
        for i in range(3):
            tmp = get_node_by_name(
                gn_node_group.nodes,
                "%s_VectorDot%s_%s" % (self.label, i, name),
                "ShaderNodeVectorMath",
            )
            tmp.operation = "DOT_PRODUCT"
            VectorDot.append(tmp)
            tmp.inputs[1].default_value = matrix[:, i]
            gn_node_group.links.new(vector_output, tmp.inputs[0])
            gn_node_group.links.new(tmp.outputs["Value"], CombineXYZ.inputs[i])
        return CombineXYZ

    def add_geometry_node(self):
        """
        add geometry node for species
        """
        return NotImplementedError("add_geometry_node is not implemented.")

    @property
    def realize_instances(self):
        return self.get_realize_instances()

    @realize_instances.setter
    def realize_instances(self, state):
        self.set_realize_instances(state)

    def get_realize_instances(self):
        return list(self.coll.batoms.realize_instances)

    def set_realize_instances(self, realize_instances):
        """Make instancing object real
        # TODO: add make real to geometry node
        """
        #
        nodes = self.gn_node_group.nodes
        RealizeInstances = get_node_by_name(
            self.gn_node_group.nodes,
            "%s_RealizeInstances" % self.label,
            "GeometryNodeRealizeInstances",
        )
        if not realize_instances:
            # switch off
            if len(RealizeInstances.outputs[0].links) > 0:
                link = RealizeInstances.outputs[0].links[0]
                self.gn_node_group.links.remove(link)
            self.gn_node_group.links.new(
                nodes["%s_JoinGeometry" % self.label].outputs[0], nodes[1].inputs[0]
            )
        else:
            self.gn_node_group.links.new(
                nodes["%s_JoinGeometry" % self.label].outputs[0],
                RealizeInstances.inputs[0],
            )
            self.gn_node_group.links.new(
                RealizeInstances.outputs[0], nodes[1].inputs[0]
            )
        self.gn_node_group.update_tag()

    def update(
        self,
    ):
        """
        update object'data.
        calculate data in all farmes
        """
        pass

    @property
    def location(self):
        return self.get_location()

    def get_location(self):
        return np.array(self.obj.location)

    @location.setter
    def location(self, location):
        self.set_location(location)

    def set_location(self, location):
        self.obj.location = location

    @property
    def arrays(self):
        return self.get_arrays()

    @arrays.setter
    def arrays(self, arrays):
        self.set_arrays(arrays)

    def get_arrays(self):
        """Get the arrays of the Batoms object.

        arrays is the combination of attributes and the  positions


        Returns:
            dict: _description_
        """
        object_mode()
        arrays = self.attributes
        arrays.update(
            {
                "positions": self.positions[0],
            }
        )
        return arrays

    @property
    def attributes(self):
        return self.get_attributes()

    @attributes.setter
    def attributes(self, attributes):
        self.set_attributes(attributes)

    def get_attributes(self, domains=["POINT"]):
        """Get all attributes of the Batoms object.


        Returns:
            dict: attributes dict
        """
        attributes = {}
        for att in self.obj.data.attributes:
            if att.name.startswith(".") or att.domain not in domains:
                continue
            attributes[att.name] = self.get_attribute(att.name)
        return attributes

    def get_attribute(self, key, index=None):
        """Get an attribute of Batoms object by name
        If the attribute is single value, just use the attribute of mesh.
        If the attribute is an array, we should gather several attributes
        from the mesh.

        Args:
            key (string): name of the attribute

        Raises:
            KeyError: _description_

        Returns:
            _type_: _description_
        """
        from batoms.utils.attribute import get_mesh_attribute, get_mesh_attribute_bmesh

        # get the mesh
        obj = self.obj
        att = obj.data.attributes[key]
        if obj.mode == "EDIT" and att.data_type in [
            "STRING",
            "INT",
            "FLOAT",
            "FLOAT_VECTOR",
            "FLOAT_COLOR",
        ]:
            attribute = get_mesh_attribute_bmesh(obj, att.name, index)
        else:
            attribute = get_mesh_attribute(obj, att.name, index)

        return attribute

    def add_attribute_from_array(self, name, data, domain="POINT"):
        """Add an attribute from array"""
        from batoms.utils import type_py_to_blender

        # check the type of the array, and convert it to blender type
        # only the first element is checked
        # check if shape of the data array
        # if 1x2, 1x3, 1x4
        if isinstance(data, list):
            data = np.array(data)
        try:
            if len(data.shape) == 1:
                array = np.asarray(data)
                type_py = type(array.flat[0])
                dtype_bl = type_py_to_blender(type_py)
                if dtype_bl is False:
                    print(
                        f"Attribute: {name} is not added. The type of the array: {type_py} is not supported."  # noqa E501
                    )
                    return False
                self.add_attribute(name=name, data_type=dtype_bl, domain=domain)
            elif len(data.shape) == 2:
                if data.shape[1] == 2:
                    self.add_attribute(name=name, data_type="FLOAT2", domain=domain)
                elif data.shape[1] == 3:
                    self.add_attribute(
                        name=name, data_type="FLOAT_VECTOR", domain=domain
                    )
                elif data.shape[1] == 4:
                    self.add_attribute(name=name, data_type="QUATERNION", domain=domain)
                else:
                    print(
                        f"Attribute: {name} is not added. The shape of the array: {data.shape} is not supported."  # noqa E501
                    )
                    return False
            else:
                # print a Warning message
                print(
                    f"Attribute: {name} is not added. The shape of the array: {data.shape} is not supported."  # noqa E501
                )
                return False
        except Exception:
            print(
                f"Attribute: {name} is not added. The shape of the array: {data.shape} is not supported."  # noqa E501
            )
            return False

    def add_attribute(self, name, data_type="FLOAT", domain="POINT"):
        """Add an attribute to the mesh"""
        self.obj.data.attributes.new(name=name, type=data_type, domain=domain)
        return True

    def set_attributes(self, attributes):
        """Set attributes

        Args:
            attributes (dict): attributes bound to every atoms
        """
        obj = self.obj
        me = obj.data
        # print(attributes)
        for key, array in attributes.items():
            # if key in self._attributes
            if key not in obj.data.attributes:
                # try to create a new attribute
                flag = self.add_attribute_from_array(key, array)
                logger.debug(f"Add new attribute: {key}")
                # if failed, do not add this attribute
                if flag is False:
                    continue
            self.set_attribute(key, array)
        me.update()

    def set_attribute(self, key, array, index=None):
        """Set one attribute

        An attribute is a generic term to describe data stored per-element
        in a geometry data-block. For example, every vertex can have an
        associated number or vector.

        class bpy.types.Attribute

        Two cases:
        1) for single value data (dimension=0), treat it as normal
        2) for array data (dimension > 0), scatter the data.

        Two mode:
        1) Edit mode, use bmesh
        2) Object mode, use data.attributes directly

        Special case:
        1) String, can not use foreach_set. Must use bmesh with encode,
            otherwise, can not read use bmesh.
        2) Boolean, does not supported by bmesh. Must use Object mode.
           We set all Boolean properties to INT.

        Args:
            key (str): name of the attribute
            array (np.array): value of the attribute
        """
        from batoms.utils.attribute import set_mesh_attribute, set_mesh_attribute_bmesh

        obj = self.obj
        me = obj.data
        if key not in me.attributes:
            raise KeyError(
                "Attribute: {} is not exist. Please add it first.".format(key)
            )
        array = np.array(array)
        if len(array.shape) == 0:
            array = np.array([array])
        if len(array) == 0:
            return
        # single value data
        att = me.attributes.get(key)
        if att.data_type == "STRING" or (
            obj.mode == "EDIT" and att.data_type in ["INT", "FLOAT"]
        ):
            set_mesh_attribute_bmesh(obj, key, array, index)
        else:
            set_mesh_attribute(obj, key, array, index)

    def set_attribute_with_indices(self, name, indices, data):
        data0 = self.get_attribute(name)
        data0[indices] = data
        self.set_attributes({name: data0})

    def get_attribute_with_indices(self, name, indices):
        """Get attribute with indices, with len(indices) > 1

        Read the whole attribute, and then slice.

        If len(indices) == 1, use get_attribute(key, index=indices[0]),
        this is faster.

        Args:
            name (str): _description_
            indices (array): _description_

        Returns:
            _type_: _description_
        """
        return self.get_attribute(name)[indices]

    @property
    def positions(self):
        return self.get_positions()

    def get_positions(self):
        """Get the local positions of the object."""
        n = len(self)
        positions = np.empty(n * 3, dtype=np.float64)
        self.obj.data.vertices.foreach_get("co", positions)
        positions = positions.reshape((n, 3))
        return positions

    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)

    def set_positions(self, positions):
        """
        Set positions to local vertices
        """
        object_mode()
        natom = len(self)
        n = len(positions)
        if n != natom:
            raise ValueError("positions has wrong shape %s != %s." % (n, natom))
        if natom == 0:
            return
        positions = positions.reshape((natom * 3, 1))
        self.obj.data.vertices.foreach_set("co", positions)
        self.obj.data.update()
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.object.mode_set(mode="OBJECT")

    @property
    def global_positions(self):
        return self.get_global_positions()

    @global_positions.setter
    def global_positions(self, global_positions):
        self.set_global_positions(global_positions)

    def get_global_positions(self):
        """
        Get global positions.
        """
        from batoms.utils import local2global

        global_positions = local2global(self.positions, np.array(self.obj.matrix_world))
        return global_positions

    def set_global_positions(self, positions):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.utils import local2global

        natom = len(self)
        nposition = len(positions)
        if nposition != natom:
            raise ValueError("positions has wrong shape %s != %s." % (nposition, natom))
        if natom == 0:
            return
        local_positions = local2global(
            positions, np.array(self.obj.matrix_world), reversed=True
        )
        self.local_positions = local_positions

    @property
    def nframe(self):
        return self.get_nframe()

    def get_nframe(self):
        if self.obj.data.shape_keys is None:
            return 0
        nframe = len(self.obj.data.shape_keys.key_blocks)
        return nframe

    def get_shape_key(self, obj, local=True):
        """
        read shape key
        """
        from batoms.utils import local2global

        n = len(self)
        nframe = self.nframe
        trajectory = np.empty((nframe, n, 3), dtype=np.float64)
        for i in range(nframe):
            positions = np.empty(n * 3, dtype=np.float64)
            sk = obj.data.shape_keys.key_blocks[i]
            sk.data.foreach_get("co", positions)
            local_positions = positions.reshape((n, 3))
            if local:
                trajectory[i] = local_positions
            else:
                global_positions = local2global(
                    local_positions, np.array(self.obj.matrix_world)
                )
                trajectory[i] = global_positions
        return trajectory

    def set_shape_key(self, name, obj, trajectory=None, frame_start=0):
        """Set shape key for trajectory

        name: str
            name of the shape key
        obj: bpy.data.objects
            object to set shape key
        trajectory: list
            list of positions
        frame_start: int
            start frame
        """
        from batoms.utils.butils import add_keyframe_to_shape_key

        if trajectory is None:
            trajectory = self.trajectory
        nframe = len(trajectory)
        if nframe == 0:
            return
        # shape key should not be add for empty mesh
        if len(trajectory[0]) == 0:
            return
        # name = '%s_bond%s'%(self.label, sp)
        # obj = bpy.data.objects.get(name)
        base_name = "Basis_%s" % (name)
        if obj.data.shape_keys is None:
            sk = obj.shape_key_add(name=base_name)
        elif base_name not in obj.data.shape_keys.key_blocks:
            sk = obj.shape_key_add(name=base_name)
        else:
            sk = obj.data.shape_keys.key_blocks.get(base_name)
        # set basis key
        nvert = len(obj.data.shape_keys.key_blocks[0].data)
        positions = trajectory[0]
        vertices = positions.reshape(nvert * 3)
        sk.data.foreach_set("co", vertices)
        # self.obj.data.update()
        for i in range(1, nframe):
            name = str(i)
            if name not in obj.data.shape_keys.key_blocks:
                sk = obj.shape_key_add(name=name)
                # Add Keytrajectory, the last one is different
                if i != nframe - 1:
                    add_keyframe_to_shape_key(
                        sk,
                        "value",
                        [0, 1, 0],
                        [frame_start + i - 1, frame_start + i, frame_start + i + 1],
                    )
                else:
                    add_keyframe_to_shape_key(
                        sk, "value", [0, 1], [frame_start + i - 1, frame_start + i]
                    )
            else:
                sk = obj.data.shape_keys.key_blocks.get(name)
            # Use the local position here
            positions = trajectory[i]
            vertices = positions.reshape((nvert * 3, 1))
            sk.data.foreach_set("co", vertices)
        self.update_mesh(obj)

    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            # if not obj.batoms.flag:
            # raise Exception(
            # "Failed, the name {} already in use and is not
            # Batom object!".format(name))
            bpy.data.objects.remove(obj, do_unlink=True)


class childObjectGN:
    """
    Slice of Object with Geometry Node.

    """

    def __init__(self, label, indices, obj_name=None, parent=None):
        self.label = label
        self.parent = parent
        if isinstance(indices, slice):
            self.indices = np.arange(len(self.parent))[indices]
        elif isinstance(indices, int):
            self.indices = np.array([indices])
        else:
            self.indices = np.array(indices)

        if obj_name is None:
            self.obj_name = label
        else:
            self.obj_name = obj_name

    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError("%s object is not exist." % self.obj_name)
        return obj

    @property
    def vertice(self):
        if self.obj.data.shape_keys is None:
            return self.obj.data.vertices[self.indices[0]]
        else:
            return self.obj.data.shape_keys.key_blocks[0].data[self.indices[0]]

    @property
    def bm(self):
        if self.obj.mode == "EDIT":
            bm = bmesh.from_edit_mesh(self.obj.data)
        else:
            bm = bmesh.new()
            bm.from_mesh(self.obj.data)
        return bm

    @property
    def attributes(self):
        return self.get_attributes()

    def get_attributes(self):
        """Get all attributes of the Batoms object.

        Returns:
            dict: attributes dict
        """
        attributes = {}
        for att in self.parent._attributes:
            attributes[att.name] = self.get_attribute(att.name)
        return attributes

    @property
    def local_position(self):
        if len(self.indices) == 1:
            return self.vertice.co
        else:
            return self.parent.local_positions[self.indices]

    @property
    def position(self):
        if len(self.indices) == 1:
            return self.get_position()
        else:
            return self.parent.positions[self.indices]

    @position.setter
    def position(self, value):
        if len(self.indices) == 1:
            self.set_position(value)
        else:
            positions = self.parent.positions
            positions[self.indices] = value
            self.parent.positions[self.indices] = positions

    def get_position(self):
        """
        Get global position.
        """
        from batoms.utils import local2global

        position = local2global(
            np.array([self.local_position]), np.array(self.obj.matrix_world)
        )
        return position[0]

    def set_position(self, value):
        """
        Set global position to local vertices
        """
        object_mode()
        from batoms.utils import local2global

        position = np.array([value])
        position = local2global(
            position, np.array(self.obj.matrix_world), reversed=True
        )
        # rashpe to (natoms*3, 1) and use forseach_set
        self.vertice.co = position[0]
        self.obj.data.update()
        update_object(self.obj)

    @property
    def scale(self):
        return self.get_attribute("scale")

    @scale.setter
    def scale(self, value):
        self.set_attribute("scale", value)

    @property
    def show(self):
        return self.get_attribute("show")

    @show.setter
    def show(self, value):
        self.set_attribute("show", value)

    def get_attribute(self, key):
        """Helper function to get attribute

        Args:
            key (str): name of the attribute

        Returns:
            _type_: _description_
        """
        if len(self.indices) == 1:
            value = self.parent.get_attribute(key, index=self.indices[0])
        else:
            value = self.parent.get_attribute_with_indices(key, self.indices)
        return value

    def set_attribute(self, key, value):
        """Helper function to set attribute

        Args:
            key (str): _description_
            value (_type_): _description_
        """
        if len(self.indices) == 1:
            self.parent.set_attribute(key, value, index=self.indices[0])
        else:
            self.parent.set_attribute_with_indices(key, self.indices, value)
