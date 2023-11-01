import bpy
import numpy as np
from batoms.utils.butils import (get_nodes_by_name, object_mode, set_look_at,
                                 update_object)

from time import time
import bmesh
import logging

# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


default_object_attributes = [
]

default_object_datas = {
}


class BaseObject():

    def __init__(self, obj_name, btype="batoms"):
        self.obj_name = obj_name
        self.btype = btype

    @property
    def obj(self):
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.' % self.obj_name)
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
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)

    def rotate(self, angle, axis='Z', orient_type='GLOBAL'):
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
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(),
                                 orient_type=orient_type)

    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink=True)

    def delete_material(self, name):
        if name in bpy.data.materials:
            obj = bpy.data.materials.get(name)
            bpy.data.materials.remove(obj, do_unlink=True)

    def update_mesh(self, obj=None):
        import bmesh
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
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode=mode)

    def add_verts(self, count, obj=None):
        import bmesh
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

    def delete_vertices_bmesh(self, index=[], obj=None,):
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
        bmesh.ops.delete(bm, geom=verts_select, context='VERTS')
        bm.to_mesh(obj.data)
        bm.clear()

    @property
    def shape_keys(self):
        base_name = "Basis_%s"%self.obj_name
        if self.obj.data.shape_keys is None and len(self) > 0:
            sk = self.obj.shape_key_add(name=base_name)
        return self.obj.data.shape_keys

    def __len__(self):
        return len(self.obj.data.vertices)


class ObjectGN(BaseObject):
    """
    Object with Geometry Node

    """

    def __init__(self, label, name=None):
        if name:
            self.name = name
            obj_name = '%s_%s' % (label, name)
        else:
            obj_name = label
        BaseObject.__init__(self, obj_name=obj_name)

    def build_object(self, arrays, attributes={}):
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis=False)

    def load(self):
        flag = True
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            flag = False
        return flag

    @property
    def gnodes(self):
        return self.get_gnodes()

    def get_gnodes(self):
        name = 'GeometryNodes_%s' % self.obj_name
        modifier = self.obj.modifiers.get(name)
        if modifier is None:
            self.build_geometry_node()
        return modifier

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        from batoms.utils.butils import build_modifier
        name = 'GeometryNodes_%s' % self.obj_name
        modifier = build_modifier(self.obj, name)

    def vectorDotMatrix(self, gn, vector_output, matrix, name):
        """
        """
        CombineXYZ = get_nodes_by_name(gn.node_group.nodes,
                                       '%s_CombineXYZ_%s' % (self.label, name),
                                       'ShaderNodeCombineXYZ')
        #
        VectorDot = []
        for i in range(3):
            tmp = get_nodes_by_name(gn.node_group.nodes,
                                    '%s_VectorDot%s_%s' % (
                                        self.label, i, name),
                                    'ShaderNodeVectorMath')
            tmp.operation = 'DOT_PRODUCT'
            VectorDot.append(tmp)
            tmp.inputs[1].default_value = matrix[:, i]
            gn.node_group.links.new(vector_output, tmp.inputs[0])
            gn.node_group.links.new(tmp.outputs['Value'], CombineXYZ.inputs[i])
        return CombineXYZ

    def add_geometry_node(self, sp):
        """
        add geometry node for each bond pair
        """
        pass

    @property
    def realize_instances(self):
        return self.get_realize_instances()

    @realize_instances.setter
    def realize_instances(self, state):
        self.set_realize_instances(state)

    def get_realize_instances(self):
        return list(self.coll.batoms.realize_instances)

    def set_realize_instances(self, realize_instances):
        """ Make instancing object real
        # TODO: add make real to geometry node
        """
        #
        nodes = self.gnodes.node_group.nodes
        RealizeInstances = get_nodes_by_name(self.gnodes.node_group.nodes,
                                             '%s_RealizeInstances' % self.label,
                                             'GeometryNodeRealizeInstances')
        if not realize_instances:
            # switch off
            if len(RealizeInstances.outputs[0].links) > 0:
                link = RealizeInstances.outputs[0].links[0]
                self.gnodes.node_group.links.remove(link)
            self.gnodes.node_group.links.new(
                nodes['%s_JoinGeometry' % self.label].outputs[0],
                nodes[1].inputs[0])
        else:
            self.gnodes.node_group.links.new(
                nodes['%s_JoinGeometry' % self.label].outputs[0],
                RealizeInstances.inputs[0])
            self.gnodes.node_group.links.new(
                RealizeInstances.outputs[0],
                nodes[1].inputs[0])
        self.gnodes.node_group.update_tag()

    def update(self, ):
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
        arrays.update({'positions': self.positions[0],
                       })
        return arrays

    @property
    def attributes(self):
        return self.get_attributes()

    @attributes.setter
    def attributes(self, attributes):
        self.set_attributes(attributes)

    def get_attributes(self):
        """Get all attributes of the Batoms object.


        Returns:
            dict: attributes dict
        """
        attributes = {}
        for att in self._attributes:
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
        from batoms.utils import type_blender_to_py
        from batoms.utils.butils import get_att_length
        from batoms.utils.attribute import get_mesh_attribute,get_mesh_attribute_bmesh
        # get the mesh
        obj = self.obj
        att = self._attributes[key]
        if att.dimension == 0:
            if obj.mode == 'EDIT' and att.data_type in ['STRING', 'INT', 'FLOAT']:
                attribute = get_mesh_attribute_bmesh(obj, att.name, index)
            else:
                attribute = get_mesh_attribute(obj, att.name, index)
        else:
            natt = att.natt
            name = "{}{}{}".format(att.name, att.delimiter, 0)
            # init a large array has the size of n*np.product(shape)
            if index is not None:
                n = 1
            else:
                mesh_att = obj.data.attributes.get(name)
                n = get_att_length(obj.data, mesh_att)
            attribute = np.zeros(natt*n, dtype=type_blender_to_py(att.data_type))
            for i in range(natt):
                name = "{}{}{}".format(att.name, att.delimiter, i)
                if obj.mode == 'EDIT' and att.data_type in ['STRING', 'INT', 'FLOAT']:
                    attribute[i*n:(i+1)*n] = get_mesh_attribute_bmesh(obj, name, index)
                else:
                    attribute[i*n:(i+1)*n] = get_mesh_attribute(obj, name, index)
            # reshape to (n, shape)
            attribute = attribute.reshape((n, ) + att.shape, order='F')
        return attribute

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
            if not self._attributes.find(key):
                # try to create a new attribute
                flag = self._attributes.from_array(key, array)
                # if failed, do not add this attribute
                if not flag:
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
        from batoms.utils.butils import get_att_length
        from batoms.utils.attribute import set_mesh_attribute, set_mesh_attribute_bmesh
        tstart = time()
        obj = self.obj
        me = obj.data
        att_coll = self._attributes[key]
        shape = att_coll.shape
        dimension = att_coll.dimension
        delimiter = att_coll.delimiter
        #
        array = np.array(array)
        if len(array.shape) == 0:
            array = np.array([array])
        if len(array) == 0:
            return
        # single value data
        if dimension == 0:
            att = me.attributes.get(key)
            if att.data_type == 'STRING' or (obj.mode == 'EDIT' and att.data_type in ['INT', 'FLOAT']):
                set_mesh_attribute_bmesh(obj, key, array, index)
            else:
                set_mesh_attribute(obj, key, array, index)
        # array data
        else:
            # M is the number of sub-array, for 2x2 array, M is 4
            M = att_coll.natt
            array = array.reshape(-1, 1, order='F')
            for i in range(M):
                sub_key = "{}{}{}".format(key, delimiter, i)
                att = me.attributes.get(sub_key)
                n = get_att_length(obj.data, att)
                if obj.mode == 'EDIT' and att.data_type in ['STRING', 'INT', 'FLOAT']:
                    set_mesh_attribute_bmesh(obj, sub_key, array[i*n:(i+1)*n], index)
                else:
                    set_mesh_attribute(obj, sub_key, array[i*n:(i+1)*n], index)
        logger.debug('Time: {:10s} {:5.2f}'.format(key, time() - tstart))

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
    def local_positions(self):
        return self.get_local_positions()

    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = np.empty(n*3, dtype=np.float64)
        self.obj.data.vertices.foreach_get('co', local_positions)
        local_positions = local_positions.reshape((n, 3))
        return local_positions

    @local_positions.setter
    def local_positions(self, local_positions):
        self.set_local_positions(local_positions)

    def set_local_positions(self, local_positions):
        """
        Set local_positions to local vertices
        """
        object_mode()
        from batoms.utils import local2global
        natom = len(self)
        n = len(local_positions)
        if n != natom:
            raise ValueError('local_positions has wrong shape %s != %s.' %
                             (n, natom))
        # no shape key for empty mesh
        if natom == 0:
            return
        local_positions = local_positions.reshape((natom*3, 1))
        self.shape_keys.key_blocks[0].data.foreach_set(
            'co', local_positions)
        self.obj.data.update()
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.object.mode_set(mode='OBJECT')

    @property
    def positions(self):
        return self.get_positions()

    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)

    def get_positions(self):
        """
        Get global positions.
        """
        from batoms.utils import local2global
        positions = local2global(self.local_positions,
                                 np.array(self.obj.matrix_world))
        return positions

    def set_positions(self, positions):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.utils import local2global
        natom = len(self)
        nposition = len(positions)
        if nposition != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                             (nposition, natom))
        # no shape key for empty mesh
        if natom == 0:
            return
        positions = local2global(positions,
                                 np.array(self.obj.matrix_world),
                                 reversed=True)
        # rashpe to (natoms*3, 1) and use forseach_set
        positions = positions.reshape((natom*3, 1))
        # I don't know why 'Basis' shape keys is not updated when editing mesh,
        # so we edit the 'Basis' shape keys directly.
        # self.obj.data.vertices.foreach_set('co', positions)
        self.shape_keys.key_blocks[0].data.foreach_set(
            'co', positions)
        self.obj.data.update()
        self.update_mesh()

    @property
    def nframe(self):
        return self.get_nframe()

    def get_nframe(self):
        if self.obj.data.shape_keys is None:
            return 0
        nframe = len(self.obj.data.shape_keys.key_blocks)
        return nframe

    @property
    def frames(self):
        return self.get_frames()

    @frames.setter
    def frames(self, frames):
        self.set_frames(frames)

    def get_obj_frames(self, obj, local=True):
        """
        read shape key
        """
        from batoms.utils import local2global
        n = len(self)
        nframe = self.nframe
        frames = np.empty((nframe, n, 3), dtype=np.float64)
        for i in range(nframe):
            positions = np.empty(n*3, dtype=np.float64)
            sk = obj.data.shape_keys.key_blocks[i]
            sk.data.foreach_get('co', positions)
            local_positions = positions.reshape((n, 3))
            if local:
                frames[i] = local_positions
            else:
                global_positions = local2global(local_positions,
                                                np.array(self.obj.matrix_world))
                frames[i] = global_positions
        return frames

    def set_frames_positions(self, frames=None, frame_start=0,
                             only_basis=False):
        name = ''
        obj = self.obj
        self.set_frames(name, obj, frames, frame_start, only_basis)

    def set_obj_frames(self, name, obj, frames=None,
                       frame_start=0, only_basis=False):
        """

        frames: list
            list of positions

        use shape_keys (faster)
        """
        from batoms.utils.butils import add_keyframe_to_shape_key
        if frames is None:
            frames = self._frames
        centers = frames
        nframe = len(centers)
        if nframe == 0:
            return
        # shape key should not be add for empty mesh
        if len(frames[0]) == 0:
            return
        # name = '%s_bond%s'%(self.label, sp)
        # obj = bpy.data.objects.get(name)
        base_name = 'Basis_%s' % (name)
        if obj.data.shape_keys is None:
            sk = obj.shape_key_add(name=base_name)
        elif base_name not in obj.data.shape_keys.key_blocks:
            sk = obj.shape_key_add(name=base_name)
        else:
            sk = obj.data.shape_keys.key_blocks.get(base_name)
        # set basis key
        nvert = len(obj.data.shape_keys.key_blocks[0].data)
        positions = frames[0]
        vertices = positions.reshape(nvert*3)
        sk.data.foreach_set('co', vertices)
        # self.obj.data.update()
        if only_basis:
            self.update_mesh(obj)
            return
        for i in range(1, nframe):
            name = str(i)
            if name not in obj.data.shape_keys.key_blocks:
                sk = obj.shape_key_add(name=name)
                # Add Keyframes, the last one is different
                if i != nframe - 1:
                    add_keyframe_to_shape_key(sk, 'value',
                                              [0, 1, 0],
                                              [frame_start + i - 1,
                                                  frame_start + i,
                                                  frame_start + i + 1])
                else:
                    add_keyframe_to_shape_key(sk, 'value',
                                              [0, 1],
                                              [frame_start + i - 1,
                                                  frame_start + i])
            else:
                sk = obj.data.shape_keys.key_blocks.get(name)
            # Use the local position here
            positions = frames[i]
            vertices = positions.reshape((nvert*3, 1))
            sk.data.foreach_set('co', vertices)
        self.update_mesh(obj)

    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            # if not obj.batoms.flag:
            # raise Exception(
            # "Failed, the name {} already in use and is not
            # Batom object!".format(name))
            bpy.data.objects.remove(obj, do_unlink=True)


class childObjectGN():
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
            raise KeyError('%s object is not exist.' % self.obj_name)
        return obj

    @property
    def vertice(self):
        return self.obj.data.shape_keys.key_blocks[0].data[self.indices[0]]

    @property
    def bm(self):
        if self.obj.mode == 'EDIT':
            bm =bmesh.from_edit_mesh(self.obj.data)
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
        position = local2global(np.array([self.local_position]),
                                np.array(self.obj.matrix_world))
        return position[0]

    def set_position(self, value):
        """
        Set global position to local vertices
        """
        object_mode()
        from batoms.utils import local2global
        position = np.array([value])
        position = local2global(position,
                                np.array(self.obj.matrix_world),
                                reversed=True)
        # rashpe to (natoms*3, 1) and use forseach_set
        self.vertice.co = position[0]
        self.obj.data.update()
        update_object(self.obj)

    @property
    def scale(self):
        return self.get_attribute('scale')

    @scale.setter
    def scale(self, value):
        self.set_attribute('scale', value)

    @property
    def show(self):
        return self.get_attribute('show')

    @show.setter
    def show(self, value):
        self.set_attribute('show', value)

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
