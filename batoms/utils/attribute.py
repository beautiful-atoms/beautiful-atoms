"""
Get and set attribute
"""
import numpy as np
import logging

# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def get_mesh_attribute_bmesh(obj, key, index=None):
    """Get the attribute of mesh by name usinb bmesh method.
    When use this function:
    1) Edit mode, use bmesh.
    2) For STRING attribute

    Args:
        key (string): name of the attribute
        index (bool, int): index of the data, used to get singe attribute value

    Raises:
        KeyError: _description_

    Returns:
        array: _description_
    """
    import bmesh
    from batoms.utils.butils import get_att_length, get_bmesh_layer, get_bmesh_domain
    from batoms.utils import type_blender_to_py

    # get the mesh
    me = obj.data
    # get attribute data type for the attribute
    att = me.attributes.get(key)
    if att is None:
        raise KeyError("{} is not exist.".format(key))
    # get type and domain
    dtype = att.data_type
    # get attribute length based on domain
    n = get_att_length(obj.data, att)
    # check mode
    if obj.mode == "EDIT":
        bm = bmesh.from_edit_mesh(obj.data)
        domain = get_bmesh_domain(bm, att)
        domain.ensure_lookup_table()
        layer = get_bmesh_layer(domain, key, dtype)
        if index is not None:
            value = domain[index][layer]
            if dtype == "STRING":
                value = value.decode()
            return np.array([value])
        # init attribute array
        attribute = np.zeros(n, dtype=type_blender_to_py(dtype, str="U20"))
        for i in range(n):
            # print(domain[i][layer])
            attribute[i] = domain[i][layer]
    else:
        bm = bmesh.new()
        bm.from_mesh(me)
        domain.ensure_lookup_table()
        domain = get_bmesh_domain(bm, att)
        if index is not None:
            value = domain[index][layer]
            if dtype == "STRING":
                value = value.decode()
            return np.array([value])
        n = len(domain)
        for i in range(n):
            attribute[i] = domain[i][layer]
        bm.free()
    return attribute


def get_mesh_attribute(obj, key, index=None):
    """Get the attribute of mesh by name usinb bmesh method.

    When use this function:
    1) Object mode, read attribute directly
    2) For BOOLEAN attribute, because bmesh does not support BOOLEAN.

    Args:
        key (string): name of the attribute
        index (bool, int): index of the data, used to get singe attribute value

    Raises:
        KeyError: _description_

    Returns:
        array: _description_
    """
    from batoms.utils.butils import get_att_length
    from batoms.utils import type_blender_to_py

    # get the mesh
    me = obj.data
    # get attribute data type for the attribute
    att = me.attributes.get(key)
    # print(f"Attribute: {att.data_type} {att.domain}")
    if att is None:
        raise KeyError("{} is not exist.".format(key))
    dtype = att.data_type
    # get single attribute value
    if index is not None:
        value = att.data[index].value
        if dtype == "STRING":
            value = value.decode()
        return np.array([value])
    else:
        # get attribute length based on domain
        n = get_att_length(obj.data, att)
        # init
        attribute = np.zeros(n, dtype=type_blender_to_py(dtype, str="U20"))
        if dtype == "STRING":
            for i in range(n):
                attribute[i] = att.data[i].value
        elif dtype in ["INT", "FLOAT", "BOOLEAN"]:
            att.data.foreach_get("value", attribute)
        elif dtype == "FLOAT2":
            attribute = np.zeros(n * 2, dtype="float")
            att.data.foreach_get("vector", attribute)
            attribute = attribute.reshape(n, 2)
        elif dtype == "FLOAT_VECTOR":
            attribute = np.zeros(n * 3, dtype="float")
            att.data.foreach_get("vector", attribute)
            attribute = attribute.reshape(n, 3)
        elif dtype in ["FLOAT_COLOR"]:
            attribute = np.zeros(n * 4, dtype="float")
            att.data.foreach_get("color", attribute)
            attribute = attribute.reshape(n, 4)
        elif dtype in ["QUATERNION"]:
            attribute = np.zeros(n * 4, dtype="float")
            att.data.foreach_get("value", attribute)
            attribute = attribute.reshape(n, 4)
        else:
            raise KeyError("Attribute type: %s is not support." % dtype)
        attribute = np.array(attribute)
    return attribute


def set_mesh_attribute_bmesh(obj, key, value, index=None):
    """Set mesh attribute using bmesh method.

    Args:
        obj (bpy.type.object): obj
        key (str): name of the attribute
        value (np.array): value of the attribute
        index (bool, int): index of the data, used to set singe attribute value
    """
    import bmesh
    from batoms.utils.butils import get_bmesh_domain, get_bmesh_layer

    me = obj.data
    # get attribute type
    att = me.attributes.get(key)
    dtype = att.data_type
    if dtype == "STRING":
        value = np.char.encode(value)
    # check object mode
    if obj.mode == "EDIT":
        bm = bmesh.from_edit_mesh(me)
        domain = get_bmesh_domain(bm, att)
        domain.ensure_lookup_table()
        layer = get_bmesh_layer(domain, key, dtype)
        if index is not None:
            domain[index][layer] = value
        else:
            n = len(domain)
            for i in range(n):
                domain[i][layer] = value[i]
        bmesh.update_edit_mesh(me)
    else:
        bm = bmesh.new()
        bm.from_mesh(me)
        domain = get_bmesh_domain(bm, att)
        domain.ensure_lookup_table()
        layer = get_bmesh_layer(domain, key, dtype)
        if index is not None:
            domain[index][layer] = value
        else:
            n = len(domain)
            for i in range(n):
                domain[i][layer] = value[i]
        bm.to_mesh(me)
        bm.free()


def set_mesh_attribute(obj, key, value, index=None):
    """Set mesh attribute using bmesh method

    Args:
        obj (bpy.type.object): obj
        key (str): name of the attribute
        value (np.array): value of the attribute
        index (bool, int): index of the data, used to set singe attribute value
    """
    from batoms.utils.butils import get_att_length

    me = obj.data
    # get attribute
    att = me.attributes.get(key)
    # print(f"Attribute: {att.data_type} {att.domain}")
    if index is not None:
        att.data[index].value = value[0]
    else:
        # get attribute domain and length
        n = get_att_length(obj.data, att)
        if att.data_type == "STRING":
            for j in range(n):
                att.data[j].value = value[j]
        elif att.data_type == "FLOAT2":
            value = value.reshape((n * 2, 1))
            att.data.foreach_set("vector", value)
        elif att.data_type == "FLOAT_VECTOR":
            value = value.reshape((n * 3, 1))
            att.data.foreach_set("vector", value)
        elif att.data_type in ["FLOAT_COLOR"]:
            value = value.reshape((n * 4, 1))
            att.data.foreach_set("color", value)
        elif att.data_type in ["QUATERNION"]:
            value = value.reshape((n * 4, 1))
            att.data.foreach_set("value", value)
        else:
            att.data.foreach_set("value", value)
