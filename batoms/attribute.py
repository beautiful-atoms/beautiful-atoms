"""

"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)



class Attributes(Setting):
    """Attribute information for bpy.data.object.

    It should be bounded to an object, because each object
    has different attributes. One collection, however, could have
    many objects.

    Args:
        Setting (_type_): _description_
    """
    def __init__(self, label, parent, obj_name = None,
                 ) -> None:
        """Attributes object
        The Attributes object store the attributes information.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        """
        Setting.__init__(self, label, obj_name=obj_name)
        self.label = label
        self.name = 'batoms'
        self.parent = parent

    def get_bpy_setting(self):
        if self.obj_name:
            coll = bpy.data.objects.get(self.obj_name)
            data = getattr(coll, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.settings_attribute

    def add(self, datas=None):
        """Add new attribute

        Args:
            datas (list, dict): attributes
        """
        if datas is None:
            return
        if isinstance(datas, dict):
            datas = [datas]
        for data in datas:
            # print(data)
            # print("Add new attribute: {}".format(data))
            # add attribute into _attributes collection
            self[data["name"]] = data

            # create attribute in mesh
            att = self[data["name"]]
            mesh = self.parent.obj.data
            # is single value, no delimiter needed
            if data['dimension'] == 0:
                mesh.attributes.new(name=data["name"],
                        type=att.data_type, domain=att.domain)
                logger.debug("add attribute: {} {} {}"
                    .format(att.name, att.data_type, att.domain))
            else:
                # is array, scatter to new attribute,
                # with name = "{}{}{}".format(name, delimiter, index)
                # default delimiter = "@". Compare new name with other attribute's name
                # add "@" recursively if needed.
                M = np.product(data["shape"][0:data["dimension"]])
                name_list = list(self.keys())
                delimiter = self.find_delimiter(name_list, data["name"], M, att.delimiter)
                # update delimiter
                att.delimiter = delimiter
                for i in range(M):
                    name = "{}{}{}".format(att.name, att.delimiter, i)
                    mesh.attributes.new(name=name,
                                type=att.data_type, domain=att.domain)

    def find_delimiter(self, name_list, name, M, delimiter):
        """find delimiter recursively

        Args:
            name_list (list):
                list of the name of attributes already in the collection
            name (string): the name used for scatter
            M (int): number of new attribute
            delimiter (string): default @

        Returns:
            string: "@"*n (n>=1)
        """
        # print("find delimiter: ", name, name_list)
        for i in range(M):
            new_name = "{}{}{}".format(name, delimiter, i)
            if new_name in name_list:
                delimiter += "@"
                delimiter = self.find_delimiter(name_list, name, M, delimiter)
                return delimiter
                # print("new delimiter: ", delimiter)
        return delimiter


    def from_array(self, name, data):
        """_summary_

        Args:
            name (_type_): _description_
            data (_type_): _description_

        Returns:
            _type_: _description_
        """
        from batoms.utils import type_py_to_blender
        # print("from_array: ", name, data)
        array = np.asarray(data)
        type_py = type(array.flat[0])
        # print(name, dtype)
        dtype_bl = type_py_to_blender(type_py)
        if dtype_bl is False:
            logger.critical('Attribute: {}, {} is not supported.'.format(name, type_py))
            return False
        shape = array[0].shape
        dimension = len(shape)
        data = {"name": name, "data_type": dtype_bl, "domain": "POINT",
                "dimension": dimension, "shape": shape}
        self.add(data)
        return True


    def set_collection(self, label):
        """
        """
        if not bpy.data.collections.get(label):
            coll = bpy.data.collections.new(label)
            self.parent.batoms.coll.children.link(coll)
            coll.batoms.type = 'ATTRIBUTE'
            coll.batoms.label = label

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "{:20s}{:10s}{:10s}{:10s}   {:20s}\n".format("Name", "Type", "Domain", "Dimension", "Shape")
        for att in self.bpy_setting:
            s += "{:20s}{:10s}{:10s}{:10d}  [".format(
                att.name, att.data_type, att.domain, att.dimension)
            for i in range(att.dimension):
                s += "  {}  ".format(att.shape[i])
            s += "] \n"
        s += "-"*60 + "\n"
        return s
