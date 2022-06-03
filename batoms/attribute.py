"""

"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting


class Attributes(Setting):
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
        self.name = 'battribute'
        self.parent = parent

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
            self[data["name"]] = data
            att = self[data["name"]]
            mesh = self.parent.obj.data
            if data['dimension'] == 0:
                mesh.attributes.new(name=data["name"],
                        type=att.type, domain=att.domain)
            else:
                M = np.product(data["shape"][0:data["dimension"]])
                for i in range(M):
                    name = "{}@{}".format(data["name"], i)
                    mesh.attributes.new(name=name,
                                type=att.type, domain=att.domain)
    
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
            print('Attribute: {}, {} is not supported.'.format(name, type_py))
            return False
        shape = array[0].shape
        dimension = len(shape)
        data = {"name": name, "type": dtype_bl, "domain": "POINT", 
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
        for att in self.collection:
            s += "{:20s}{:10s}{:10s}{:10d}  [".format(
                att.name, att.type, att.domain, att.dimension)
            for i in range(att.dimension):
                s += "  {}  ".format(att.shape[i])
            s += "] \n"
        s += "-"*60 + "\n"
        return s
