"""

"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting


class Attributes(Setting):
    def __init__(self, label, parent
                 ) -> None:
        """Attributes object
        The Attributes object store the attributes information.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        """
        Setting.__init__(self, label)
        self.label = label
        self.name = 'battributes'
        self.parent = parent

    def add(self, name, datas={}):
        self[name] = datas

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
        s = "name  type domain dimension    shape\n"
        for att in self.collection:
            s += "{:4s}   {:6s}   {:6s}   {:}   [".format(
                att.name, att.type, att.domain, att.dimension)
            for i in range(att.dimension):
                s += "  {}  ".format(att.shape[i])
            s += "]"
        s += "-"*60 + "\n"
        return s
    

