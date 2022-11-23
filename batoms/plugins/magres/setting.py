"""
Visualising tensors using ellipsoids.

1) magnetic shielding tensor
2) chemical shift tensor

"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)



class MagresSettings(Setting):
    def __init__(self, label, parent, scale=1.4,
                 resolution=0.5,
                 ) -> None:
        """Marges object
        The Marges object store the Marges information.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        scale: float
            length of the largest tensor principal component [angstroms], default is 3.
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Bmagres'
        self.scale = scale
        self.parent = parent
        self.ms_name = '%s_ms' % self.label
        self.cs_name = '%s_cs' % self.label
        self.resolution = resolution
        if len(self) == 0:
            self['1'] = {'select': 'all'}

    def add(self, name, datas={}):
        self[name] = datas

    def set_collection(self, label):
        """
        """
        if not bpy.data.collections.get(label):
            coll = bpy.data.collections.new(label)
            self.parent.batoms.coll.children.link(coll)
            coll.batoms.type = 'MAGRES'
            coll.batoms.label = label

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "name  type select  scale   resolution    color  \n"
        for magres in self.bpy_setting:
            s += "{:4s}   {:6s}   {:6s}   {:1.3f}  {:1.3f} ".format(
                magres.name, magres.type, magres.select, magres.scale, magres.resolution)
            s += "[{:1.1f}  {:1.1f}  {:1.1f}  {:1.1f}] \n".format(
                magres.color[0], magres.color[1], magres.color[2], magres.color[3])
        s += "-"*60 + "\n"
        return s
