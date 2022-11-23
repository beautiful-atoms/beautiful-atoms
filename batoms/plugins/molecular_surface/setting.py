"""
Molecular surface:
1) van der Waals surface
2) Accessible surface area (SAS)
3) Solvent-excluded surface (SES) or Connolly surface
"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)



class MolecularSurfaceSettings(Setting):
    def __init__(self, label, parent, probe=1.4,
                 resolution=0.5,
                 ) -> None:
        """MS object
        The MS object store the MS information.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        probe: float
            probe radius, in [0,inf], default is 1.4
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Bmolecularsurface'
        self.probe = probe
        self.parent = parent
        self.sas_name = '%s_sas' % self.label
        self.ses_name = '%s_ses' % self.label
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
            coll.batoms.type = 'MS'
            coll.batoms.label = label

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "name  type select  probe   resolution    color  \n"
        for ms in self.bpy_setting:
            s += "{:4s}   {:6s}   {:6s}   {:1.3f}  {:1.3f} ".format(
                ms.name, ms.type, ms.select, ms.probe, ms.resolution)
            s += "[{:1.1f}  {:1.1f}  {:1.1f}  {:1.1f}] \n".format(
                ms.color[0], ms.color[1], ms.color[2], ms.color[3])
        s += "-"*60 + "\n"
        return s
