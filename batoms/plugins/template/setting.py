"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting
import logging
logger = logging.getLogger(__name__)



class TemplateSettings(Setting):
    def __init__(self, label, parent) -> None:
        """Template object
        This object store the plugin setting.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Btemplate'
        self.parent = parent
        # initialize a default setting
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
            coll.batoms.type = 'PTemplate'
            coll.batoms.label = label

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "name  type select\n"
        for setting in self.bpy_setting:
            s += "{:4s}   {:6s}   {:6s}".format(
                setting.name, setting.type, setting.select)
        s += "-"*60 + "\n"
        return s
    
    