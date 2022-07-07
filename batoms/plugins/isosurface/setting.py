"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting
from time import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class IsosurfaceSettings(Setting):
    def __init__(self, label, batoms=None, parent=None,
                 isosurfacesetting=None) -> None:
        """IsosurfaceSetting object

            The IsosurfaceSetting object store the isosurface information.

            Parameters:

            label: str
                The label define the batoms object that a Setting belong to.


        Args:
            label (_type_):
                _description_
            batoms (_type_, optional):
                _description_. Defaults to None.
            isosurfacesetting (_type_, optional):
                _description_. Defaults to None.
        """
        Setting.__init__(self, label, coll_name='%s' % label)
        self.label = label
        self.batoms = batoms
        self.parent = parent
        self.name = 'Bisosurface'
        # add a default level
        if isosurfacesetting is not None:
            for key, data in isosurfacesetting.items():
                self[key] = data

    def add(self, isosurface):
        if isinstance(isosurface, str):
            self[isosurface] = {}
        elif isinstance(isosurface, dict):
            self[isosurface['name']] = isosurface

    def remove_isosurfaces(self, isosurface):
        for key in isosurface:
            name = '%s-%s' % (key[0], key[1])
            i = self.bpy_setting.find(name)
            if i != -1:
                self.bpy_setting.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center     level                color  \n"
        for iso in self.bpy_setting:
            s += "{0:10s}   {1:1.6f}  ".format(iso.name, iso.level)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}] \n".format(
                iso.color[0], iso.color[1], iso.color[2], iso.color[3])
        s += "-"*60 + "\n"
        return s
