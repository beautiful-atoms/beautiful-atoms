"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting
from time import time

default_colors = [(1, 1, 0, 0.8), (0.0, 0.0, 1.0, 0.8)]


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
        Setting.__init__(self, label, coll_name='%s_isosurface' % label)
        self.label = label
        self.batoms = batoms
        self.parent = parent
        self.name = 'bisosurface'
        # add a default level
        if isosurfacesetting is not None:
            for key, data in isosurfacesetting.items():
                self[key] = data
        volume = self.batoms.volume
        if volume is not None:
            if len(self) == 0:
                self['1'] = {'level': volume.max()/8, 'color': [1, 1, 0, 0.8]}

    def set_default(self):
        """
        """
        for sp, data in self.species.items():
            self[sp] = {'color': np.append(
                data['color'][:3], 0.3), 'level': 0.005}

    def add(self, isosurfacepair):
        for key in isosurfacepair:
            self.set_default(key)

    def remove_isosurfaces(self, isosurfacepair):
        for key in isosurfacepair:
            name = '%s-%s' % (key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center     level                color  \n"
        for iso in self.collection:
            s += "{0:10s}   {1:1.6f}  ".format(iso.name, iso.level)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}] \n".format(
                iso.color[0], iso.color[1], iso.color[2], iso.color[3])
        s += "-"*60 + "\n"
        return s
