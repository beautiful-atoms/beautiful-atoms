"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting


class CavitySettings(Setting):
    def __init__(self, label, batoms=None, parent=None,
                 cavitysetting=None) -> None:
        """CavitySetting object

            The CavitySetting object store the cavity information.

            Parameters:

            label: str
                The label define the batoms object that a Setting belong to.


        Args:
            label (_type_):
                _description_
            batoms (_type_, optional):
                _description_. Defaults to None.
            cavitysetting (_type_, optional):
                _description_. Defaults to None.
        """
        Setting.__init__(self, label, coll_name='%s_cavity' % label)
        self.label = label
        self.batoms = batoms
        self.parent = parent
        self.name = 'bcavity'
        # add a default level
        if cavitysetting is not None:
            for key, data in cavitysetting.items():
                self[key] = data
        volume = self.batoms.volume
        if volume is not None:
            if len(self) == 0:
                self['1'] = {'level': volume.max()/8, 'color': [1, 1, 0, 0.8]}

    def add(self, cavity):
        if isinstance(cavity, str):
            self[cavity] = {}
        elif isinstance(cavity, dict):
            self[cavity['name']] = cavity

    def remove_cavitys(self, cavity):
        for key in cavity:
            name = '%s-%s' % (key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center     level                color  \n"
        for cav in self.collection:
            s += "{:10s}   {:1.6f}   {:1.6f}".format(cav.name, cav.min, cav.max)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}] \n".format(
                cav.color[0], cav.color[1], cav.color[2], cav.color[3])
        s += "-"*60 + "\n"
        return s
