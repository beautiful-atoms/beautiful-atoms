"""
"""
import bpy
from batoms.base.collection import Setting, tuple2string
import numpy as np
from time import time


class CrystalShapeSettings(Setting):
    """
    PlaneSetting object

    The PlaneSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """

    def __init__(self, label, parent=None, plane=None) -> None:
        Setting.__init__(self, label, coll_name='%s' % label)
        self.name = 'Bcrystalshape'
        self.parent = parent
        if plane is not None:
            for key, data in plane.items():
                self[key] = data

    def get_collection(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            collection = getattr(coll, self.name)
        elif self.obj_name:
            obj = bpy.data.objects.get(self.obj_name)
            collection = getattr(obj, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return collection.settings

    def __setitem__(self, key, setdict):
        """
        Set properties
        """
        subset = self.find(key)
        if subset is None:
            subset = self.bpy_setting.add()
            name = tuple2string(key)
            subset.name = name
            self.ui_list_index = len(self) - 1
        subset.indices = key
        subset.name = name
        subset.flag = True
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label

    def add(self, indices):
        self[indices] = {'indices': indices}

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Indices   distance  symmetry  \n"
        for p in self.bpy_setting:
            s += "{0:10s}   {1:1.3f}   ".format(p.name, p.distance)
            s += "{:10s} \n".format(
                str(p.symmetry))
        s += "-"*60 + "\n"
        return s

    def get_symmetry_indices(self):
        from batoms.utils import get_equivalent_indices
        if self.no == 1:
            return
        for p in self:
            if p.symmetry:
                indices = get_equivalent_indices(self.no, p.indices)
                setdict = p.as_dict()
                for index in indices:
                    name = tuple2string(index)
                    p1 = self.find(name)
                    if p1 is None:
                        p1 = self.bpy_setting.add()
                        for key, value in setdict.items():
                            setattr(p1, key, value)
                        p1.name = name
                        p1.indices = index
                        p1.label = self.label
                        p1.flag = True

    @property
    def no(self, ):
        return self.parent.batoms.get_spacegroup_number()

    @no.setter
    def no(self, no):
        self.no = no
