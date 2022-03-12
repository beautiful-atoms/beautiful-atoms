"""

Lattice Planes

To insert lattice planes in structural models.
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

    def __init__(self, label, parent = None, plane=None) -> None:
        Setting.__init__(self, label, coll_name='%s_plane' % label)
        self.name = 'bcrystalshape'
        self.parent = parent
        if plane is not None:
            for key, data in plane.items():
                self[key] = data

    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        p = self.find(name)
        if p is None:
            p = self.collection.add()
        p.indices = index
        p.name = name
        p.flag = True
        for key, value in setdict.items():
            setattr(p, key, value)
        p.label = self.label

    def add(self, indices):
        self[indices] = {'indices': indices}

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Indices   distance  symmetry  \n"
        for p in self.collection:
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
                        p1 = self.collection.add()
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