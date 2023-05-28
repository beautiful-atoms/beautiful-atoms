"""
Add Volumetric data.
"""
import bpy
import numpy as np
from time import time
from batoms.base.collection import Setting
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)



class VolumetricData(Setting):
    def __init__(self, label, volume = None, parent = None) -> None:
        """VolumetricData object
        The VolumetricData object store the Volumetric Data information.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'batoms'
        self.parent = parent
        if volume is not None:
            for key, data in volume.items():
                self[key] = data

    def get_ui_list_index(self):
        return self.bpy_data.ui_list_index_volumetric_data

    def set_ui_list_index(self, value):
        self.bpy_data.ui_list_index_volumetric_data = value

    def get_bpy_setting(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            data = getattr(coll, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.settings_volumetric_data


    def __getitem__(self, key):
        item = self.get_volume(key)
        if item is None:
            raise Exception('%s not in %s setting' % (key, self.name))
        return item

    def __setitem__(self, key, volume):
        """
        Set properties
        """
        subset = self.find(key)
        if subset is None:
            subset = self.bpy_setting.add()
            subset.name = key
            self.ui_list_index = len(self) - 1
        subset.label = self.label
        subset.flag = True
        self.build_object(subset, volume)

    def add(self, name, datas={}):
        self[name] = datas

    def set_collection(self, label):
        """
        """
        if not bpy.data.collections.get(label):
            coll = bpy.data.collections.new(label)
            self.parent.batoms.coll.children.link(coll)
            coll.batoms.type = 'VOLUME'
            coll.batoms.label = label

    def build_object(self, setting, volume):
        """Save volumetric data as a mesh

        Args:
            volume (array):
                volumetric data, e.g. electron density
        """
        # remove old volume point
        # tstart = time()
        if volume is None:
            return
        name = "{}_volume_{}".format(self.label, setting.name)
        if name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[name], do_unlink=True)
        if name in bpy.data.meshes:
            bpy.data.meshes.remove(bpy.data.meshes[name], do_unlink=True)
        shape = volume.shape
        volume = volume.reshape(-1, 1)
        npoint = len(volume)
        setting.shape = shape
        setting.npoint = npoint
        dn = 3 - npoint % 3
        verts = np.append(volume, np.zeros((dn, 1)), axis=0)
        verts = verts.reshape(-1, 3)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(verts, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        obj.data = mesh
        obj.parent = self.parent.obj
        obj.batoms.type = 'VOLUME'
        obj.batoms.volume.shape = shape
        self.coll.objects.link(obj)
        obj.hide_set(True)
        obj.hide_render = True
        # print('Draw volume: {0:1.2f}'.format(time() - tstart))

    def get_volume(self, name):
        """Retrieve volume data from a mesh

        Returns:
            array: Volumetric data
        """
        # tstart = time()
        setting = self.find(name)
        if setting is None:
            return None
        obj = bpy.data.objects.get('{}_volume_{}'.format(self.label, setting.name))
        if obj is None:
            return None
        n = len(obj.data.vertices)
        volume = np.empty(n*3, dtype=np.float64)
        obj.data.vertices.foreach_get('co', volume)
        volume = volume.reshape(-1, 1)
        shape = setting.shape
        npoint = np.product(shape)
        volume = volume[:npoint]
        volume = volume.reshape(shape)
        # print('Read volume: {0:1.2f}'.format(time() - tstart))
        return volume

    def __imul__(self, m):
        import numpy as np
        for volume in self.bpy_setting:
            data = self[volume.name]
            self[volume.name] = np.tile(data, m)
        return self

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "name          npoint        shape  \n"
        for volume in self.bpy_setting:
            s += "{:4s}   {:10d} ".format(
                volume.name, volume.npoint)
            s += "[{:5d}  {:5d}  {:5d}] \n".format(
                volume.shape[0], volume.shape[1], volume.shape[2])
        s += "-"*60 + "\n"
        return s
