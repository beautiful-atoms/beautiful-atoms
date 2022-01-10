from ase import data
import bpy
import numpy as np
from batoms.butils import object_mode, set_look_at

class BaseObject():
    
    def __init__(self, obj_name, bobj_name):
        self.obj_name = obj_name
        self.bobj_name = bobj_name
    
    @property
    def obj(self):
        return self.get_obj()
    
    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.'%self.obj_name)
        return obj
    
    @property
    def bobj(self):
        return self.get_bobj()
    
    def get_bobj(self):
        bobj = getattr(self.obj.batoms, self.bobj_name)
        return bobj
    
    @property
    def location(self):
        return self.get_location()
    
    def get_location(self):
        return np.array(self.obj.location)
    
    @location.setter    
    def location(self, location):
        self.set_location(location)
    
    def set_location(self, location):
        self.obj.location = location
    
    @property
    def type(self):
        return self.get_type()
    # this is weak, because not work for mesh object
    
    @type.setter
    def type(self, type):
        self.set_type(type)
    
    def get_type(self):
        return self.obj.data.type
    
    def set_type(self, type):
        self.obj.data.type = type.upper()
    
    @property
    def hide(self):
        return self.get_hide()
    
    @hide.setter    
    def hide(self, state):
        self.set_hide(state)
    
    def get_hide(self):
        return self.obj.hide_get()
    
    def set_hide(self, state):
        self.obj.hide_render = state
        self.obj.hide_set(state)
    
    @property
    def select(self):
        return self.get_select()
    
    @select.setter    
    def select(self, state):
        self.set_select(state)
    
    def get_select(self):
        return self.obj.select_get()
    
    def set_select(self, state):
        self.obj.select_set(state)
    
    @property
    def scene(self):
        return self.get_scene()
    
    def get_scene(self):
        return bpy.data.scenes['Scene']
    
    @property
    def look_at(self):
        return self.get_look_at()
    
    @look_at.setter    
    def look_at(self, look_at):
        self.set_look_at(look_at)
    
    def get_look_at(self):
        return np.array(self.bobj.look_at)
    
    def set_look_at(self, look_at):
        self.bobj.look_at = look_at
        set_look_at(self.obj, look_at, roll = 0.0)
    
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h.translate([0, 0, 5])
        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), 
                        orient_type = orient_type)    

class BaseCollection():
    
    def __init__(self, coll_name):
        self.coll_name = coll_name
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            coll = bpy.data.collections.new(self.coll_name)
        return coll
    
    @property
    def scene(self):
        return self.get_scene()
    
    def get_scene(self):
        return bpy.data.scenes['Scene']
    
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h.translate([0, 0, 5])
        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), 
                        orient_type = orient_type)    


def tuple2string(index):
    if isinstance(index, (str, int, float)):
        return str(index)
    name = str(index[0])
    for key in index[1:]:
        name += '-%s'%key
    return name

class Setting():
    """
    Setting object

    The Setting object store the other information for batoms.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.
    """
    
    def __init__(self, label, coll_name = None) -> None:
        self.label = label
        self.name = 'base'
        if coll_name is None:
            self.coll_name = self.label
        else:
            self.coll_name = coll_name
    
    def get_data(self):
        data = {}
        for b in self.collection:
            data[b.name] = b
        return data
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            coll = bpy.data.collections.new(self.coll_name)
            bpy.data.scenes['Scene'].collection.children.link(coll)
        return coll

    @property
    def collection(self):
        return self.get_collection()
    
    def get_collection(self):
        collection = getattr(self.coll.batoms, self.name)
        return collection
    
    @property
    def data(self):
        return self.get_data()
    
    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting'%(name, self.name))
        return item
    
    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
    
    def from_dict(self, datas):
        if isinstance(datas, dict):
            datas = [datas]
        for data in datas:
            subset = self.find(data['name'])
            if subset is None:
                subset = self.collection.add()
            for key, value in data.items():
                setattr(subset, key, value)
            subset.label = self.label
            subset.flag = True
    
    def copy(self, label):
        object_mode()
        setting = self.__class__(label)
        for key, b in self.data.items():
            setting[key] = b.as_dict()
        return setting
    
    def remove(self, index):
        """
        index: list of tuple
        """
        if isinstance(index, (str, int, float)):
            index = [(index)]
        if isinstance(index[0], (str, int, float)):
            index = [index]
        for key in index:
            name = tuple2string(key)
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
            else:
                raise Exception('%s is not in %s'%(name, self.name))
    
    def __delitem__(self, index):
        self.remove(index)
    
    def __add__(self, other):
        self += other
        return self
    
    def __iadd__(self, other):
        self.extend(other)
        return self
    
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_dict()
        return self
    
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    
    def __len__(self):
        return len(self.collection)
    
    def find(self, name):
        i = self.collection.find(str(name))
        if i == -1:
            # print('%s is not in %s'%(name, self.name))
            return None
        else:
            return self.collection[i]
    
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name\n'
        for b in self.collection:
            s += '{0:10s}  \n'.format(\
                b.name)
        s += '-'*60 + '\n'
        return s
