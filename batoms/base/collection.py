import bpy
import numpy as np
from batoms.utils.butils import object_mode, set_look_at, get_nodes_by_name
# from time import time

class BaseCollection():

    def __init__(self, coll_name):
        self.coll_name = coll_name

    def set_collection(self, label):
        """
        build main collection and its child collections.
        """
        if bpy.data.collections.get(label):
            raise Exception("Failed, the name %s already in use!" % label)
        coll = bpy.data.collections.new(label)
        return coll

    @property
    def coll(self):
        return self.get_coll()

    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            raise('No collection: %s' % self.coll_name)
        return coll

    @property
    def scene(self):
        return self.get_scene()

    def get_scene(self):
        return bpy.data.scenes['Scene']

    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H2o molecule by a vector [0, 0, 5]

        >>> h2o.translate([0, 0, 5])
        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)

    def rotate(self, angle, axis='Z', orient_type='GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:

        >>> h2o.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        angle = angle/180.0*np.pi
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(),
                                 orient_type=orient_type)

    def mirror(self, axis='Z', orient_type='GLOBAL'):
        """mirror atomic based on a axis.

        Parameters:

        angle: float
            Angle that the atoms is mirrord around the axis.
        axis: str
            Constraint Axis: 'X', 'Y' or 'Z'.

        For example, mirror h2o using 'YZ' plane:

        >>> h2o.mirror('X')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        constraint_axis = [False, False, False]
        constraint_axis['XYZ'.index(axis.upper())] = True
        bpy.ops.transform.mirror(constraint_axis=constraint_axis,
                                 orient_type=orient_type)


def tuple2string(index):
    if isinstance(index, (str, int, float)):
        return str(index)
    name = str(index[0])
    for key in index[1:]:
        name += '-%s' % key
    return name


class Setting():
    """
    Setting object

    The Setting object store the additional information for batoms.
    The inforation is a colleciton of data. For example, species, bond pairs.
    In Batoms, collection properties are added to a bpy.type.collection 
    or a bpy.type.object.

    Setting has "add", "remove", "find" and "extend" operators.

    
    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.
        
    """

    def __init__(self, label, coll_name=None, obj_name=None) -> None:
        self.label = label
        self.name = 'base'
        self.coll_name = coll_name
        self.obj_name = obj_name

    def get_data(self):
        """Return a dict of all data in the collection.

        Returns:
            dict: dict of properties.
        """
        data = {}
        for b in self.collection:
            data[b.name] = b
        return data

    @property
    def coll(self):
        """collection

        Returns:
            bpy.type.collection: collection which setting atached to.
        """
        return self.get_coll()

    def get_coll(self):
        coll = bpy.data.collections.get(self.coll_name)
        if coll is None:
            coll = bpy.data.collections.new(self.coll_name)
            bpy.data.scenes['Scene'].collection.children.link(coll)
        return coll

    @property
    def obj(self):
        """object

        Returns:
            bpy.type.object: object which setting atached to.
        """
        return self.get_obj()

    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.' % self.obj_name)
        return obj

    @property
    def collection(self):
        """Collection properties

        Returns:
            bpy.props.CollectionProperty: colleciton of the properties
        """
        return self.get_collection()

    def get_collection(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            collection = getattr(coll.batoms, self.name)
        elif self.obj_name:
            obj = bpy.data.objects.get(self.obj_name)
            collection = getattr(obj.batoms, self.name)
        else:
            raise KeyError("The collection property {}not exist!".format(self.name))
        return collection

    @property
    def data(self):
        return self.get_data()

    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting' % (name, self.name))
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
                raise Exception('%s is not in %s' % (name, self.name))

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
            s += '{0:10s}  \n'.format(
                b.name)
        s += '-'*60 + '\n'
        return s

    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink=True)
