import bpy
import numpy as np
from time import time
from batoms.butils import object_mode, set_look_at

default_object_attributes = [
        ]

default_object_datas = {
        }

class ObjectGN():
    """
    Object with Geometry Node

    """
    def __init__(self, label, name = None):
        if name:
            self.name = name
            self.obj_name = '%s_%s'%(label, name)
        else:
            self.obj_name = label
    
    def build_object(self, arrays, attributes = {}):
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis = False)
        #
    
    def load(self):
        flag = True
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            flag = False
        return flag
    
    @property
    def obj(self):
        return self.get_obj()
    
    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.'%self.obj_name)
        return obj

    @property
    def gnodes(self):
        return self.get_gnodes()
    
    def get_gnodes(self):
        name = 'GeometryNodes_%s'%self.obj_name
        modifier = self.obj.modifiers.get(name)
        if modifier is None:
            self.build_geometry_node()
        return modifier

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        name = 'GeometryNodes_%s'%self.obj_name
        modifier = self.obj.modifiers.new(name = name, type = 'NODES')
        modifier.node_group.name = name
    
    def add_geometry_node(self, sp):
        """
        add geometry node for each bond pair
        """
        gn = self.gnodes
    
    def update(self, ):
        """
        update object'data.
        calculate data in all farmes
        """
        pass

    @property
    def arrays(self):
        return self.get_arrays()
    
    @arrays.setter
    def arrays(self, arrays):
        self.set_arrays(arrays)

    def get_arrays(self):
        """
        """
        object_mode()
        arrays = self.attributes
        arrays.update({'positions': self.positions[0],
                        })
        return arrays
    
    @property
    def attributes(self):
        return self.get_attributes()
    
    @attributes.setter
    def attributes(self, attributes):
        self.set_attributes(attributes)

    def get_attributes(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        # attributes
        me = self.obj.data
        nvert = len(me.vertices)
        attributes = {}
        for key in me.attributes.keys():
            att = me.attributes.get(key)
            dtype = att.data_type
            if dtype == 'STRING':
                attributes[key] = np.zeros(nvert, dtype = 'U20')
                for i in range(nvert):
                    attributes[key][i] = att.data[i].value
            elif dtype == 'INT':
                attributes[key] = np.zeros(nvert, dtype = int)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'FLOAT':
                attributes[key] = np.zeros(nvert, dtype = float)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'BOOLEAN':
                attributes[key] = np.zeros(nvert, dtype = bool)
                att.data.foreach_get("value", attributes[key])
            else:
                raise KeyError('%s is not support.'%dtype)
            attributes[key] = np.array(attributes[key])
        return attributes
    
    def set_attributes(self, attributes):
        tstart = time()
        me = self.obj.data
        for key, data in attributes.items():
            # print(key)
            att = me.attributes.get(key)
            if att is None:
                dtype = type(attributes[key][0])
                if np.issubdtype(dtype, int):
                    dtype = 'INT'
                elif np.issubdtype(dtype, float):
                    dtype = 'FLOAT'
                elif np.issubdtype(dtype, str):
                    dtype = 'STRING'
                else:
                    raise KeyError('%s is not supported.'%dtype)
                att = me.attributes.new(name = key, type = dtype, domain = 'POINT')
            if att.data_type == 'STRING':
                nvert = len(me.vertices)
                for i in range(nvert):
                    att.data[i].value = data[i]
            else:
                att.data.foreach_set("value", data)
        me.update()
        # print('set_attributes: %s'%(time() - tstart))

    def set_attribute_with_indices(self, name, indices, data):
        data0 = self.attributes[name]
        data0[indices] = data
        self.set_attributes({name: data0})
    
    @property
    def local_positions(self):
        return self.get_local_positions()
    
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = np.empty(n*3, dtype=np.float64)
        self.obj.data.vertices.foreach_get('co', local_positions)  
        local_positions = local_positions.reshape((n, 3))
        return local_positions
    
    @property
    def positions(self):
        return self.get_positions()
    
    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)
    
    def get_positions(self):
        """
        Get global positions.
        """
        from batoms.tools import local2global
        positions = local2global(self.local_positions, 
                np.array(self.obj.matrix_world))
        return positions
    
    def set_positions(self, positions):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        natom = len(self)
        if len(positions) != natom:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natom))
        positions = local2global(positions, 
                np.array(self.obj.matrix_world), reversed = True)
        # rashpe to (natoms*3, 1) and use forseach_set
        positions = positions.reshape((natom*3, 1))
        # I don't know why 'Basis' shape keys is not updated when editing mesh,
        # so we edit the 'Basis' shape keys directly.
        # self.obj.data.vertices.foreach_set('co', positions)
        self.obj.data.shape_keys.key_blocks[0].data.foreach_set('co', positions)
        self.obj.data.update()
        # bpy.context.view_layer.update()
        # I don't why this is need to update the mesh positions
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
        
    @property    
    def nframe(self):
        return self.get_nframe()
    
    def get_nframe(self):
        if self.obj.data.shape_keys is None:
            return 0
        nframe = len(self.obj.data.shape_keys.key_blocks)
        return nframe
    
    def set_frames_positions(self, frames = None, frame_start = 0, only_basis = False):
        name = ''
        obj = self.obj
        self.set_frames(name, obj, frames, frame_start, only_basis)

    def set_obj_frames(self, name, obj, frames = None, frame_start = 0, only_basis = False):
        """

        frames: list
            list of positions
        
        use shape_keys (faster)
        """
        from batoms.butils import add_keyframe_to_shape_key
        if frames is None:
            frames = self._frames
        centers = frames
        nframe = len(centers)
        if nframe == 0 : return
        sp = ''
        # name = '%s_bond%s'%(self.label, sp)
        # obj = bpy.data.objects.get(name)
        base_name = 'Basis_%s'%(name)
        if obj.data.shape_keys is None:
            obj.shape_key_add(name = base_name)
        elif base_name not in obj.data.shape_keys.key_blocks:
            obj.shape_key_add(name = base_name)
        if only_basis:
            return
        nvert = len(obj.data.vertices)
        for i in range(1, nframe):
            sk = obj.shape_key_add(name = str(i))
            # Use the local position here
            positions = frames[i]
            positions = positions.reshape((nvert*3, 1))
            sk.data.foreach_set('co', positions)
            # Add Keyframes, the last one is different
            if i != nframe - 1:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1, 0], [frame_start + i - 1, 
                    frame_start + i, frame_start + i + 1])
            else:
                add_keyframe_to_shape_key(sk, 'value', 
                    [0, 1], [frame_start + i - 1, frame_start + i])

    def __len__(self):
        return len(self.obj.data.vertices)
    
    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
    
    
    

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

        For example, move H2o molecule by a vector [0, 0, 5]

        >>> h2o.translate([0, 0, 5])
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
        
        >>> h2o.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        angle = angle/180.0*np.pi
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), 
                        orient_type = orient_type)    
    
    def mirror(self, axis = 'Z', orient_type = 'GLOBAL'):
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
        bpy.ops.transform.mirror(constraint_axis = constraint_axis,
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
    
    def delete_obj(self, name):
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
    
