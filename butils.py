import bpy
from ase import Atom, Atoms
from mathutils import Vector, Matrix
import numpy as np


def object_mode():
    for object in bpy.data.objects:
            if object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode = 'OBJECT')
def eidt_mode():
    try:
        bpy.ops.object.mode_set(mode = 'EDIT')
    except:
        pass

def read_batoms_list():
    """
    Read all batoms collection 
    """
    items = [col.name for col in bpy.data.collections if col.batoms.flag]
    return items

def get_selected_batoms():
    """
    """
    batoms_list = []
    for obj in bpy.context.selected_objects:
        for p in ['batom', 'bbond', 'bisosurface', 'bpolyhedra', 'bplane']:
            if getattr(obj.batoms, p).flag:
                batoms_list.append(getattr(obj.batoms, p).label)
    batoms_list = list(set(batoms_list))
    return batoms_list

def get_selected_objects(attr):
    """
    """
    objs = []
    for obj in bpy.context.selected_objects:
        if getattr(obj.batoms, attr).flag:
            objs.append(obj.name)
    return objs

def get_selected_vertices():
    """
    change to 'Object' mode
    in order
    """
    import numpy as np
    eidt_mode()
    selected_vertices = []
    objs = []
    for obj in bpy.context.objects_in_mode:
        if obj.batoms.batom.flag:
            objs.append(obj)
    if len(objs) == 0:
        print('Warning: No atoms is selected, please switch to Edit mode.')
        return selected_vertices
    dict_mode = {}
    for obj in objs:
        dict_mode[obj.name] = obj.mode
    bpy.ops.object.mode_set(mode='OBJECT')
    for obj in objs:
        count = len(obj.data.vertices)
        sel = np.zeros(count, dtype=np.bool)
        obj.data.vertices.foreach_get('select', sel)
        index = np.where(sel)[0]
        if len(index) > 0:
            selected_vertices.append((obj.batoms.batom.label, obj.batoms.batom.species, obj.name, index))
    # back to whatever mode we were in
    bpy.ops.object.mode_set(mode='EDIT')
    return selected_vertices

def remove_collection(name, keep_batom = True):
    """
    """
    collection = bpy.data.collections.get(name)
    objs = collection.all_objects.keys()
    for obj in objs:
        obj = bpy.data.objects.get(obj)
        if keep_batom and obj.batoms.batom: continue
        bpy.data.objects.remove(obj, do_unlink=True)
    collection = bpy.data.collections.get(name)
    colls = collection.children.keys()
    for coll in colls:
        coll = bpy.data.collections.get(coll)
        bpy.data.collections.remove(coll)
    collection = bpy.data.collections.get(name)
    bpy.data.collections.remove(collection)
        
def clean_objects_by_name(name):
    obj = bpy.data.objects.get(name)
    bpy.data.objects.remove(obj)

def clean_coll_objects(coll, names = None):
    """
    remove all bond object in the bond collection
    """
    if not names:
        for obj in coll.all_objects:
            bpy.data.objects.remove(obj, do_unlink = True)
    else:
        for name in names:
            for obj in coll.all_objects:
                if name in obj.name:
                    bpy.data.objects.remove(obj, do_unlink = True)    

def removeAll():
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)
    for obj in bpy.data.objects:
        bpy.data.objects.remove(obj)
    for cam in bpy.data.cameras:
        bpy.data.cameras.remove(cam)
    for light in bpy.data.lights:
        bpy.data.lights.remove(light)
    for mat in bpy.data.materials:
        bpy.data.materials.remove(mat)
    for coll in bpy.data.collections:
        if coll.name == 'Collection': continue
        bpy.data.collections.remove(coll)

def get_keyframes_of_batoms(batoms):
    """
    get keyframes of a batoms
    """
    keyframes = []
    for sp, ba in batoms.items():
        anim = ba.obj.data.animation_data
        if anim is not None and anim.action is not None:
            for fcu in anim.action.fcurves:
                for keyframe in fcu.keyframe_points:
                    x, y = keyframe.co
                    if x not in keyframes:
                        keyframes.append((int(x)))
    return keyframes
def get_shape_key_of_batoms(batoms):
    """
    get shape keys of a batoms
    """
    keyframes = []
    for sp, ba in batoms.items():
        anim = ba.obj.data.animation_data
        if anim is not None and anim.action is not None:
            for fcu in anim.action.fcurves:
                for keyframe in fcu.keyframe_points:
                    x, y = keyframe.co
                    if x not in keyframes:
                        keyframes.append((int(x)))
    return keyframes
def add_keyframe_to_shape_key(obj, attr, values, frames):
    """
    """
    for value, frame in zip(values, frames):
        setattr(obj, attr, value)
        obj.keyframe_insert(data_path = attr, frame = frame)

def add_keyframe_visibility(obj, nframe, frame):
    # key as visible on the current frame
    obj.animation_data_create()
    ac = bpy.data.actions.new('Visibility Action')
    obj.animation_data.action = ac
    for data_path in ['hide_viewport', 'hide_render']:
        fc = ac.fcurves.new(data_path=data_path)
        fc.keyframe_points.add(nframe)
        value = [0]*2*nframe
        for i in range(nframe):
            value[2*i] = i
            if i == frame:
                value[2*i + 1] = 0
            else:
                value[2*i + 1] = 1
        fc.keyframe_points.foreach_set('co', value)

def set_look_at(obj, target, roll=0):
        """
        Rotate obj to look at target
        """
        if not isinstance(target, Vector):
            target = Vector(target)
        loc = obj.location
        direction = target - loc
        quat = direction.to_track_quat('-Z', 'Y')
        quat = quat.to_matrix().to_4x4()
        rollMatrix = Matrix.Rotation(roll, 4, 'Z')
        loc = loc.to_tuple()
        obj.matrix_world = quat @ rollMatrix
        obj.location = loc

def lock_to(obj, target = None, location = True, rotation = True):
        """
        track to obj
        """
        if target is not None:
            if location:
                obj.constraints.new(type = 'COPY_LOCATION')
                obj.constraints["Copy Location"].target = target
            if rotation:
                obj.constraints.new(type = 'COPY_ROTATION')
                obj.constraints["Copy Rotation"].target = target
        else:
            for c in obj.constraints:
                obj.constraints.remove(c)

def set_world(color = [0.2, 0.2, 0.2, 1.0]):
        """
        """
        world = bpy.data.scenes['Scene'].world
        world.use_nodes = True
        node_tree = world.node_tree
        node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
        node_tree.nodes["Background"].inputs["Color"].default_value = color

def get_nodes_by_name(nodes, name, type = None):
    node = nodes.get(name)
    if node is None:
        node = nodes.new(type)
        node.name = name
    return node

def clean_default(camera = False, light = True):
        if 'Cube' in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects["Cube"], do_unlink=True)
        if camera and 'Camera' in bpy.data.cameras:
            bpy.data.cameras.remove(bpy.data.cameras['Camera'])
        if light and 'Light' in bpy.data.lights:
            bpy.data.lights.remove(bpy.data.lights['Light'])

def show_index():
        """
        show index of vertices in Edit mode
        """
        bpy.context.preferences.view.show_developer_ui = True
        for a in bpy.context.screen.areas:
            if a.type == 'VIEW_3D':
                overlay = a.spaces.active.overlay
                overlay.show_extra_indices = True
    
def get_area(me):
    npoly = len(me.polygons)
    areas = np.zeros(npoly)
    me.polygons.foreach_get('area', areas)
    total_area = np.sum(areas)
    return total_area

def get_volume(me):
    import bmesh
    bm = bmesh.new()
    bm.from_mesh(me)
    volume = bm.calc_volume()
    return volume