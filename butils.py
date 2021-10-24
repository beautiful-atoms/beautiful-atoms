import bpy
from ase import Atom, Atoms



def object_mode():
    for object in bpy.data.objects:
            if object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode = 'OBJECT')

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
            if getattr(obj, p).flag:
                batoms_list.append(getattr(obj, p).label)
    batoms_list = list(set(batoms_list))
    return batoms_list

def get_selected_objects(attr):
    """
    """
    bond_list = []
    for obj in bpy.context.selected_objects:
        if getattr(obj, attr).flag:
            bond_list.append(obj.name)
    return bond_list

def get_selected_vertices():
    """
    change to 'Object' mode
    in order
    """
    import numpy as np
    selected_vertices = []
    objs = []
    for obj in bpy.context.objects_in_mode:
        if obj.batom.flag:
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
            selected_vertices.append((obj.name, index))
    # back to whatever mode we were in
    bpy.ops.object.mode_set(mode='EDIT')
    return selected_vertices

def remove_collection(name):
    """
    """
    collection = bpy.data.collections.get(name)
    objs = collection.all_objects.keys()
    for obj in objs:
        obj = bpy.data.objects.get(obj)
        bpy.data.objects.remove(obj, do_unlink=True)
    collection = bpy.data.collections.get(name)
    colls = collection.children.keys()
    for coll in colls:
        coll = bpy.data.collections.get(coll)
        bpy.data.collections.remove(coll)
    collection = bpy.data.collections.get(name)
    bpy.data.collections.remove(collection)
        
def clean_objects():
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)

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

def lock_camera_to_view(switch):
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.lock_camera = switch

def get_keyframes_of_batoms(batoms):
    """
    get keyframes of a batoms
    """
    keyframes = []
    for sp, ba in batoms.items():
        anim = ba.batom.data.animation_data
        if anim is not None and anim.action is not None:
            for fcu in anim.action.fcurves:
                for keyframe in fcu.keyframe_points:
                    x, y = keyframe.co
                    if x not in keyframes:
                        keyframes.append((int(x)))
    return keyframes