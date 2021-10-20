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
        for p in ['batom', 'bbond', 'bisosurface', 'bpolyhedra']:
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
    objs = []
    for obj in bpy.context.objects_in_mode:
        if obj.batom.flag:
            objs.append(obj)
    dict_mode = {}
    for obj in objs:
        dict_mode[obj.name] = obj.mode
    bpy.ops.object.mode_set(mode='OBJECT')
    selected_vertices = []
    for obj in objs:
        count = len(obj.data.vertices)
        sel = np.zeros(count, dtype=np.bool)
        obj.data.vertices.foreach_get('select', sel)
        selected_vertices.append((obj.name, sel))
    # back to whatever mode we were in
    # bpy.ops.object.mode_set(mode=mode)
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