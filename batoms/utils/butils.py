from shutil import ExecError
import bpy
from mathutils import Vector, Matrix
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

def get_selected_vertices_bmesh(obj):
    """_summary_

    Warnning: bmesh.select_history does not support
     box/lasso/circle/other selected vertices.

    Args:
        obj (bpy.type.object): _description_

    Returns:
        _type_: _description_
    """
    import bmesh
    if obj.batoms.type == 'BATOMS' and obj.mode == 'EDIT':
        data = obj.data
        bm = bmesh.from_edit_mesh(data)
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        return v
    else:
        return []

def get_all_consoles():
    """Get all consoles

    Returns:
        list: a list of consoles
    """
    from console_python import get_console
    consoles = []
    for area in bpy.context.screen.areas:
        if area.type == 'CONSOLE':
            # areas.append(area)
            for region in area.regions:
                if region.type == 'WINDOW':
                    console, stdout, stderr = get_console(hash(region))
                    consoles.append(console)
    return consoles

def object_mode():
    for object in bpy.data.objects:
        if object.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')


def eidt_mode():
    try:
        bpy.ops.object.mode_set(mode='EDIT')
    except ExecError:
        logger.critical('Can not change to EDIT mode')


def read_batoms_list():
    """
    Read all batoms collection
    """
    items = [col.name for col in bpy.data.collections if col.batoms.type != 'OTHER']
    return items


def get_selected_batoms():
    """
    """
    batoms_list = []
    for obj in bpy.context.selected_objects:
        # for p in ['batom', 'bbond', 'bisosurface', 'bpolyhedra', 'bplane']:
            if obj.batoms.type != 'OTHER':
                batoms_list.append(obj.batoms.label)
    batoms_list = list(set(batoms_list))
    # print(batoms_list)
    return batoms_list


def get_selected_objects(attr):
    """
    """
    objs = []
    for obj in bpy.context.selected_objects:
        if obj.batoms.type == attr.upper():
            objs.append(obj.name)
    return objs

def get_selected_vertices(obj):
    """
    """
    import numpy as np
    selected_vertices = []
    # if obj.mode != "EDIT":
        # logger.warning('Warning: Please switch to Edit mode.')
        # return selected_vertices
    mode = obj.mode
    bpy.ops.object.mode_set(mode='OBJECT')
    count = len(obj.data.vertices)
    sel = np.zeros(count, dtype=np.bool)
    obj.data.vertices.foreach_get('select', sel)
    selected_vertices = np.where(sel)[0]
    # back to whatever mode we were in
    bpy.ops.object.mode_set(mode=mode)
    return selected_vertices.tolist()

def get_selected_vertices_all():
    """
    change to 'Object' mode
    in order
    """
    import numpy as np
    eidt_mode()
    selected_vertices = []
    objs = []
    for obj in bpy.context.objects_in_mode:
        if obj.batoms.type != 'OTHER':
            objs.append(obj)
    if len(objs) == 0:
        logger.warning('Warning: No atoms is selected, please switch to Edit mode.')
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
            selected_vertices.append((obj.batoms.label, index))
    # back to whatever mode we were in
    bpy.ops.object.mode_set(mode='EDIT')
    return selected_vertices


def remove_collection(name, keep_batom=True):
    """

    Note: to avoid crash in macOS, read the name first,
    and then use get method
    """
    collection = bpy.data.collections.get(name)
    objs = collection.all_objects.keys()
    for obj in objs:
        obj = bpy.data.objects.get(obj)
        if keep_batom and obj.batoms.type != "OTHER":
            continue
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


def clean_coll_objects(coll, names=None):
    """
    remove all bond object in the bond collection
    """
    if not names:
        for obj in coll.all_objects:
            bpy.data.objects.remove(obj, do_unlink=True)
    else:
        for name in names:
            for obj in coll.all_objects:
                if name in obj.name:
                    bpy.data.objects.remove(obj, do_unlink=True)


def clean_coll_object_by_type(coll, type):
    """
    remove all objects which have certain type in the given collection.
    """
    objs_name = []
    for obj in coll.all_objects:
        if obj.batoms.type == type:
            objs_name.append(obj.name)
    for obj in objs_name:
        obj = bpy.data.objects.get(obj)
        bpy.data.objects.remove(obj, do_unlink=True)


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
        if coll.name == 'Collection':
            continue
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
        obj.keyframe_insert(data_path=attr, frame=frame)


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


def lock_to(obj, target=None, location=True, rotation=True):
    """
    track to obj
    """
    if target is not None:
        if location:
            obj.constraints.new(type='COPY_LOCATION')
            obj.constraints["Copy Location"].target = target
        if rotation:
            obj.constraints.new(type='COPY_ROTATION')
            obj.constraints["Copy Rotation"].target = target
    else:
        for c in obj.constraints:
            obj.constraints.remove(c)


def set_world(color=[0.2, 0.2, 0.2, 1.0]):
    """
    """
    world = bpy.context.scene.world
    world.use_nodes = True
    node_tree = world.node_tree
    node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
    node_tree.nodes["Background"].inputs["Color"].default_value = color


def get_nodes_by_name(nodes, name, type=None):
    node = nodes.get(name)
    if node is None:
        node = nodes.new(type)
        node.name = name
    return node


def clean_default(camera=False, light=True):
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


def update_object(obj):
    mode = obj.mode
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.mode_set(mode=mode)


def hideOneLevel():
    screen = bpy.context.screen
    outliners = [a for a in screen.areas if a.type == 'OUTLINER']
    c = bpy.context.copy()
    for ol in outliners:
        c["area"] = ol
        ol.tag_redraw()
        bpy.ops.outliner.show_one_level(c, open=False)
        ol.tag_redraw()

def build_modifier(obj, name):
    from bl_operators.geometry_nodes import geometry_node_group_empty_new
    modifier = obj.modifiers.new(name=name, type='NODES')
    if bpy.app.version_string >= '3.2.0':
        # bpy.context.view_layer.objects.active = obj
        # bpy.context.object.modifiers.active = modifier
        # bpy.ops.node.new_geometry_node_group_assign()
        group = geometry_node_group_empty_new()
        modifier.node_group = group
    modifier.node_group.name = name
    return modifier

def get_att_length(mesh, att):
    """get attribute length based on domain

    Args:
        att (bpy.types.Attribute): attribute of a mesh
    Returns:
        int: length of the attribute
    """
    domain = att.domain
    if domain == 'POINT':
        n = len(mesh.vertices)
    elif domain == 'EDGE':
        n = len(mesh.edges)
    else:
        n = len(mesh.polygons)
    return n

def get_bmesh_domain(bm, att):
    """Get bmesh domain

    Args:
        bm (_type_): _description_
        att (_type_): _description_
    """
    if att.domain == 'POINT':
        domain = getattr(bm, "verts")
    elif att.domain == 'EDGE':
        domain = getattr(bm, "edges")
    elif att.domain == 'FACE':
        domain = getattr(bm, "faces")
    return domain

def get_bmesh_layer(domain, key, dtype):
    """Get bmesh layer

    Args:
        domain (bm.verts, bm.edges, bm.faces): _description_
        key (str): _description_
        dtype (str): _description_
    """
    if dtype == 'STRING':
        layer = domain.layers.string.get(key)
    elif dtype == 'INT':
        layer = domain.layers.int.get(key)
    elif dtype == 'FLOAT':
        layer = domain.layers.float.get(key)

    return layer


def set_vertex_color(obj, name, color):
    """Set vertex color

    Args:
        obj (bpy.types.Object): Object to be set
        name (str): name of the color attribute
        color (array): value of color for each vertex
    """
    npoint = len(color)
    mesh = obj.data
    if bpy.app.version_string >= '3.2.0':
        color = color.reshape((npoint*4, 1))
        mesh.color_attributes.new(name, 'FLOAT_COLOR', 'POINT')
        mesh.color_attributes[name].data.foreach_set('color', color)
    else:
        import bmesh
        if obj.mode == 'EDIT':
            bm =bmesh.from_edit_mesh(obj.data)
            volume_layer = bm.loops.layers.color.new(name)
            for v in bm.verts:
                for loop in v.link_loops:
                    loop[volume_layer] = color[v.index]
            bmesh.update_edit_mesh(obj.data)
        else:
            bm = bmesh.new()
            bm.from_mesh(obj.data)
            volume_layer = bm.loops.layers.color.new(name)
            for v in bm.verts:
                for loop in v.link_loops:
                    loop[volume_layer] = color[v.index]
            bm.to_mesh(obj.data)
            bm.free()

# ========================================================
if bpy.app.version_string >= '3.1.0':
    compareNodeType = 'FunctionNodeCompare'
else:
    compareNodeType = 'FunctionNodeCompareFloats'
