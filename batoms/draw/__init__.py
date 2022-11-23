import bpy
from batoms.material import create_material
import numpy as np
from time import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)
# ========================================================


def draw_cell_curve(coll, verts, label=None):
    """
    Draw unit cell by edge, however, can not be rendered.
    """
    if verts is not None:
        edges = [0, 4, 6, 2, 0, 1, 5, 7, 3, 1]
        crv = bpy.data.curves.new("edge_cell", 'CURVE')
        crv.dimensions = '3D'
        spline = crv.splines.new(type='NURBS')
        spline.points.add(len(edges)-1)
        for p, i in zip(spline.points, edges):
            p.co = np.append(verts[i], [1.0])  # (add nurbs weight)
        cell = bpy.data.objects.new("cell_%s_edge" % label, crv)
        coll.objects.link(cell)


def draw_cylinder(
        name=None,
        datas=[],
        coll=None,
        node_type='Principled BSDF',
        node_inputs=None,
        material_style='plastic'):
    """
    Draw cylinder.
    """
    from batoms.data.source_data import bond_source
    # materials
    material = create_material(name,
                               datas['color'],
                               node_type=node_type,
                               node_inputs=node_inputs,
                               material_style=material_style)
    #
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    obj.data.materials.append(material)
    if len(datas['centers']) == 0:
        return obj

    source = bond_source[datas['vertices']]
    # tstart = time()
    verts, faces = cylinder_mesh_from_vec(
        datas['centers'], datas['normals'],
        datas['lengths'], datas['width'], source)
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.polygons.foreach_set('use_smooth', [True]*len(mesh.polygons))
    obj.data = mesh
    obj.data.materials.append(material)
    #
    for name, inputs in datas['battr_inputs'].items():
        battr = getattr(obj, name)
        for key, value in inputs.items():
            setattr(battr, key, value)
    # bpy.ops.object.shade_smooth()
    coll.objects.link(obj)
    return obj


def draw_surface_from_vertices(name,
                               datas,
                               coll=None,
                               use_smooth=True):
    if len(datas['vertices']) == 0:
        return
    #
    mesh = bpy.data.meshes.new(name)
    if len(datas['edges']) > 0:
        mesh.from_pydata(datas['vertices'], datas['edges'], datas['faces'])
    else:
        mesh.from_pydata(datas['vertices'], [], datas['faces'])
    mesh.update()
    mesh.polygons.foreach_set('use_smooth',
                              [use_smooth]*len(mesh.polygons))
    obj = bpy.data.objects.new(name, mesh)
    obj.data = mesh
    #
    for name, inputs in datas['battr_inputs'].items():
        battr = getattr(obj, name)
        for key, value in inputs.items():
            setattr(battr, key, value)
    bpy.ops.object.shade_smooth()
    if coll is not None:
        coll.objects.link(obj)
    return obj
    # print('bonds: {0}   {1:10.2f} s'.format(name, time() - tstart))


def draw_vertices(name, vertices):
    datas = {'vertices': vertices,
             'edges': [],
             'faces': [],
             'color': [0.5, 0.5, 0.5, 1.0],
             'battr_inputs': {}}
    coll = bpy.data.collections['Collection']
    draw_surface_from_vertices(name, datas, coll)


def draw_text(coll_text=None, atoms=None, type=None):
    tstart = time()
    positions = atoms.positions
    n = len(positions)
    for i in range(n):
        location = positions[i] + np.array([0, 0, 1.0])
        FontCurve = bpy.data.curves.new(type="FONT", name="myFontCurve")
        ob = bpy.data.objects.new("myFontOb", FontCurve)
        if type == 0:
            ob.data.body = "%s" % i
        elif type == 1:
            ob.data.body = "%s" % atoms[i].symbol
        ob.location = location
        coll_text.objects.link(ob)
    logger.debug('text: {0:10.2f} s'.format(time() - tstart))


def draw_2d_slicing(name,
                    datas,
                    coll=None
                    ):
    """
    using texture image
    """
    from bpy_extras.image_utils import load_image
    image = load_image(datas['imagename'], check_existing=True)
    bpy.ops.mesh.primitive_plane_add(size=1, location=datas['location'],
                                     rotation=datas['rotation'],
                                     )
    plane = bpy.context.active_object
    bpy.ops.object.mode_set(mode='OBJECT')
    plane.scale = datas['size']
    plane.data.name = name
    plane.name = name
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    material.blend_method = 'BLEND'
    node_tree = material.node_tree
    material_output = node_tree.nodes.get("Material Output")
    core_shader = node_tree.nodes.get('Principled BSDF')
    node_texture = node_tree.nodes.new('ShaderNodeTexImage')
    node_texture.image = image
    node_texture.show_texture = True
    node_tree.links.new(core_shader.inputs[0],
                        node_texture.outputs['Color'])
    node_tree.links.new(
        core_shader.inputs['Alpha'], node_texture.outputs['Alpha'])
    node_tree.links.new(
        material_output.inputs['Surface'], core_shader.outputs[0])
    plane.data.materials.append(material)
    coll.objects.link(plane)


def clean_node_tree(node_tree):
    """Clear all nodes in a shader node tree except the output.

    Returns the output node
    """
    nodes = node_tree.nodes
    # copy to avoid altering the loop's data source
    for node in list(nodes):
        if not node.type == 'OUTPUT_MATERIAL':
            nodes.remove(node)

    return node_tree.nodes[0]


# draw bonds
def bond_source(vertices=12, depth=1.0):
    bpy.ops.mesh.primitive_cylinder_add(vertices=vertices, depth=depth)
    cyli = bpy.context.view_layer.objects.active
    me = cyli.data
    faces = []
    n = len(me.vertices)
    vertices = np.empty(n*3, dtype=np.float64)
    me.vertices.foreach_get('co', vertices)
    vertices = vertices.reshape((n, 3))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start,
                                poly.loop_start + poly.loop_total):
            # print("    Vertex: %d" % me.loops[loop_index].vertex_index)
            face.append(me.loops[loop_index].vertex_index)
        faces.append(face)
    cyli.select_set(True)
    bpy.data.objects.remove(cyli, do_unlink=True)
    n = len(faces[0])
    faces1 = [faces[i] for i in range(len(faces)) if len(faces[i]) == n]
    faces2 = [faces[i] for i in range(len(faces)) if len(faces[i]) != n]
    return vertices, faces1, faces2

# draw bonds


def cylinder2(vertices=16, depth=1.0):
    bpy.ops.mesh.primitive_cylinder_add(vertices=32, depth=1)
    obj = bpy.context.view_layer.objects.active
    me = obj.data
    # select edges for subdivde
    n = len(me.vertices)
    vertices = np.zeros(n*3, dtype=np.float64)
    me.vertices.foreach_get('co', vertices)
    vertices = vertices.reshape((n, 3))
    #
    me.update()
    m = len(me.edges)
    selects = np.zeros(m, dtype=bool)
    for i in range(m):
        v0 = me.edges[i].vertices[0]
        v1 = me.edges[i].vertices[1]
        center = (vertices[v0] + vertices[v1])/2
        # print(i, v0, v1, center)
        if np.isclose(center[2], 0):
            selects[i] = True

    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    me.edges.foreach_set('select', selects)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.subdivide(number_cuts=1, smoothness=0,
                           fractal_along_normal=0)
    bpy.ops.object.mode_set(mode='OBJECT')
    #
    m = len(me.polygons)
    centers = np.zeros(m*3, dtype=np.float64)
    me.polygons.foreach_get('center', centers)
    centers = centers.reshape((m, 3))
    return obj


def atom_source():
    bpy.ops.mesh.primitive_uv_sphere_add()  # , segments=32, ring_count=16)
    # bpy.ops.mesh.primitive_cylinder_add()
    sphe = bpy.context.view_layer.objects.active
    me = sphe.data
    verts = []
    faces = []
    for vertices in me.vertices:
        verts.append(np.array(vertices.co))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start,
                                poly.loop_start + poly.loop_total):
            # print("    Vertex: %d" % me.loops[loop_index].vertex_index)
            face.append(me.loops[loop_index].vertex_index)
        faces.append(face)
    sphe.select_set(True)
    bpy.ops.object.delete()
    return [verts, faces]


def sphere_mesh_from_instance(centers, radius, source):
    verts = []
    faces = []
    vert0, face0 = source
    nvert = len(vert0)
    nb = len(centers)
    for i in range(nb):
        center = centers[i]
        # r = radius[i]
        # normal = normal/np.linalg.norm(normal)
        for vert in vert0:
            vert = vert*[radius, radius, radius]
            vert += center
            verts.append(vert)
        for face in face0:
            face = [x + i*nvert for x in face]
            faces.append(face)
    return verts, faces


def cylinder_mesh_from_vec(centers, normals, lengths, scale, source):
    from scipy.spatial.transform import Rotation as R
    tstart = time()
    vert0, face1, face2 = source
    nb = len(centers)
    nvert = len(vert0)
    centers = np.array(centers)
    normals = np.array(normals)
    vert0 = np.array(vert0)
    face1 = np.array(face1)
    face2 = np.array(face2)
    scale = np.array([[scale, scale]]*len(centers))
    scale = np.append(scale, np.array(lengths).reshape(-1, 1), axis=1)
    # print(np.cross([0.0000014159, 0.000001951, 1], normals[0]))
    vec = np.cross([0.0, 0.0, 1], normals) + np.array([1e-6, 0, 0])
    vec = vec/np.linalg.norm(vec, axis=1)[:, None]
    # print(np.arccos(normals[0, 2]*0.999999))
    ang = np.arccos(normals[:, 2]*0.999999)
    # print(-1*ang[0]*vec[0])
    vec = -1*(vec.T*ang).T
    # print(R.from_rotvec(vec[0]).as_matrix())
    r = R.from_rotvec(vec)
    matrix = r.as_matrix()
    # print('vert1 0: ', vert0*scale[0])
    verts = np.tile(vert0, (nb, 1, 1))
    verts = verts*scale[:, None]
    # print('matrix: ', verts[0].dot(matrix[0]))
    verts = np.matmul(verts, matrix)
    # print('matrix: ', verts)
    # print(verts.shape)
    # centers = np.tile(centers, (nb, 1))
    # print('center: ', verts[0] + centers[0])
    verts += centers[:, None]
    # print('center: ', verts)
    # faces.append(face)
    verts = verts.reshape(-1, 3)
    # ----------------------------
    nf1 = len(face1[0])
    nf2 = len(face2[0])
    faces1 = np.tile(face1, (nb, 1, 1))
    faces2 = np.tile(face2, (nb, 1, 1))
    offset = np.arange(nb)*nvert
    offset = offset.reshape(-1, 1, 1)
    faces1 = faces1 + offset
    faces1 = faces1.reshape(-1, nf1)
    faces2 = faces2 + offset
    faces2 = faces2.reshape(-1, nf2)
    faces = list(faces1) + list(faces2)
    return verts, faces


def draw_plane(name='plane',
               location=(0, 0, -0.01),
               rotation=(0, 0, 0),
               color=(0.2, 0.2, 1.0, 1.0),
               size=200,
               node_type=None,
               node_inputs=None,
               material_style='default'):
    """
    Draw a plane.
    location: array
    color: array
    size: float
    node_inputs: dict
    material_style: str
    """
    # build materials
    material = create_material(name,
                               color,
                               node_type=node_type,
                               node_inputs=node_inputs,
                               material_style=material_style)
    # Instantiate a floor plane
    bpy.ops.mesh.primitive_plane_add(
        size=size, location=location, rotation=rotation)
    obj = bpy.context.object
    obj.data.materials.append(material)
    return obj
