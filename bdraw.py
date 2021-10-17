import bpy
import numpy as np
from batoms.data import material_styles_dict
import time
#========================================================
def draw_cell_curve(coll, verts, label = None):
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
            p.co = np.append(verts[i], [1.0]) # (add nurbs weight)
        cell = bpy.data.objects.new("cell_%s_edge"%label, crv)
        coll.objects.link(cell)
def draw_cell_cylinder(coll_cell, cell_vertices, label = None, celllinewidth = 0.05):
    """
    Draw unit cell using cylinder.
    """
    if cell_vertices is not None:
        # build materials
        material = bpy.data.materials.new('cell_{0}'.format(label))
        # material.label = 'cell'
        material.diffuse_color = (0.0, 0.0, 0.0, 1.0)
        # draw points
        bpy.ops.mesh.primitive_uv_sphere_add(radius = celllinewidth) #, segments=32, ring_count=16)
        sphere = bpy.context.view_layer.objects.active
        sphere.name = 'instancer_cell_%s_sphere'%label
        sphere.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        sphere.hide_set(True)
        mesh = bpy.data.meshes.new('point_cell' )
        obj_cell = bpy.data.objects.new('cell_%s_point'%label, mesh )
        # Associate the vertices
        obj_cell.data.from_pydata(cell_vertices, [], [])
        sphere.parent = obj_cell
        obj_cell.instance_type = 'VERTS'
        coll_cell.objects.link(sphere)
        coll_cell.objects.link(obj_cell)
        #
        # edges
        edges = [[3, 0], [3, 1], [4, 0], [4, 1],
                    [2, 5], [2, 6], [7, 5], [7, 6], 
                    [3, 2], [0, 6], [1, 5], [4, 7]
            ]
        cell_edges = {'lengths': [], 
                      'centers': [],
                      'normals': []}
        for e in edges:
            center = (cell_vertices[e[0]] + cell_vertices[e[1]])/2.0
            vec = cell_vertices[e[0]] - cell_vertices[e[1]]
            length = np.linalg.norm(vec)
            nvec = vec/length
            # print(center, nvec, length)
            cell_edges['lengths'].append(length)
            cell_edges['centers'].append(center)
            cell_edges['normals'].append(nvec)
        #
        source = bond_source(vertices=6)
        verts, faces = cylinder_mesh_from_instance_vec(cell_edges['centers'], cell_edges['normals'], cell_edges['lengths'], celllinewidth, source)
        # print(verts)
        mesh = bpy.data.meshes.new("cell_cylinder")
        mesh.from_pydata(verts, [], faces)  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_edge = bpy.data.objects.new("cell_%s_cylinder"%label, mesh)
        obj_edge.data = mesh
        obj_edge.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        coll_cell.objects.link(obj_edge)


def draw_text(coll_text = None, atoms = None, type = None):
    tstart = time.time()
    positions = atoms.positions
    n = len(positions)
    for i in range(n):
        location = positions[i] + np.array([0, 0, 1.0])
        FontCurve = bpy.data.curves.new(type="FONT",name="myFontCurve")
        ob = bpy.data.objects.new("myFontOb",FontCurve)
        if type == 0:
            ob.data.body = "%s"%i
        elif type == 1:
            ob.data.body = "%s"%atoms[i].symbol
        ob.location = location
        coll_text.objects.link(ob)
    print('text: {0:10.2f} s'.format(time.time() - tstart))


def draw_bond_kind(kind, 
                   datas, 
                   label = None,
                   coll = None,
                   source = None, 
                   bsdf_inputs = None, 
                   material_style = 'plastic'):
    if len(datas['centers']) == 0:
        return
    if not bsdf_inputs:
        bsdf_inputs = material_styles_dict[material_style]
    if datas['style'] in ['0', '1']:
        vertices = 16
    elif datas['style'] in ['2', '3']:
        vertices = 6
    source = bond_source(vertices = vertices)
    tstart = time.time()
    material = bpy.data.materials.new('bond_kind_{0}'.format(kind))
    material.diffuse_color = np.append(datas['color'], datas['transmit'])
    material.metallic = bsdf_inputs['Metallic']
    material.roughness = bsdf_inputs['Roughness']
    material.blend_method = 'BLEND'
    material.use_nodes = True
    principled_node = material.node_tree.nodes['Principled BSDF']
    principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
    principled_node.inputs['Alpha'].default_value = datas['transmit']
    for key, value in bsdf_inputs.items():
        principled_node.inputs[key].default_value = value
    datas['materials'] = material
    #
    verts, faces = cylinder_mesh_from_instance_vec(datas['centers'], datas['normals'], datas['lengths'], datas['width'], source)
    mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
    mesh.from_pydata(verts, [], faces)  
    mesh.update()
    for f in mesh.polygons:
        f.use_smooth = True
    obj_bond = bpy.data.objects.new("bond_{0}_{1}".format(label, kind), mesh)
    obj_bond.data = mesh
    obj_bond.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    coll.objects.link(obj_bond)
    # print('bonds: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))
    
def draw_cavity(kind, 
                   datas, 
                   label = None,
                   coll = None,
                   source = None, 
                   bsdf_inputs = None, 
                   material_style = 'plastic'):
    if len(datas['edges']) == 0:
        return
    if not bsdf_inputs:
        bsdf_inputs = material_styles_dict[material_style]
    tstart = time.time()
    material = bpy.data.materials.new('bond_kind_{0}'.format(kind))
    material.diffuse_color = np.append(datas['color'], datas['transmit'])
    material.metallic = bsdf_inputs['Metallic']
    material.roughness = bsdf_inputs['Roughness']
    material.blend_method = 'BLEND'
    material.use_nodes = True
    principled_node = material.node_tree.nodes['Principled BSDF']
    principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
    principled_node.inputs['Alpha'].default_value = datas['transmit']
    for key, value in bsdf_inputs.items():
        principled_node.inputs[key].default_value = value
    datas['materials'] = material
    #
    mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
    mesh.from_pydata(datas['vertices'], datas['edges'], [])
    mesh.update()
    name = "cavity_{0}_{1}".format(label, kind)
    obj_bond = bpy.data.objects.new(name, mesh)
    obj_bond.data = mesh
    obj_bond.data.materials.append(material)
    coll.objects.link(obj_bond)
    bpy.context.view_layer.objects.active = bpy.context.scene.objects.get(name)
    bpy.ops.object.mode_set(mode = 'EDIT')
    # fill edge with face
    bpy.ops.mesh.edge_face_add()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.shade_smooth()
    # print('bonds: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))
    


def draw_polyhedra_kind(kind, 
                        datas, 
                        label = None,
                        coll = None,
                        source = None, 
                        show_edge = True,
                        bsdf_inputs = None, 
                        material_style = 'default'):
        """
        """
        if not bsdf_inputs:
            bsdf_inputs = material_styles_dict[material_style]
        tstart = time.time()
        material = bpy.data.materials.new('polyhedra_kind_{0}'.format(kind))
        material.diffuse_color = np.append(datas['color'], datas['transmit'])
        material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
        principled_node.inputs['Alpha'].default_value = datas['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['materials'] = material
        #
        # create new mesh structure
        mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
        # mesh.from_pydata(datas['vertices'], datas['edges'], datas['faces'])  
        mesh.from_pydata(datas['vertices'], [], datas['faces'])  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = False
        obj_polyhedra = bpy.data.objects.new("polyhedra_{0}_{1}_face".format(label, kind), mesh)
        obj_polyhedra.data = mesh
        obj_polyhedra.data.materials.append(material)
        # bpy.ops.object.shade_smooth()
        # bpy.ops.object.shade_flat()
        coll.objects.link(obj_polyhedra)
        #---------------------------------------------------
        if not show_edge: return
        source = bond_source(vertices=6)
        material = bpy.data.materials.new('polyhedra_edge_kind_{0}'.format(kind))
        material.diffuse_color = np.append(datas['edge_cylinder']['color'], datas['edge_cylinder']['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['edge_cylinder']['color'], datas['edge_cylinder']['transmit'])
        principled_node.inputs['Alpha'].default_value = datas['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        verts, faces = cylinder_mesh_from_instance_vec(datas['edge_cylinder']['centers'], datas['edge_cylinder']['normals'], datas['edge_cylinder']['lengths'], datas['edgewidth'], source)
        # print(verts)
        mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
        mesh.from_pydata(verts, [], faces)  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_edge = bpy.data.objects.new("polyhedra_{0}_{1}_edge".format(label, kind), mesh)
        obj_edge.data = mesh
        obj_edge.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        # STRUCTURE.append(obj_polyhedra)
        coll.objects.link(obj_edge)
        # print('polyhedras: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))


def draw_isosurface(coll_isosurface, verts, faces, color,
                    bsdf_inputs = None, material_style = 'default'):
    """Computes an isosurface from a volume grid.
    
    Parameters:

    
    """
    #material
    if not bsdf_inputs:
        bsdf_inputs = material_styles_dict[material_style]
    material = bpy.data.materials.new('isosurface')
    material.name = 'isosurface'
    material.diffuse_color = color
    material.blend_method = 'BLEND'
    material.use_nodes = True
    principled_node = material.node_tree.nodes['Principled BSDF']
    principled_node.inputs['Base Color'].default_value = color
    principled_node.inputs['Alpha'].default_value = color[3]
    for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
    # create new mesh structure
    isosurface = bpy.data.meshes.new("isosurface")
    isosurface.from_pydata(verts, [], faces)  
    isosurface.update()
    for f in isosurface.polygons:
        f.use_smooth = True
    iso_object = bpy.data.objects.new("isosurface", isosurface)
    iso_object.data = isosurface
    iso_object.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    coll_isosurface.objects.link(iso_object)


# draw bonds
def bond_source(vertices = 12, depth = 1.0):
    bpy.ops.mesh.primitive_cylinder_add(vertices = vertices, depth = depth)
    cyli = bpy.context.view_layer.objects.active
    me = cyli.data
    verts = []
    faces = []
    for vertices in me.vertices:
        verts.append(np.array(vertices.co))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
            # print("    Vertex: %d" % me.loops[loop_index].vertex_index)
            face.append(me.loops[loop_index].vertex_index)
        faces.append(face)
    cyli.select_set(True)
    bpy.data.objects.remove(cyli, do_unlink = True)
    n = len(faces[0])
    faces1 = [faces[i] for i in range(len(faces)) if len(faces[i]) == n]
    faces2 = [faces[i] for i in range(len(faces)) if len(faces[i]) != n]
    return verts, faces1, faces2
# draw atoms
def atom_source():
    bpy.ops.mesh.primitive_uv_sphere_add() #, segments=32, ring_count=16)
    # bpy.ops.mesh.primitive_cylinder_add()
    sphe = bpy.context.view_layer.objects.active
    me = sphe.data
    verts = []
    faces = []
    for vertices in me.vertices:
        verts.append(np.array(vertices.co))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
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
            face = [x+ i*nvert for x in face]
            faces.append(face)
    return verts, faces

def cylinder_mesh_from_instance_vec(centers, normals, lengths, scale, source):
    from scipy.spatial.transform import Rotation as R
    tstart = time.time()
    vert0, face1, face2 = source
    nb = len(centers)
    nvert = len(vert0)
    centers = np.array(centers)
    normals = np.array(normals)
    vert0 = np.array(vert0)
    face1 = np.array(face1)
    face2 = np.array(face2)
    scale = np.array([[scale, scale]]*len(centers))
    scale = np.append(scale, np.array(lengths).reshape(-1, 1), axis = 1)
    # print(np.cross([0.0000014159, 0.000001951, 1], normals[0]))
    vec = np.cross([0.0, 0.0, 1], normals) + np.array([0.000000001, 0, 0])
    vec = vec/np.linalg.norm(vec, axis = 1)[:, None]
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
    #----------------------------
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
    faces = list(faces1) + list(face2)
    # print('cylinder_mesh_from_instance: {0:10.2f} s'.format( time.time() - tstart))
    return verts, faces

def draw_plane(location = (0, 0, -0.01), color = (0.2, 0.2, 1.0, 1.0), size = 200, bsdf_inputs = None, material_style = 'default'):
    """
    Draw a plane.
    location: array
    color: array
    size: float
    bsdf_inputs: dict
    material_style: str
    """
    # build materials
    if not bsdf_inputs:
        bsdf_inputs = material_styles_dict[material_style]
    material = bpy.data.materials.new('plane')
    material.name = 'plane'
    material.diffuse_color = color
    # material.blend_method = 'BLEND'
    material.use_nodes = True
    principled_node = material.node_tree.nodes['Principled BSDF']
    principled_node.inputs['Alpha'].default_value = color[3]
    for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
    # Instantiate a floor plane
    bpy.ops.mesh.primitive_plane_add(size=size, location=location)
    current_object = bpy.context.object
    current_object.data.materials.append(material)