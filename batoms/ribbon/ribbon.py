"""
https://en.wikipedia.org/wiki/Ribbon_diagram

https://behreajj.medium.com/scripting-curves-in-blender-with-python-c487097efd13


Carson, M.; Bugg, C. E. (1986),
"Algorithm for Ribbon Models of Proteins",
Journal of Molecular Graphics, 4 (2): 121-122,
doi:10.1016/0263-7855(86)80010-8.

"""

import bpy
from time import time
import numpy as np
from batoms.ribbon.protein import Protein
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def draw_curve_from_vertices_bezier(name, data,
                                    coll, extrude=None,
                                    bevel_depth=None,
                                    bevel_control=None,
                                    node_type='Principled BSDF',
                                    node_inputs=None,
                                    material_style='metallic',
                                    backface_culling=True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 10
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    spline = crv.splines.new(type='BEZIER')
    nvert = len(vertices)
    spline.bezier_points.add(nvert-1)
    # vertices = np.append(vertices, np.zeros((nvert, 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    spline.bezier_points.foreach_set('co', vertices)
    for i in range(nvert):
        spline.bezier_points[i].handle_right_type = 'AUTO'
        spline.bezier_points[i].handle_left_type = 'AUTO'
    if extrude:
        crv.extrude = extrude
    if bevel_depth:
        crv.bevel_depth = bevel_depth
    # Create bevel control curve.
    if bevel_control:
        crv.bevel_object = bevel_control
        crv.bevel_mode = 'OBJECT'
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name,
                               data['color'],
                               node_type=node_type,
                               node_inputs=node_inputs,
                               material_style=material_style,
                               backface_culling=backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)


def draw_rope_from_vertices_nurbs(name, data, coll,
                                  extrude=None,
                                  bevel_depth=None,
                                  bevel_control=None,
                                  node_type='Principled BSDF',
                                  node_inputs=None,
                                  material_style='metallic',
                                  backface_culling=True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 10
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP'  # 'TANGENT' #'Z_UP'
    crv.twist_smooth = 1
    spline = crv.splines.new(type='NURBS')
    # spline.use_endpoint_u = True
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis=1)
    vertices = vertices.reshape(-1, 1)
    # tilts = data['tilts']
    spline.points.foreach_set('co', vertices)
    # spline.points.foreach_set('tilt', tilts)
    # spline.order_u = 4
    # Set the main curve's bevel control to the bevel control curve.
    if extrude:
        crv.extrude = extrude
    if bevel_depth:
        crv.bevel_depth = bevel_depth
    # Create bevel control curve.
    if bevel_control:
        crv.bevel_object = bevel_control
        crv.bevel_mode = 'OBJECT'
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name,
                               data['color'],
                               node_type=node_type,
                               node_inputs=node_inputs,
                               material_style=material_style,
                               backface_culling=backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)


def curve2mesh(name, vertices, resolution_u=20, coll=None):
    """
    Use Blender's curve to get b-spline
    resolution: number of vertices = (n(control points) - 1)*preview_U
    default preview_U  is 12
    """
    # tstart = time()
    crv = bpy.data.curves.new('%s-curve' % name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = resolution_u
    crv.fill_mode = 'FULL'
    spline = crv.splines.new(type='NURBS')
    # spline.use_endpoint_u = True
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis=1)
    vertices = vertices.reshape(-1, 1)
    spline.points.foreach_set('co', vertices)
    obj = bpy.data.objects.new(name, crv)
    # if coll:
    # coll.objects.link(obj)
    me = obj.to_mesh()
    n = len(me.vertices)
    vertices = np.zeros(3*n)
    # print('vertices: ', nvert, n)
    me.vertices.foreach_get('co', vertices)
    # print('curve2mesh: %s vertices, time: %s'%(nvert, (time() - tstart)))
    return vertices.reshape(-1, 3)


def draw_sheet_from_vertices_spline(name, data, coll,
                                    node_type='Principled BSDF',
                                    node_inputs=None,
                                    material_style='metallic',
                                    backface_culling=True,
                                    shade_smooth=True,):
    """
    """
    from batoms.material import create_material
    from batoms.ribbon.profile import build_mesh
    vertices = curve2mesh('%s-vertices' % name,
                          data['vertices'], data['resolution'], coll)
    sides = curve2mesh('%s-sides' % name,
                       data['sides'], data['resolution'])
    normals = curve2mesh('%s-normals' % name,
                         data['normals'], data['resolution'])
    profiles = data['profiles']
    scales = data['scales']
    vertices, faces = build_mesh(vertices, normals, sides, profiles, scales)
    me = bpy.data.meshes.new(name)
    me.from_pydata(vertices, [], faces)
    if shade_smooth:
        me.polygons.foreach_set('use_smooth', [True]*len(me.polygons))
    obj = bpy.data.objects.new(name, me)
    #
    material = create_material(name,
                               data['color'],
                               node_type=node_type,
                               node_inputs=node_inputs,
                               material_style=material_style,
                               backface_culling=backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)


class Ribbon():
    """

    """

    def __init__(self, label, batoms=None, datas={}, update=False) -> None:
        self.label = label
        self.batoms = batoms
        self.protein = Protein(batoms)
        self.protein.setting(datas)
        if update:
            self.protein.update()

    @property
    def coll(self):
        return self.get_coll()

    def get_coll(self):
        return bpy.data.collections.get('%s_ribbon' % self.label)

    def draw_sheet(self):
        tstart = time()
        for name, sheet in self.protein.sheets.items():
            draw_sheet_from_vertices_spline('sheet-%s' % name,
                                            sheet.as_dict(),
                                            self.coll,
                                            shade_smooth=False)
        logger.debug('draw sheet: %s' % (time() - tstart))
        self.batoms.selects['all'].show = False

    def draw_helix(self):
        tstart = time()
        for name, helix in self.protein.helixs.items():
            draw_sheet_from_vertices_spline('helix-%s' % name,
                                            helix.as_dict(),
                                            self.coll,
                                            shade_smooth=True)
        logger.debug('draw helix: %s' % (time() - tstart))
        self.batoms.selects['all'].show = False

    def draw_turn(self):
        # Create bevel control curve.
        tstart = time()
        bpy.ops.curve.primitive_bezier_circle_add(
            radius=0.25, enter_editmode=False)
        bevel_control = bpy.context.active_object
        bevel_control.hide_set(True)
        bevel_control.hide_render = True
        bevel_control.data.name = bevel_control.name = '%s_turn_bevel' % \
            self.label
        self.coll.objects.link(bevel_control)
        for name, turn in self.protein.turns.items():
            draw_rope_from_vertices_nurbs('turn-%s' % name,
                                          turn.as_dict(),
                                          self.coll,
                                          bevel_control=bevel_control)
        logger.debug('draw turn: %s' % (time() - tstart))
        self.batoms.selects['all'].show = False

    def draw(self):
        tstart = time()
        self.draw_sheet()
        self.draw_helix()
        self.draw_turn()
        logger.debug('draw ribbon: %s' % (time() - tstart))
