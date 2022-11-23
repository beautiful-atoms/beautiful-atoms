"""Definition of the plane class.

This module defines the plane object in the Batoms package.

"""

import bpy
import bmesh
from time import time
import numpy as np
from batoms.base.object import BaseObject
from .setting import MolecularSurfaceSettings
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class MolecularSurface(BaseObject):
    def __init__(self,
                 label=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 ):
        """Plane Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = 'plane'
        BaseObject.__init__(self, label, name)
        self.settings = MolecularSurfaceSettings(
            self.label, parent=self)
        self.settings.bpy_data.active = True

    def build_materials(self, name, color, node_inputs=None,
                        material_style='default',
                        vertex_color = None,
                        color_by_attribute = None):
        """
        """
        from batoms.material import create_material
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                              color=color,
                              node_inputs=node_inputs,
                              material_style=material_style,
                              backface_culling=False,
                              vertex_color=vertex_color,
                              color_by_attribute=color_by_attribute)
        return mat

    def draw(self, ms_name="ALL"):
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old plane
        clean_coll_object_by_type(self.batoms.coll, 'MS')
        for ms in self.settings.bpy_setting:
            if ms_name.upper() != "ALL" and ms.name != ms_name:
                continue
            if ms.type == 'SAS':
                self.draw_SAS(ms)
            elif ms.type == 'SES':
                self.draw_SES(ms)

    @property
    def sas_objs(self):
        sas_objs = {}
        for ms in self.settings.bpy_setting:
            sas_name = '%s_%s_sas' % (self.label, ms.name)
            obj = bpy.data.objects.get(sas_name)
            sas_objs[ms.name] = obj
        return sas_objs

    @property
    def sas_mb(self):
        return bpy.data.metaballs.get(self.sas_name)

    @property
    def sas_mat(self):
        return bpy.data.materials.get(self.sas_name)

    @property
    def ses_objs(self):
        ses_objs = {}
        for ms in self.settings.bpy_setting:
            sas_name = '%s_%s_ses' % (self.label, ms.name)
            obj = bpy.data.objects.get(sas_name)
            ses_objs[ms.name] = obj
        return ses_objs

    def get_space(self, resolution):
        return [resolution]*3

    def get_box(self, vertices, padding=4):
        """
        Find a minmum box contains all atoms
        """
        from batoms.utils import get_box
        self.box = get_box(vertices, padding=padding)
        self.box_origin = self.box[:, 0]
        return self.box, self.box_origin

    def build_grid(self, resolution):
        """
        generate gridpoints
        """
        from batoms.utils import build_grid
        self.meshgrids, self.shape = build_grid(self.box, resolution)

    def calc_isosurface(self, volume, level, spacing, origin=[0, 0, 0],
                        color=(0.2, 0.1, 0.9, 1.0)):
        """
        Computes an isosurface from a volume grid.
        Parameters:
        volume: np.array
        cell: np.array
        level: float
        """
        from skimage import measure
        tstart = time()
        verts, faces, normals, values = \
            measure.marching_cubes(volume,
                                   level=level,
                                   spacing=spacing)
        verts += origin
        faces = list(faces)
        isosurface = {'vertices': verts,
                      'edges': [],
                      'faces': faces,
                      'color': color,
                      'battr_inputs': {'Bmolecularsurface': {}}
                      }
        logger.debug('Marching_cube: %s' % (time() - tstart))
        return isosurface

    def draw_SAS(self, ms, parallel=1):
        """
        1) Generate meshgrids
        2) Calculate power distance for grids
        3) Marching cube find isosurface = 5
        """
        from batoms.draw import draw_surface_from_vertices
        resolution = ms.resolution
        probe = ms.probe
        logger.debug('Resolution: {:1.3f}, Probe: {:1.3}'.format(
            resolution, probe))
        tstart = time()
        indices = self.batoms.selects[ms.select].indices
        if len(indices) == 0:
            return
        radii = np.array(self.batoms.radii_vdw[indices]) + probe
        positions = self.batoms.local_positions[indices]
        self.get_box(positions, padding=max(radii) + resolution)
        self.build_grid(resolution=resolution)
        logger.debug('Grid Points: %s %s %s' % self.shape)
        indices_sas, volume_sas = \
            self.calc_power_distance(self.meshgrids,
                                     positions,
                                     radii,
                                     parallel=parallel)
        volume = volume_sas.reshape(self.shape)
        isosurface = self.calc_isosurface(
            volume, 5, spacing=self.get_space(resolution),
            origin=self.box_origin)
        isosurface['color'] = ms.color
        isosurface['material_style'] = ms.material_style
        logger.debug('Vertices: %s' % len(isosurface['vertices']))
        sas_name = '%s_%s_sas' % (self.label, ms.name)
        self.delete_obj(sas_name)
        coll = self.batoms.coll.children['%s_surface' %
                                         self.batoms.coll_name]
        obj = draw_surface_from_vertices(sas_name,
                                         datas=isosurface,
                                         coll=coll,
                                         )
        #
        if ms.color_by != "None":
            from ase.cell import Cell
            from batoms.utils import map_volumetric_data
            from batoms.utils.attribute import set_mesh_attribute
            if ms.color_by.upper() == "ELECTROSTATIC_POTENTIAL":
                if 'charges' not in self.batoms._attributes:
                    self.batoms.auto_assign_charge()
                data = self.batoms.calc_electrostatic_potential(isosurface['vertices'])
            else:
                # using volumetric_data
                # scaled positions
                scaled_verts = Cell(self.batoms.cell).scaled_positions(isosurface['vertices'])
                data = map_volumetric_data(
                                self.batoms.volumetric_data[ms.color_by],
                                scaled_verts
                                )
            # normalize
            data = (data - np.min(data))/(np.max(data) - np.min(data))
            # color_attribute = map_color(data, [1, 0, 0, ms.transparency],
                                    # [0, 0, 1, ms.transparency])
            # set_vertex_color(obj, ms.color_by, color_attribute)
            obj.data.attributes.new(name='{}_data'.format(ms.color_by),
                                type='FLOAT', domain='POINT')
            set_mesh_attribute(obj, '{}_data'.format(ms.color_by), data)
            color_by_attribute = {'attribute_name': '{}_data'.format(ms.color_by),
                              'ValToRGB':[ms.color1[:],
                                        ms.color2[:]]
                                        }
        else:
            color_by_attribute = None
        mat = self.build_materials(sas_name, color=isosurface['color'],
                                   material_style=isosurface['material_style'],
                                #    vertex_color = ms.color_by,
                                color_by_attribute = color_by_attribute,
                                   )
        obj.data.materials.append(mat)

        obj.parent = self.batoms.obj
        obj.batoms.type = 'MS'
        obj.batoms.label = self.batoms.label
        logger.debug('Draw SAS: %s' % (time() - tstart))

    def draw_SES(self, ms, parallel=1):
        """
        1) Get SAS
        2) Point probe sphere on vertices of SAS mesh
        3) Calculate power distance of grids
        4) value outside SAS set a value smaller than 5
        5) Marching cube find isosurface = 5
        """
        from batoms.draw import draw_surface_from_vertices
        resolution = ms.resolution
        probe = ms.probe
        logger.debug('Resolution: {:1.3f}, Probe: {:1.3}'.format(
            resolution, probe))
        tstart = time()
        indices = self.batoms.selects[ms.select].indices
        if len(indices) == 0:
            return
        radii_vdw = np.array(self.batoms.radii_vdw[indices])
        radii = radii_vdw + probe
        positions = self.batoms.local_positions[indices]
        self.get_box(positions, padding=max(radii) + resolution + probe)
        self.build_grid(resolution=resolution)
        # draw_vertices('meshgrid', self.meshgrids)
        logger.debug('Grid Points: %s %s %s' % self.shape)
        # build SAS
        indices_sas, volume_sas = \
            self.calc_power_distance(self.meshgrids,
                                     positions,
                                     radii,
                                     parallel=parallel)
        volume = volume_sas.reshape(self.shape)
        isosurface = self.calc_isosurface(
            volume, 5, self.get_space(resolution), origin=self.box_origin)
        # build SAS with probe = -self.resolution
        radii = radii_vdw - resolution
        volume_sas1 = np.where(volume_sas > 5, 6, volume_sas)
        indices = np.where(volume_sas < 5)[0]
        volume_sas1[indices]
        # select grid points inside sas for query
        indices_sas1, tmp = \
            self.calc_power_distance(self.meshgrids[indices],
                                     positions, radii,
                                     parallel=parallel)
        volume_sas1[indices] = tmp
        # ----------------------------
        # tstart1 = time()
        vertices = isosurface['vertices']
        nvert = len(vertices)
        radii = np.ones(nvert)*probe
        # select grid points between sas and sas1
        self.volume_ses = np.where(volume_sas > 5, 4, volume_sas)
        self.volume_ses = np.where(volume_sas1 < 5, 6, self.volume_ses)
        indices = np.where((volume_sas < 5) & (volume_sas1 > 5))[0]
        indices_ses, volume_ses = \
            self.calc_power_distance(self.meshgrids[indices],
                                     vertices, radii,
                                     parallel=parallel)
        self.volume_ses[indices] = volume_ses
        volume = self.volume_ses.reshape(self.shape)
        isosurface = self.calc_isosurface(
            volume, 5, self.get_space(resolution), origin=self.box_origin)
        isosurface['color'] = ms.color
        isosurface['material_style'] = ms.material_style
        logger.debug('Vertices: %s' % len(isosurface['vertices']))
        ses_name = '%s_%s_ses' % (self.label, ms.name)
        self.delete_obj(ses_name)
        coll = self.batoms.coll.children['%s_surface' %
                                         self.batoms.coll_name]
        obj = draw_surface_from_vertices(ses_name,
                                         datas=isosurface,
                                         coll=coll,
                                         )
        mat = self.build_materials(ses_name, color=isosurface['color'],
                                   material_style=isosurface['material_style'],
                                   )
        obj.data.materials.append(mat)
        obj.parent = self.batoms.obj
        obj.batoms.type = 'MS'
        obj.batoms.label = self.batoms.label
        logger.debug('Time SES: %s' % (time() - tstart))
        self.get_sesa(ms.name)

    def get_sasa(self, name):
        """
        """
        from batoms.utils.butils import get_area, get_volume
        me = self.sas_objs[name].data
        area = get_area(me)
        volume = get_volume(me)
        print('Area: {:5.3f},    Volume: {:5.3f}'.format(area, volume))
        return area, volume

    def get_psasa(self):
        """
        """
        me = self.sas_objs.data
        npoly = len(me.polygons)
        areas = np.zeros(npoly)
        me.polygons.foreach_get('area', areas)
        pareas = {}
        arrays = self.batoms.arrays
        positions = arrays['positions']
        symbols = arrays['species']
        natom = len(positions)
        radii = np.array(self.batoms.radii_vdw) + self.probe
        atom_indices, distance = self.map_face_to_atom(me, positions, radii)
        for j in range(natom):
            parea = np.sum(areas[atom_indices == j])
            pareas[j] = [symbols[j], parea]
            print('species: {}, area: {:1.3f}'.format(symbols[j], parea))
        # print(pareas_list)
        return pareas

    def get_sesa(self, name):
        """
        """
        from batoms.utils.butils import get_area, get_volume
        me = self.ses_objs[name].data
        area = get_area(me)
        volume = get_volume(me)
        print('SES: area: {:1.3f}, volume: {:1.3f}'.format(area, volume))
        return area, volume

    def to_mesh(self):
        """

        """
        # tstart = time()
        depsgraph = bpy.context.evaluated_depsgraph_get()
        me = self.sas_objs.evaluated_get(depsgraph).to_mesh()
        nvert = len(me.vertices)
        npoly = len(me.polygons)
        logger.debug('vertices: %s, polygons: %s.' % (nvert, npoly))
        # print('SAS to mesh evaluated: %s' % (time() - tstart))
        return me

    def get_sasa_mb(self, frame_indices=[]):
        """
        """
        from batoms.utils.butils import get_area, get_volume
        if isinstance(frame_indices, int):
            frame_indices = [frame_indices]
        if frame_indices == []:
            frame_indices = range(self.nframe)
        areas_list = []
        pareas_list = []
        tstart = time()
        print('Frame    SAS Area    SAS Volume')
        for i in frame_indices:
            bpy.context.scene.frame_set(i)
            me = self.to_mesh()
            area = get_area(me)
            volume = get_volume(me)
            print('{:4d}  {:11.3f}  {:11.3f}'.format(i, area, volume))
            areas_list.append(area)
        logger.debug('Time for SAS area: %s' % (time() - tstart))
        return areas_list, pareas_list

    def get_psasa_mb(self, frame_indices=[]):
        """
        """
        if isinstance(frame_indices, int):
            frame_indices = [frame_indices]
        if frame_indices == []:
            frame_indices = range(self.nframe)
        pareas_list = []
        radii = np.array(self.batoms.radii_vdw) + self.probe
        for i in frame_indices:
            bpy.context.scene.frame_set(i)
            me = self.to_mesh()
            npoly = len(me.polygons)
            areas = np.zeros(npoly)
            me.polygons.foreach_get('area', areas)
            pareas = {}
            atoms = self.batoms.atoms
            positions = atoms.positions
            symbols = atoms.get_chemical_symbols()
            natom = len(positions)
            atom_indices, distance = self.map_face_to_atom(
                me, positions, radii)
            for j in range(natom):
                parea = np.sum(areas[atom_indices == j])
                pareas[j] = [symbols[j], parea]
            pareas_list.append(pareas)
        # print(pareas_list)
        return pareas_list

    def map_face_to_atom(self, me, positions, radii, k=1):
        """
        """
        npoly = len(me.polygons)
        centers = np.zeros(npoly*3)
        me.polygons.foreach_get('center', centers)
        centers = centers.reshape(-1, 3)
        indices, distance = self.calc_power_distance(
            centers, positions, radii, k=k)
        return indices, distance

    def map_vertice_to_atom(self, me, positions, radii, k=1):
        """
        """
        nvertice = len(me.vertice)
        vertices = np.zeros(nvertice*3)
        me.vertices.foreach_get('co', vertices)
        vertices = vertices.reshape(-1, 3)
        indices, distance = self.calc_power_distance(
            vertices, positions, radii, k=k)
        return indices, distance

    def query_radius(self, points, positions, radii, k=1):
        """
        Algorithm:
        Use KDTree to query the tree for neighbors within a radius r.
        """
        from scipy import spatial
        tstart = time()
        # ----------------------------------------------------
        tstart = time()
        tree = spatial.KDTree(points)
        logger.debug('KDTree point: %s' % (time() - tstart))
        tstart = time()
        indices = tree.query_ball_point(positions, radii)
        logger.debug('KDTree query: %s' % (time() - tstart))
        tstart = time()
        indices = np.unique(np.concatenate(indices).astype(int))
        logger.debug('array indices: %s' % (time() - tstart))
        return indices

    def calc_power_distance(self, points, positions, radii,
                            k=1, parallel=1):
        """
        Algorithm:
        Use KDTree to find the nearest atom for all vertices with
        the power distance as metric.
        """
        from scipy import spatial
        tstart = time()
        npoint = len(points)
        # ----------------------------------------------------
        # change to 4-dimension to set power distance
        points = np.append(points, np.zeros((npoint, 1)), axis=1)
        # maxr = max(max(radii), 5)
        maxr = 5
        radii = np.sqrt(maxr*maxr - np.square(radii)).reshape(-1, 1)
        positions = np.append(positions, radii, axis=1)
        tstart = time()
        tree = spatial.KDTree(positions)
        logger.debug('KDTree positions: %s' % (time() - tstart))
        tstart = time()
        distance, indices = tree.query(points, k=k, workers=parallel)
        logger.debug('KDTree query: %s' % (time() - tstart))
        return indices, distance

    def to_mesh_object(self):
        """
        """
        import bmesh
        mbmesh = self.to_mesh()
        bpy.data.objects.remove(self.sas_objs, do_unlink=True)
        bpy.data.metaballs.remove(self.sas_mb, do_unlink=True)
        bm = bmesh.new()
        bm.from_mesh(mbmesh)
        tmp = bpy.data.meshes.new(name='Mesh')
        bm.to_mesh(tmp)
        obj = bpy.data.objects.new(self.sas_name, tmp)
        obj.data.materials.append(self.sas_mat)
        # print('SAS to mesh evaluated: %s'%(time() - tstart))
        coll = self.batoms.coll.children['%s_surface' % self.batoms.coll_name]
        coll.objects.link(obj)
        return obj

    def draw_SES_mb(self, resolution=0.4,
                    threshold=1e-4,
                    stiffness=1,
                    area=True,
                    smooth=None,
                    subdivide=0,
                    refine_normal=True,
                    refine=True,
                    steps=5,
                    ):
        """
        Time:
        sas_mesh: 1.1
        shift_vertices: 0.09
        vertices_smooth: 0.8
        total: 2.55
        """
        name = '%s_ses' % self.label
        tstart0 = time()
        positions = self.batoms.positions
        radii = np.array(self.batoms.radii_vdw) + self.probe
        obj = self.build_SAS(resolution=resolution,
                             threshold=threshold,
                             stiffness=stiffness)
        obj = self.to_mesh_object()
        logger.debug('sas_mesh: %s' % (time() - tstart0))
        obj.name = name
        eps = resolution/2
        # -------------------------------------
        self.SES_subdivide(positions, radii, eps, subdivide)
        # return
        # ----------------------
        me = self.ses_obj.data
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        normals = np.zeros(nvert*3)
        me.vertices.foreach_get('co', vertices)
        me.vertices.foreach_get('normal', normals)
        vertices = vertices.reshape(-1, 3)
        normals = normals.reshape(-1, 3)
        indices, mask1, mask2, mask3 = self.SES_map_to_atom(vertices,
                                                            positions,
                                                            radii,
                                                            eps)
        # return
        # refine_normal = False
        if refine_normal:
            normals = self.SES_refine_normal(vertices, normals, positions,
                                             indices, mask1, mask2, mask3)
        vertices = vertices - normals*self.probe
        vertices = vertices.reshape(-1, 1)
        me.vertices.foreach_set('co', vertices)
        # return
        # ----------------------------------
        mask = np.ones(nvert, dtype=bool)
        if refine:
            mask = ~mask1
            for i in range(steps):
                self.SES_vertices_smooth(
                    mask, resolution=resolution, subdivide=0, smooth=3)
                self.SES_refine_position(
                    positions, indices, mask1, mask2, mask3)
        self.SES_vertices_smooth(
            mask, resolution=resolution, subdivide=subdivide, smooth=smooth)
        logger.debug('Time: %s' % (time() - tstart0))
        if area:
            self.get_sesa()
        return obj

    def SES_map_to_atom(self, vertices, positions, radii, eps, obj=None):
        """
        map vertices to atoms
        1) three sphere
        2) two sphere
        1) one sphere: contact surface
        """
        from batoms.utils import check_origin_3, check_origin_2
        from batoms.draw import draw_vertices
        logger.debug('eps: ', eps)
        n = len(vertices)
        indices, distance = self.calc_power_distance(
            vertices, positions, radii, k=3)
        delta0 = abs(distance[:, 0] - distance[:, 1])
        delta1 = abs(distance[:, 1] - distance[:, 2])
        delta2 = abs(distance[:, 0] - distance[:, 2])
        # three spheres
        indices3 = np.where((delta0 < eps) & (
            delta1 < eps) & (delta2 < eps))[0]
        origins0 = positions[indices[:, 0][indices3]]
        origins1 = positions[indices[:, 1][indices3]]
        origins2 = positions[indices[:, 2][indices3]]
        vertices3 = vertices[indices3]
        indices31, origins_probe3 = check_origin_3(vertices3,
                                                   origins0,
                                                   origins1,
                                                   origins2,
                                                   radii[indices[:, 0]
                                                         [indices3]],
                                                   radii[indices[:, 1]
                                                         [indices3]],
                                                   radii[indices[:, 2]
                                                         [indices3]],
                                                   0,
                                                   5*eps)
        indices3 = indices3[indices31]
        mask3 = np.zeros(n, dtype=bool)
        mask3[indices3] = True
        draw_vertices('origin3', origins_probe3)
        # two spheres
        indices2 = np.where((delta0 < eps))[0]
        origins0 = positions[indices[:, 0][indices2]]
        origins1 = positions[indices[:, 1][indices2]]
        vertices2 = vertices[indices2]
        indices21, origins_probe2 = check_origin_2(vertices2,
                                                   origins0, origins1,
                                                   radii[indices[:, 0]
                                                         [indices2]],
                                                   radii[indices[:, 1]
                                                         [indices2]],
                                                   0,
                                                   5*eps)
        indices2 = indices2[indices21]
        mask2 = np.zeros(n, dtype=bool)
        mask2[indices2] = True
        draw_vertices('origin2', origins_probe2)
        # ---------------------------
        mask1 = np.where(delta0 > eps, True, False)
        logger.debug('indices: 1, 2, 3 ', len(indices), np.count_nonzero(mask1),
              np.count_nonzero(mask2), np.count_nonzero(mask3))
        # ----------------
        self.origins_probe2 = origins_probe2
        self.origins_probe3 = origins_probe3
        if obj is not None:
            self.batoms_mesh(obj, mask3, mode='VERT')
        return indices, mask1, mask2, mask3

    def select_mesh(self, obj, mask, mode='VERT'):
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        if mode == 'VERT':
            obj.data.vertices.foreach_set('select', mask)
        elif mode == 'EDGE':
            obj.data.edges.foreach_set('select', mask)
        elif mode == 'FACE':
            obj.data.polygons.foreach_set('select', mask)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_mode(type=mode)

    def SES_subdivide(self, positions, radii, eps, number_cuts=2):
        """
        subdivide edge
        """
        if number_cuts < 1:
            return
        me = self.ses_obj.data
        nedge = len(me.edges)
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        iverts = np.zeros(nedge*2, dtype='int')
        me.vertices.foreach_get('co', vertices)
        me.edges.foreach_get('vertices', iverts)
        vertices = vertices.reshape(-1, 3)
        iverts = iverts.reshape(-1, 2)
        centers = (vertices[iverts[:, 0]] + vertices[iverts[:, 1]])/2
        indices, distance = self.calc_power_distance(
            centers, positions, radii, k=3)
        delta0 = abs(distance[:, 0] - distance[:, 1])
        delta1 = abs(distance[:, 0] - distance[:, 2])
        delta0 = np.absolute(distance[:, 0] - distance[:, 1])
        normals1 = vertices[iverts[:, 0]] - vertices[iverts[:, 1]]
        lengths = np.linalg.norm(normals1, axis=1)
        normals1 = normals1/lengths[:, None]
        normals2 = positions[indices[:, 0]] - positions[indices[:, 1]]
        normals2 = normals2/np.linalg.norm(normals2, axis=1)[:, None]
        # proj = np.absolute(np.einsum('ij,ij->i', normals1, normals2))
        mask = np.where((((delta0 < eps)) | (delta1 < eps)) &
                        (lengths > eps/2/number_cuts), True, False)
        # print(eps, lengths[2608], mask[2608])
        self.batoms_mesh(self.ses_obj, mask, mode='EDGE')
        bpy.ops.mesh.subdivide(number_cuts=int(number_cuts),
                               smoothness=0, fractal_along_normal=0)
        # return
        bpy.ops.object.mode_set(mode='OBJECT')

    def SES_refine_normal(self, vertices, normals, positions,
                          indices, indices1, indices2, indices3):
        # -----------------------------------------
        origins0 = positions[indices[:, 0][indices1]]
        vertices0 = vertices[indices1]
        vecs0 = vertices0 - origins0
        vecs0 = vecs0/np.linalg.norm(vecs0, axis=1)[:, None]
        normals[indices1] = vecs0
        # -------------------
        origins0 = positions[indices[:, 0][indices2]]
        origins1 = positions[indices[:, 1][indices2]]
        vertices1 = vertices[indices2]
        xaxis = origins1 - origins0
        l2 = np.linalg.norm(xaxis, axis=1)
        # Heron's formula
        xaxis = xaxis/l2[:, None]
        vecs1 = vertices1 - origins0
        zaxis = np.cross(xaxis, vecs1)
        zaxis = zaxis/np.linalg.norm(zaxis, axis=1)[:, None]
        projz = np.einsum('ij,ij->i', normals[indices2], zaxis)[:, None]*zaxis
        normals[indices2] = normals[indices2] - projz
        return normals

    def SES_vertices_smooth(self, mask, resolution=0.4,
                            subdivide=0, smooth=None):
        if smooth is None:
            repeat = int(1/resolution)
        else:
            repeat = smooth
        if repeat == 0:
            return
        repeat += subdivide
        tstart = time()
        bpy.context.view_layer.objects.active = self.ses_obj
        self.batoms_mesh(self.ses_obj, mask, mode='VERT')
        if repeat > 0:
            bpy.ops.mesh.vertices_smooth(factor=0.5, repeat=repeat)
        bpy.ops.object.mode_set(mode='OBJECT')
        logger.debug('smooth_vertices: %s' % (time() - tstart))

    def SES_refine_position(self, positions, indices, indices1,
                            indices2, indices3):
        radii_vdw = np.array(self.batoms.radii_vdw)
        obj = self.ses_obj
        me = obj.data
        # -------------------------------------------------------
        tstart = time()
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        # print('vertices: %s'%(nvert))
        me.vertices.foreach_get('co', vertices)
        vertices = vertices.reshape(-1, 3)
        # -----------------------------------------
        # print('indices1: ', len(indices1))
        if True:
            origins0 = positions[indices[:, 0][indices1]]
            vertices1 = vertices[indices1]
            vecs0 = vertices1 - origins0
            dist = np.linalg.norm(vecs0, axis=1)[:, None]
            radii_vdw0 = radii_vdw[indices[:, 0][indices1]]
            normals0 = vecs0/dist
            vertices[indices1] = origins0 + radii_vdw0[:, None]*normals0
        # ----------
        # print('indices2: ', len(indices2))
        vertices2 = vertices[indices2]
        vecs2 = vertices2 - self.origins_probe2
        dist = np.linalg.norm(vecs2, axis=1)[:, None]
        normals2 = vecs2/dist
        vertices[indices2] = self.origins_probe2 + self.probe*normals2
        # draw_vertices('origin2', origins_probe)
        # ----------------------------------------------------------
        vertices3 = vertices[indices3]
        vecs3 = vertices3 - self.origins_probe3
        dist = np.linalg.norm(vecs3, axis=1)[:, None]
        normals3 = vecs3/dist
        vertices[indices3] = self.origins_probe3 + self.probe*normals3
        # draw_vertices('origin3', origins_probe)
        # ----------------------------------------------------------
        vertices = vertices.reshape(-1, 1)
        me.vertices.foreach_set('co', vertices)
        me.update()
        logger.debug('refine vertices: %s' % (time() - tstart))

    def get_sesa_mb(self):
        """
        """
        from batoms.utils.butils import get_area, get_volume
        me = self.ses_obj.data
        area = get_area(me)
        volume = get_volume(me)
        print('SES: area: {:1.3f}, volume: {:1.3f}'.format(area, volume))
        return area, volume

    def build_SAS_Shrake_Rupley(self, probe=1.4, subdivisions=2):
        """
        Algorithm: ShrakeRupley
        """
        default_colors = [(0.2, 0.1, 0.9, 1.0), (0.0, 0.0, 1.0, 1.0)]
        #
        tstart = time()
        self.probe = probe
        atoms = self.batoms
        positions = atoms.positions
        radii = np.array(self.batoms.radii_vdw) + self.probe
        na = len(positions)
        #
        tstart = time()
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=subdivisions)
        # bpy.ops.mesh.primitive_cylinder_add()
        source = bpy.context.view_layer.objects.active
        me = source.data
        n = len(me.vertices)
        sphere0 = np.empty(n*3, dtype=np.float64)
        me.vertices.foreach_get('co', sphere0)
        sphere0 = sphere0.reshape((n, 3))
        spheres = np.tile(sphere0, (na, 1, 1))
        #
        radii = np.array(radii).reshape(-1, 1)
        scales = radii
        scales = np.concatenate((scales, scales, scales), axis=1)
        spheres = spheres*scales[:, None]
        #
        origini = positions
        spheres += origini[:, None]
        logger.debug('build_spheres: {0:10.2f} s'.format(time() - tstart))
        #
        spheres = spheres.reshape(-1, 3)
        indices, distances = self.calc_power_distance(
            spheres, positions, radii)
        mask = np.where(distances > 4.98)
        spheres = spheres[mask]
        #
        # tstart1 = time()
        faces = []
        '''
        # tri = spatial.Delaunay(spheres)
        # faces = tri.simplices
        tree = spatial.KDTree(spheres)
        logger.debug('KDTree positions: %s'%(time() - tstart))
        distance, indices = tree.query(spheres, k = 3)
        indices0 = np.arange(len(spheres))
        indices0 = indices0.reshape(-1, 1)
        faces = np.concatenate((indices0, indices[:, 1:]), axis = 1)
        faces = faces.tolist()
        logger.debug('build_faces: %s'%(time() - tstart1))
        '''
        #
        SAS_Shrake_Rupley = {'vertices': spheres.reshape(-1, 3),
                             'edges': [],
                             # 'faces': [],
                             'faces': faces,
                             # 'faces': list(tri.simplices),
                             'color': default_colors[0],
                             'battr_inputs': {},
                             }
        logger.debug('vertices: ', len(SAS_Shrake_Rupley['vertices']))
        logger.debug('build_SAS: %s' % (time() - tstart))
        self.SAS_Shrake_Rupley = SAS_Shrake_Rupley
        bpy.data.objects.remove(source, do_unlink=True)
        # return 0

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data
