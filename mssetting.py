"""
Molecular surface:
1) van der Waals surface
2) Accessible surface area (SAS)
3) Solvent-excluded surface (SES) or Connolly surface
"""
import bpy
import numpy as np
from time import time
from batoms.base import Setting

default_colors = [(0.2, 0.1, 0.9, 1.0), (0.0, 0.0, 1.0, 1.0)]

class MSsetting(Setting):
    """
    MS object

    The MS object store the MS information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.
    probe: float
        probe radius, in [0,inf], default is 1.4
    """
    
    def __init__(self, label, probe = 1.4, batoms = None,
                render_resolution = 0.2,
                resolution = 0.4,
                threshold = 1e-4,
                update_method = 'FAST'
                ) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bMS'
        self.probe = probe
        self.batoms = batoms
        self.sas_name = '%s_sas'%self.label
        self.ses_name = '%s_ses'%self.label
        self.resolution = resolution
    
    def set_collection(self, label):
        """
        """
        if not bpy.data.collections.get(label):
            coll = bpy.data.collections.new(label)
            self.batoms.coll.children.link(coll)
            coll.batoms.flag = True
            coll.batoms.label = label
        
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'probe = %s \n'%(self.probe)
        s += '-'*60 + '\n'
        return s

    @property    
    def sas_obj(self):
        return bpy.data.objects.get(self.sas_name)
    
    @property    
    def sas_mb(self):
        return bpy.data.metaballs.get(self.sas_name)
    
    @property    
    def sas_mat(self):
        return bpy.data.materials.get(self.sas_name)
    
    @property    
    def ses_obj(self):
        return bpy.data.objects.get(self.ses_name)
    
    @property
    def spacing(self):
        return [self.resolution]*3
    
    def get_box(self, vertices, padding = 4):
        """
        Find a minmum box contains all atoms
        """
        from batoms.tools import get_box
        self.box = get_box(vertices, padding = padding)
        self.box_origin = self.box[:, 0]

    def build_grid(self):
        """
        generate gridpoints
        """
        from batoms.tools import build_grid
        self.meshgrids, self.shape = build_grid(self.box, self.resolution)
    
    def calc_isosurface(self, volume, level, spacing, origin = [0, 0, 0]):
        """
        Computes an isosurface from a volume grid.
        Parameters:
        volume: np.array
        cell: np.array
        level: float
        """
        from skimage import measure
        tstart = time()
        verts, faces, normals, values = measure.marching_cubes(volume, level = level,
                        spacing = spacing)
        verts += origin
        faces = list(faces)
        color = default_colors[0]
        isosurface = {'vertices': verts, 
                                'edges': [], 
                                'faces': faces, 
                                'color': color,
                                'battr_inputs': {'bisosurface': {}}
                            }
        print('Marching_cube: %s'%(time() - tstart))
        return isosurface
    
    def draw_SAS(self, resolution = 0.4, probe = 1.4, parallel = 1):
        """
        1) Generate meshgrids
        2) Calculate power distance for grids
        3) Marching cube find isosurface = 5
        """
        from batoms.bdraw import draw_surface_from_vertices, draw_vertices
        self.resolution = resolution
        self.probe = probe
        print('Resolution: {:1.3f}, Probe: {:1.3}'.format(self.resolution, self.probe))
        tstart = time()
        radii = np.array(self.batoms.radii_vdw) + self.probe
        self.get_box(self.batoms.positions, padding = max(radii) + self.resolution)
        self.build_grid()
        # draw_vertices('meshgrid', self.meshgrids)
        print('Grid Points: %s %s %s'%self.shape)
        positions = self.batoms.positions
        # ngrid = np.product(self.shape)
        # self.volume_sas = np.ones(ngrid)*4
        # self.indices_sas = self.query_radius(self.meshgrids, positions, radii)
        # self.volume_sas[self.indices_sas] = 6
        indices_sas, volume_sas = self.calc_power_distance(self.meshgrids, 
                                        positions, radii, parallel = parallel)
        volume = volume_sas.reshape(self.shape)
        isosurface = self.calc_isosurface(volume, 5, spacing = self.spacing, origin = self.box_origin)
        print('Vertices: %s'%len(isosurface['vertices']))
        obj = bpy.data.objects.get(self.sas_name)
        if obj is not None:
            bpy.data.objects.remove(obj, do_unlink = True)
        draw_surface_from_vertices(self.sas_name, 
                            datas = isosurface,
                            coll = self.batoms.coll.children['%s_surface'%self.batoms.coll_name],
                            backface_culling = False,
                        )
        print('Draw SAS: %s'%(time() - tstart))

    def draw_SES(self, resolution = 0.4, probe = 1.4, parallel = 1):
        """
        1) Get SAS
        2) Point probe sphere on vertices of SAS mesh
        3) Calculate power distance of grids
        4) value outside SAS set a value smaller than 5
        5) Marching cube find isosurface = 5
        """
        from batoms.bdraw import draw_surface_from_vertices, draw_vertices
        self.resolution = resolution
        self.probe = probe
        print('Resolution: {:1.3f}, Probe: {:1.3}'.format(self.resolution, self.probe))
        tstart = time()
        radii = np.array(self.batoms.radii_vdw) + self.probe
        self.get_box(self.batoms.positions, padding = max(radii) + self.resolution + self.probe)
        self.build_grid()
        # draw_vertices('meshgrid', self.meshgrids)
        print('Grid Points: %s %s %s'%self.shape)
        positions = self.batoms.positions
        # build SAS
        indices_sas, volume_sas = self.calc_power_distance(self.meshgrids, 
                                        positions, radii, parallel = parallel)
        volume = volume_sas.reshape(self.shape)
        spacing = [self.resolution]*3
        isosurface = self.calc_isosurface(volume, 5, self.spacing, origin = self.box_origin)
        # build SAS with probe = -self.resolution
        radii = np.array(self.batoms.radii_vdw) - self.resolution
        volume_sas1 = np.where(volume_sas > 5, 6, volume_sas)
        indices = np.where(volume_sas < 5)[0]
        volume_sas1[indices]
        # select grid points inside sas for query
        indices_sas1, tmp = self.calc_power_distance(self.meshgrids[indices], 
                                        positions, radii, parallel = parallel)
        volume_sas1[indices] = tmp
        #----------------------------
        tstart1 = time()
        vertices = isosurface['vertices']
        nvert = len(vertices)
        radii = np.ones(nvert)*self.probe
        # select grid points between sas and sas1
        self.volume_ses = np.where(volume_sas > 5, 4, volume_sas)
        self.volume_ses = np.where(volume_sas1 < 5, 6, self.volume_ses)
        indices = np.where((volume_sas < 5) & (volume_sas1 > 5))[0]
        indices_ses, volume_ses = self.calc_power_distance(self.meshgrids[indices], 
                                        vertices, radii, parallel = parallel)
        self.volume_ses[indices] = volume_ses
        volume = self.volume_ses.reshape(self.shape)
        isosurface = self.calc_isosurface(volume, 5, self.spacing, origin = self.box_origin)
        print('Vertices: %s'%len(isosurface['vertices']))
        obj = bpy.data.objects.get(self.ses_name)
        if obj is not None:
            bpy.data.objects.remove(obj, do_unlink = True)
        draw_surface_from_vertices(self.ses_name, 
                            datas = isosurface,
                            coll = self.batoms.coll.children['%s_surface'%self.batoms.coll_name],
                            backface_culling = False,
                        )
        print('Time SES: %s'%(time() - tstart))
        self.get_sesa()
    
    def get_sasa(self):
        """
        """
        from batoms.butils import get_area, get_volume
        me = self.sas_obj.data
        area = get_area(me)
        volume = get_volume(me)
        print('Area: {:5.3f},    Volume: {:5.3f}'.format(area, volume))
        return area, volume
    
    def get_psasa(self):
        """
        """
        me = self.sas_obj.data
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
            parea = np.sum(areas[atom_indices==j])
            pareas[j] = [symbols[j], parea]
            print('species: {}, area: {:1.3f}'.format(symbols[j], parea))
        # print(pareas_list)
        return pareas
    
    def get_sesa(self):
        """
        """
        from batoms.butils import get_area, get_volume
        me = self.ses_obj.data
        area = get_area(me)
        volume = get_volume(me)
        print('SES: area: {:1.3f}, volume: {:1.3f}'.format(area, volume))
        return area, volume
        
    def to_mesh(self):
        """
        
        """
        tstart = time()
        depsgraph = bpy.context.evaluated_depsgraph_get()
        me = self.sas_obj.evaluated_get(depsgraph).to_mesh()
        nvert = len(me.vertices)
        npoly = len(me.polygons)
        print('vertices: %s, polygons: %s.'%(nvert, npoly))
        print('SAS to mesh evaluated: %s'%(time() - tstart))
        return me
    
    def get_sasa_mb(self, frame_indices = []):
        """
        """
        from batoms.butils import get_area, get_volume
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
        print('Time for SAS area: %s'%(time() - tstart))
        return areas_list, pareas_list
    
    def get_psasa_mb(self, frame_indices = []):
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
            atom_indices, distance = self.map_face_to_atom(me, positions, radii)
            for j in range(natom):
                parea = np.sum(areas[atom_indices==j])
                pareas[j] = [symbols[j], parea]
                # print('{} atom: {}, area: {:1.3f}'.format(j, symbols[j], parea))
            pareas_list.append(pareas)
        # print(pareas_list)
        return pareas_list
    
    def map_face_to_atom(self, me, positions, radii, k = 1):
        """
        """
        tstart = time()
        npoly = len(me.polygons)
        centers = np.zeros(npoly*3)
        me.polygons.foreach_get('center', centers)
        centers = centers.reshape(-1, 3)
        indices, distance = self.calc_power_distance(centers, positions, radii, k=k)
        return indices, distance

    def map_vertice_to_atom(self, me, positions, radii, k = 1):
        """
        """
        tstart = time()
        nvertice = len(me.vertice)
        vertices = np.zeros(nvertice*3)
        me.vertices.foreach_get('co', vertices)
        vertices = vertices.reshape(-1, 3)
        indices, distance = self.calc_power_distance(vertices, positions, radii, k=k)
        return indices, distance

    def query_radius(self, points, positions, radii, k = 1):
        """
        Algorithm:
        Use KDTree to query the tree for neighbors within a radius r.
        """
        from scipy import spatial
        tstart = time()
        npoint = len(points)
        #----------------------------------------------------
        tstart = time()
        tree = spatial.KDTree(points)
        print('KDTree point: %s'%(time() - tstart))
        tstart = time()
        indices = tree.query_ball_point(positions, radii)
        print('KDTree query: %s'%(time() - tstart))
        tstart = time()
        indices = np.unique(np.concatenate(indices).astype(int))
        print('array indices: %s'%(time() - tstart))
        return indices

    def calc_power_distance(self, points, positions, radii, 
                            k = 1, parallel = 1):
        """
        Algorithm:
        Use KDTree to find the nearest atom for all vertices with
        the power distance as metric.
        """
        from scipy import spatial
        tstart = time()
        npoint = len(points)
        #----------------------------------------------------
        # change to 4-dimension to set power distance
        points = np.append(points, np.zeros((npoint, 1)), axis = 1)
        # maxr = max(max(radii), 5)
        maxr = 5
        radii = np.sqrt(maxr*maxr - np.square(radii)).reshape(-1, 1)
        positions = np.append(positions, radii, axis = 1)
        tstart = time()
        tree = spatial.KDTree(positions)
        print('KDTree positions: %s'%(time() - tstart))
        tstart = time()
        distance, indices = tree.query(points, k = k, workers=parallel)
        print('KDTree query: %s'%(time() - tstart))
        return indices, distance

    def to_mesh_object(self):
        """
        """
        import bmesh
        mbmesh = self.to_mesh()
        bpy.data.objects.remove(self.sas_obj, do_unlink = True)
        bpy.data.metaballs.remove(self.sas_mb, do_unlink = True)
        bm = bmesh.new()
        bm.from_mesh(mbmesh)
        tmp=bpy.data.meshes.new(name='Mesh')
        bm.to_mesh(tmp)
        obj = bpy.data.objects.new(self.sas_name, tmp)
        obj.data.materials.append(self.sas_mat)
        # print('SAS to mesh evaluated: %s'%(time() - tstart))
        coll = self.batoms.coll.children['%s_surface'%self.batoms.coll_name]
        coll.objects.link(obj)
        return obj
    
    def draw_SES_mb(self, resolution = 0.4,
                threshold = 1e-4,
                stiffness = 1,
                area = True,
                smooth = None,
                subdivide = 0,
                refine_normal = True,
                refine = True,
                steps = 5,
                ):
        """
        Time: 
        sas_mesh: 1.1
        shift_vertices: 0.09
        vertices_smooth: 0.8
        total: 2.55
        """
        name = '%s_ses'%self.label
        tstart0 = time()
        positions = self.batoms.positions
        radii = np.array(self.batoms.radii_vdw) + self.probe
        obj = self.build_SAS(resolution = resolution,
                threshold = threshold,
                stiffness = stiffness)
        obj = self.to_mesh_object()
        print('sas_mesh: %s'%(time() - tstart0))
        obj.name = name
        eps = resolution/2
        #-------------------------------------
        self.SES_subdivide(positions, radii, eps, subdivide)
        # return
        #----------------------
        me = self.ses_obj.data
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        normals = np.zeros(nvert*3)
        me.vertices.foreach_get('co', vertices)
        me.vertices.foreach_get('normal', normals)
        vertices = vertices.reshape(-1, 3)
        normals = normals.reshape(-1, 3)
        indices, mask1, mask2, mask3 = self.SES_map_to_atom(vertices, 
                                    positions, radii, eps)#, obj = self.ses_obj)
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
        # self.SES_vertices_smooth(mask, resolution = resolution, subdivide = subdivide, smooth = smooth)
        if refine:
            mask = ~mask1
            for i in range(steps):     
                self.SES_vertices_smooth(mask, resolution = resolution, subdivide = 0, smooth = 3)
                self.SES_refine_position(positions, indices, mask1, mask2, mask3)
        self.SES_vertices_smooth(mask, resolution = resolution, subdivide = subdivide, smooth = smooth)
        print('Time: %s'%(time() - tstart0))
        if area:
            self.get_sesa()
        return obj
    
    def SES_map_to_atom(self, vertices, positions, radii, eps, obj = None):
        """
        map vertices to atoms
        1) three sphere
        2) two sphere
        1) one sphere: contact surface
        """
        from batoms.tools import calc_origin_2, check_origin_3, calc_origin_3, check_origin_2
        from batoms.bdraw import draw_vertices
        print('eps: ', eps)
        n = len(vertices)
        indices, distance = self.calc_power_distance(vertices, positions, radii, k = 3)
        delta0 = abs(distance[:, 0] - distance[:, 1])
        delta1 = abs(distance[:, 1] - distance[:, 2])
        delta2 = abs(distance[:, 0] - distance[:, 2])
        # three spheres
        indices3 = np.where((delta0 < eps) & (delta1 < eps) & (delta2 < eps))[0]
        origins0 = positions[indices[:, 0][indices3]]
        origins1 = positions[indices[:, 1][indices3]]
        origins2 = positions[indices[:, 2][indices3]]
        vertices3 = vertices[indices3]
        indices31, origins_probe3 = check_origin_3(vertices3,    
                    origins0, origins1, origins2,
                    radii[indices[:, 0][indices3]],
                    radii[indices[:, 1][indices3]],
                    radii[indices[:, 2][indices3]],
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
        vertices2= vertices[indices2]
        indices21, origins_probe2 = check_origin_2(vertices2, 
                    origins0, origins1, 
                    radii[indices[:, 0][indices2]],
                    radii[indices[:, 1][indices2]],
                    0,
                    5*eps)
        indices2 = indices2[indices21]
        mask2 = np.zeros(n, dtype=bool)
        mask2[indices2] = True
        draw_vertices('origin2', origins_probe2)
        #---------------------------
        mask1 = np.where(delta0 > eps, True, False)
        print('indices: 1, 2, 3 ', len(indices), np.count_nonzero(mask1),
                     np.count_nonzero(mask2), np.count_nonzero(mask3))
        #----------------
        self.origins_probe2 = origins_probe2
        self.origins_probe3 = origins_probe3
        if obj is not None:
            self.batoms_mesh(obj, mask3, mode = 'VERT')
        return indices, mask1, mask2, mask3
    
    def select_mesh(self, obj, mask, mode = 'VERT'):
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
        if mode == 'VERT':
            obj.data.vertices.foreach_set('select', mask)
        elif mode == 'EDGE':
            obj.data.edges.foreach_set('select', mask)
        elif mode == 'FACE':
            obj.data.polygons.foreach_set('select', mask)
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.select_mode(type=mode)
    
    def SES_subdivide(self, positions, radii, eps, number_cuts = 2):
        """
        subdivide edge
        """
        if number_cuts < 1:
            return
        me = self.ses_obj.data
        nedge = len(me.edges)
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        iverts = np.zeros(nedge*2, dtype = 'int')
        me.vertices.foreach_get('co', vertices)
        me.edges.foreach_get('vertices', iverts)
        vertices = vertices.reshape(-1, 3)
        iverts = iverts.reshape(-1, 2)
        centers = (vertices[iverts[:, 0]] + vertices[iverts[:, 1]])/2
        indices, distance = self.calc_power_distance(centers, positions, radii, k = 3)
        delta0 = abs(distance[:, 0] - distance[:, 1])
        delta1 = abs(distance[:, 0] - distance[:, 2])
        delta0 = np.absolute(distance[:, 0] - distance[:, 1])
        normals1 = vertices[iverts[:, 0]] - vertices[iverts[:, 1]]
        lengths = np.linalg.norm(normals1, axis = 1)
        normals1 = normals1/lengths[:, None]
        normals2 = positions[indices[:, 0]] - positions[indices[:, 1]]
        normals2 = normals2/np.linalg.norm(normals2, axis = 1)[:, None]
        proj = np.absolute(np.einsum('ij,ij->i', normals1, normals2))
        # mask = np.where((((delta0 < eps) & (proj > 0.1)) | (delta1 < eps)) & (lengths > eps/2/number_cuts), True, False)
        mask = np.where((((delta0 < eps)) | (delta1 < eps)) & (lengths > eps/2/number_cuts), True, False)
        # print(eps, lengths[2608], mask[2608])
        self.batoms_mesh(self.ses_obj, mask, mode = 'EDGE')
        bpy.ops.mesh.subdivide(number_cuts = number_cuts, smoothness = 0, fractal_along_normal = 0)
        # return
        bpy.ops.object.mode_set(mode = 'OBJECT')

    def SES_refine_normal(self, vertices, normals, positions, 
                    indices, indices1, indices2, indices3):
        #-----------------------------------------
        origins0 = positions[indices[:, 0][indices1]]
        vertices0 = vertices[indices1]
        vecs0 = vertices0 - origins0
        vecs0 = vecs0/np.linalg.norm(vecs0, axis = 1)[:, None]
        normals[indices1] = vecs0
        #-------------------
        origins0 = positions[indices[:, 0][indices2]]
        origins1 = positions[indices[:, 1][indices2]]
        vertices1 = vertices[indices2]
        xaxis = origins1 - origins0
        l2 = np.linalg.norm(xaxis, axis = 1)
        # Heron's formula
        xaxis = xaxis/l2[:, None]
        vecs1 = vertices1 - origins0
        zaxis = np.cross(xaxis, vecs1)
        zaxis = zaxis/np.linalg.norm(zaxis, axis = 1)[:, None]
        projz = np.einsum('ij,ij->i', normals[indices2], zaxis)[:, None]*zaxis
        normals[indices2] = normals[indices2] - projz
        return normals
    
    def SES_vertices_smooth(self, mask, resolution = 0.4, subdivide = 0, smooth = None):
        if smooth is None:
            repeat = int(1/resolution)
        else:
            repeat = smooth
        if repeat == 0:
            return
        repeat += subdivide
        tstart = time()
        bpy.context.view_layer.objects.active = self.ses_obj
        self.batoms_mesh(self.ses_obj, mask, mode = 'VERT')
        if repeat > 0:
            bpy.ops.mesh.vertices_smooth(factor = 0.5, repeat = repeat)
        bpy.ops.object.mode_set(mode = 'OBJECT')
        print('smooth_vertices: %s'%(time() - tstart))

    def SES_refine_position(self, positions, indices, indices1, indices2, indices3):
        from batoms.bdraw import draw_vertices
        radii_vdw = np.array(self.batoms.radii_vdw)
        obj = self.ses_obj
        me = obj.data
        #-------------------------------------------------------
        tstart = time()
        nvert = len(me.vertices)
        vertices = np.zeros(nvert*3)
        # print('vertices: %s'%(nvert))
        me.vertices.foreach_get('co', vertices)
        vertices = vertices.reshape(-1, 3)
        #-----------------------------------------
        # print('indices1: ', len(indices1))
        if True:
            origins0 = positions[indices[:, 0][indices1]]
            vertices1 = vertices[indices1]
            vecs0 = vertices1 - origins0
            dist = np.linalg.norm(vecs0, axis = 1)[:, None]
            radii_vdw0 = radii_vdw[indices[:, 0][indices1]]
            normals0 = vecs0/dist
            vertices[indices1] = origins0 + radii_vdw0[:, None]*normals0
        #----------
        # print('indices2: ', len(indices2))
        vertices2= vertices[indices2]
        vecs2 = vertices2 - self.origins_probe2
        dist = np.linalg.norm(vecs2, axis = 1)[:, None]
        normals2 = vecs2/dist
        vertices[indices2] = self.origins_probe2 + self.probe*normals2
        # draw_vertices('origin2', origins_probe)
        #----------------------------------------------------------
        vertices3 = vertices[indices3]
        vecs3 = vertices3 - self.origins_probe3
        dist = np.linalg.norm(vecs3, axis = 1)[:, None]
        normals3 = vecs3/dist
        vertices[indices3] = self.origins_probe3 + self.probe*normals3
        # draw_vertices('origin3', origins_probe)
        #----------------------------------------------------------
        vertices = vertices.reshape(-1, 1)
        me.vertices.foreach_set('co', vertices)
        me.update()
        print('refine vertices: %s'%(time() - tstart))

    def get_sesa_mb(self):
        """
        """
        from batoms.butils import get_area, get_volume
        me = self.ses_obj.data
        area = get_area(me)
        volume = get_volume(me)
        print('SES: area: {:1.3f}, volume: {:1.3f}'.format(area, volume))
        return area, volume

    def build_SAS_Shrake_Rupley(self, probe = 1.4, subdivisions = 2):
        """
        Algorithm: ShrakeRupley
        """
        from scipy import spatial
        #
        tstart = time()
        self.probe = probe
        atoms = self.batoms
        pbc = atoms.pbc
        cell = atoms.cell
        positions = atoms.positions
        radii = np.array(self.batoms.radii_vdw) + self.probe
        # nli, nlj = primitive_neighbor_list('ij', pbc, cell, positions, radii, self_interaction=False)
        # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
        # bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
        # nb = len(nli)
        na = len(positions)
        #
        tstart = time()
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = subdivisions)
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
        scales = np.concatenate((scales, scales, scales), axis = 1)
        spheres = spheres*scales[:, None]
        #
        origini = positions
        spheres += origini[:, None]
        print('build_spheres: {0:10.2f} s'.format(time() - tstart))
        #
        spheres = spheres.reshape(-1, 3)
        indices, distances = self.calc_power_distance(spheres, positions, radii)
        mask = np.where(distances > 4.98)
        spheres = spheres[mask]
        #
        tstart1 = time()
        faces = []
        '''
        # tri = spatial.Delaunay(spheres)
        # faces = tri.simplices
        tree = spatial.KDTree(spheres)
        print('KDTree positions: %s'%(time() - tstart))
        distance, indices = tree.query(spheres, k = 3)
        indices0 = np.arange(len(spheres))
        indices0 = indices0.reshape(-1, 1)
        faces = np.concatenate((indices0, indices[:, 1:]), axis = 1)
        faces = faces.tolist()
        print('build_faces: %s'%(time() - tstart1))
        '''
        #
        color = default_colors[0]
        SAS_Shrake_Rupley = {'vertices': spheres.reshape(-1, 3), 
                            'edges': [], 
                            # 'faces': [],
                            'faces': faces,
                            # 'faces': list(tri.simplices),
                            'color': default_colors[0],
                            'battr_inputs': {},
                        }
        print('vertices: ', len(SAS_Shrake_Rupley['vertices']))
        print('build_SAS: %s'%(time() - tstart))
        self.SAS_Shrake_Rupley = SAS_Shrake_Rupley
        bpy.data.objects.remove(source, do_unlink=True)
        # return 0
    
    def draw_SAS_Shrake_Rupley(self):
        from batoms.bdraw import draw_surface_from_vertices
        name = '%s_%s'%(self.label, 'SAS')
        draw_surface_from_vertices(name, 
                                datas = self.SAS_Shrake_Rupley,
                                coll = self.batoms.coll.children['%s_volume'%self.label],
                            )
'''

    def build_SAS_Voronoi(self, positions, radii, cell, subdivisions = 3):
        """
        Algorithm: voronoi
        """
        natom = len(positions)
        datas = calc_voronoi(positions, radii, cell)
        tstart = time()
        SAS = {}
        verts = []
        edges = []
        faces = []
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = subdivisions)
        # bpy.ops.mesh.primitive_cylinder_add()
        source = bpy.context.view_layer.objects.active
        me = source.data
        n = len(me.vertices)
        sphere0 = np.empty(n*3, dtype=np.float64)
        me.vertices.foreach_get('co', sphere0)  
        sphere0 = sphere0.reshape((n, 3))
        for i in range(natom):
            sphere = sphere0.copy()
            sphere *=radii[i]
            sphere += positions[i]
            indices = []
            for face in datas[i]['faces']:
                if face['adjacent_cell'] < 0: continue
                points = [datas[i]['vertices'][j] for j in face['vertices'][0:3]]
                v1 = points[1] - points[0]
                v2 = points[2] - points[0]
                normal = np.cross(v1, v2)
                index = sameside(normal, points[0], positions[i], sphere)
                indices.extend(index)
            # print(indices)
            verts.extend(np.delete(sphere, indices, axis = 0))
        # print(verts)
        #
        color = default_colors[0]
        SAS = {'vertices': verts, 
                            'edges': [], 
                            'faces': [],
                            # 'faces': list(tri.simplices),
                            'color': default_colors[0],
                            'battr_inputs': {},
                        }
        print('build_SAS: %s'%(time() - tstart))
        bpy.data.objects.remove(source, do_unlink=True)
        return SAS
    
    def build_voronoi(self, positions, radii, cell):
        """
        """
        tstart = time()
        datas = calc_voronoi(positions, radii, cell)
        faces = []
        vertices = []
        for data in datas:
            nvert = len(vertices)
            vertices.extend(data['vertices'])
            faces1 = [[x + nvert for x in face['vertices']] for face in data['faces']]
            faces.extend(faces1)
        voronoi = {}
        probe = self.probe
        color = default_colors[0]
        voronoi = {'vertices': vertices, 
                            'edges': [], 
                            'faces': faces, 
                            'color': color,
                            'battr_inputs': {},
                        }
        print('calc_voronoi: %s'%(time() - tstart))
        return voronoi

def calc_voronoi(positions, radii, cell, block_size = 4):
    """
    Computes an radical voronoi cell from a points.
    
    Parameters:

    volume: np.array

    cell: np.array

    level: float
    
    """
    from pyvoro import compute_voronoi
    tstart = time()
    datas = compute_voronoi(
                positions,
                cell,
                block_size, # block size
                radii,
    )
    print('calc_voronoi: %s'%(time() - tstart))
    return datas

def sameside(normal, point, center, sphere):
    """
    find points at the same side of origin for all planes
    Todo: '+' or '-'
    """
    d = np.dot(point, normal)
    x1 = np.dot(center, normal) - d
    x2 = np.dot(sphere, normal) - d
    x = x1*x2
    index = np.where(x<-1e-6)[0]
    # print(index)
    return index



'''
