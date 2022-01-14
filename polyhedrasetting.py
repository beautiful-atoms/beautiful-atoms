"""
"""
from batoms.base import Setting, tuple2string
import numpy as np
from time import time
from batoms.butils import object_mode, clean_coll_objects
from batoms.bdraw import draw_surface_from_vertices, draw_cylinder
from batoms.bondsetting import build_bondlists


class PolyhedraSetting(Setting):
    """
    PolyhedraSetting object

    The PolyhedraSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    
    def __init__(self, label, batoms = None, polyhedrasetting = None) -> None:
        Setting.__init__(self, label, coll_name='%s_polyhedra'%label)
        self.name = 'bpolyhedra'
        self.batoms = batoms
        if len(self) == 0:
            self.set_default(self.batoms.species.species_props)
        if polyhedrasetting is not None:
            for key, data in polyhedrasetting.items():
                self[key] = data
    
    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.species = index
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
    
    def set_default(self, species_props):
        """
        """
        for sel, data in species_props.items():
            for sp, data in data.items():
                self[sp] = {
                    'flag': True,
                    'label': self.label,
                    'species': sp,
                    'color': np.append(data['color'][:3], 0.3),
                    'width': 0.005,
                }
    
    def add(self, polyhedras):
        if isinstance(polyhedras, str):
            polyhedras = [polyhedras]
        species = {sp: self.batoms.species_props[sp] for sp in polyhedras}
        self.set_default(species)
    
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Center                color         width \n'
        for p in self.collection:
            s += '{0:10s}   [{1:1.2f}  {2:1.2f}  {3:1.2f}  {4:1.2f}]   {5:1.3f} \n'.format(\
                p.species, p.color[0], p.color[1], p.color[2], p.color[3], p.width)
        s += '-'*60 + '\n'
        return s
    
    def draw_polyhedra(self, mask, bondlist = None):
        """Draw polyhedra.
        Parameters:
        bondlist: dict
        """
        object_mode()
        # clean_coll_objects(self.coll, 'polyhedra')
        frames = self.batoms.get_frames()
        arrays = self.batoms.arrays
        size = arrays['radius'][mask]*arrays['scale'][mask]
        species = arrays['species'][mask]
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        bond_datas = {}
        tstart = time()
        for f in range(nframe):
            print('update polyhedra: ', f)
            positions = frames[f][mask]
        if bondlist is None:
            bondlist = build_bondlists(species, positions, 
                        self.batoms.cell, self.batoms.pbc, self.batoms.bondsetting.cutoff_dict)
        polyhedra_kinds = build_polyhedralists(species, positions, self.batoms.cell, bondlist, 
                          self.batoms.bondsetting, self)
        for species, polyhedra_data in polyhedra_kinds.items():
            name = '%s_%s_polyhedra'%(self.label, species)
            draw_surface_from_vertices(name, 
                                datas = polyhedra_data,
                                use_smooth = False,
                                coll = self.coll,
                                )
            if polyhedra_data['show_edge']:
                name = '%s_%s_polyhedra_edge'%(self.label, species)
                draw_cylinder(name = name, 
                            datas = polyhedra_data['edge'],
                            coll = self.coll)


def build_polyhedralists(speciesarray, positions, cell,
                bondlists, bondsetting, polyhedrasetting):
        """
        """
        from scipy.spatial import ConvexHull
        from batoms.tools import get_polyhedra_kind
        tstart = time()
        polyhedra_kinds = {}
        polyhedra_dict = {}
        for b in bondsetting:
            spi = b.species1
            spj = b.species2
            if b.polyhedra:
                if spi not in polyhedra_dict: polyhedra_dict[spi] = []
                polyhedra_dict[spi].append(spj)
        # loop center atoms
        npositions = positions[bondlists[:, 1]] + np.dot(bondlists[:, 2:5], cell)
        for spi, spjs in polyhedra_dict.items():
            poly = polyhedrasetting[spi]
            polyhedra_kind = get_polyhedra_kind(color = poly.color, 
                                width = poly.width, show_edge = poly.show_edge)
            indis = np.where(speciesarray == spi)[0]
            for indi in indis:
                vertices = []
                indjs = np.where(bondlists[:, 0] == indi)[0]
                for indj in indjs:
                    if speciesarray[bondlists[indj, 1]] in spjs:
                        vertices.append(npositions[indj])
                nverts = len(vertices)
                if nverts >= 4:
                    # search convex polyhedra
                    hull = ConvexHull(vertices)
                    face = hull.simplices
                    nverts = len(polyhedra_kind['vertices'])
                    face = face + nverts
                    edge = []
                    for f in face:
                        edge.append([f[0], f[1]])
                        edge.append([f[0], f[2]])
                        edge.append([f[1], f[2]])
                    polyhedra_kind['vertices'] = polyhedra_kind['vertices'] + list(vertices)
                    polyhedra_kind['edges'] = polyhedra_kind['edges'] + list(edge)
                    polyhedra_kind['faces'] = polyhedra_kind['faces'] + list(face)
                    polyhedra_kind['battr_inputs'] = {'bpolyhedra': poly.as_dict()}
                    #------------------------------------------
                    # print('edge: ', edge)
                    for e in edge:
                        # print(e)
                        center = (polyhedra_kind['vertices'][e[0]] + polyhedra_kind['vertices'][e[1]])/2.0
                        vec = polyhedra_kind['vertices'][e[0]] - polyhedra_kind['vertices'][e[1]]
                        length = np.linalg.norm(vec)
                        nvec = vec/length
                        # print(center, nvec, length)
                        polyhedra_kind['edge']['lengths'].append(length)
                        polyhedra_kind['edge']['centers'].append(center)
                        polyhedra_kind['edge']['normals'].append(nvec)
                    polyhedra_kind['edge']['vertices'] = 6
                    polyhedra_kind['edge']['battr_inputs'] = {}
            if len(polyhedra_kind['vertices']) > 0:
                polyhedra_kinds[spi] = polyhedra_kind
        # print('build_polyhedralists: {0:10.2f} s'.format(time() - tstart))
        return polyhedra_kinds
