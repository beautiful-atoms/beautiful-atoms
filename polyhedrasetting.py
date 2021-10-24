"""
"""
from batoms.bondsetting import Setting
import numpy as np
from time import time

class PolyhedraSetting(Setting):
    """
    PolyhedraSetting object

    The PolyhedraSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    def __init__(self, label, polyhedrasetting = None) -> None:
        Setting.__init__(self, label)
        self.name = 'bpolyhedra'
        if len(self) == 0:
            self.set_default(self.species)
        if polyhedrasetting is not None:
            for key, data in polyhedrasetting.items():
                self[key] = data
    def set_default(self, species):
        """
        """
        for sp, data in species.items():
            self[sp] = {
                'species': sp,
                'color': np.append(data['color'][:3], 0.3),
                'edgewidth': 0.005,
            }
    def add(self, polyhedras):
        if isinstance(polyhedras, str):
            polyhedras = [polyhedras]
        species = {sp: self.species[sp] for sp in polyhedras}
        self.set_default(species)
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Center                color         edgewidth \n'
        for p in self.collection:
            s += '{0:10s}   [{1:1.2f}  {2:1.2f}  {3:1.2f}  {4:1.2f}]   {5:1.3f} \n'.format(\
                p.species, p.color[0], p.color[1], p.color[2], p.color[3], p.edgewidth)
        s += '-'*60 + '\n'
        return s


def build_polyhedralists(atoms, bondlists, bondsetting, polyhedrasetting):
        """
        """
        from scipy.spatial import ConvexHull
        from batoms.tools import get_polyhedra_kind
        tstart = time()
        if 'species' not in atoms.arrays:
            atoms.new_array('species', np.array(atoms.get_chemical_symbols()))
        speciesarray = np.array(atoms.arrays['species'])
        positions = atoms.positions
        polyhedra_kinds = {}
        polyhedra_dict = {}
        for b in bondsetting:
            spi = b.species1
            spj = b.species2
            if b.polyhedra:
                if spi not in polyhedra_dict: polyhedra_dict[spi] = []
                polyhedra_dict[spi].append(spj)
        # loop center atoms
        npositions = positions[bondlists[:, 1]] + np.dot(bondlists[:, 2:5], atoms.cell)
        for spi, spjs in polyhedra_dict.items():
            polyhedra_kind = get_polyhedra_kind(color = polyhedrasetting[spi].color, 
                                edgewidth = polyhedrasetting[spi].edgewidth)
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
                    # print('edge: ', edge)
                    for e in edge:
                        # print(e)
                        center = (polyhedra_kind['vertices'][e[0]] + polyhedra_kind['vertices'][e[1]])/2.0
                        vec = polyhedra_kind['vertices'][e[0]] - polyhedra_kind['vertices'][e[1]]
                        length = np.linalg.norm(vec)
                        nvec = vec/length
                        # print(center, nvec, length)
                        polyhedra_kind['edge_cylinder']['lengths'].append(length)
                        polyhedra_kind['edge_cylinder']['centers'].append(center)
                        polyhedra_kind['edge_cylinder']['normals'].append(nvec)
            if len(polyhedra_kind['vertices']) > 0:
                polyhedra_kinds[spi] = polyhedra_kind
        # print('build_polyhedralists: {0:10.2f} s'.format(time() - tstart))
        return polyhedra_kinds
