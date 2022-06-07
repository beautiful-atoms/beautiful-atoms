import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       )
import bmesh
from batoms import Batoms
from ase import Atoms
from ase.build.rotate import rotation_matrix_from_points
from ase.build import rotate
from ase.data import covalent_radii, chemical_symbols
from ase.geometry import get_distances
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


molecules = {
    'C': Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
               positions=[[0, 0, 0],
                          [0, 0, -1.09],
                          [0, 1.024, 0.36],
                          [-0.89, -0.51, 0.36],
                          [0.89, -0.51, 0.36]]),
    'N': Atoms(symbols=["N", "H", "H", "H"],
               positions=[[0, 0, 0],
                          [0, 0, 1.02],
                          [0.975, 0, 0.289],
                          [-0.385, 0.896, 0.289]]),
    'O': Atoms(symbols=["O", "H", "H"],
               positions=[[0, 0, 0],
                          [0, 0, -0.969],
                          [0.94, 0, 0.234]]),
    'S': Atoms(symbols=["S", "H", "H"],
               positions=[[0, 0, 0],
                          [0, 0.001, -1.364],
                          [1.364, 0.003, 0]]),
    'F': Atoms(symbols=["F"],
               positions=[[0, 0, 0]]),
    'Cl': Atoms(symbols=["Cl"],
                positions=[[0, 0, 0]]),
    'Br': Atoms(symbols=["Br"],
                positions=[[0, 0, 0]]),
}


def edit_atom(batoms, indices, element):
    """Replace by element

    Args:
        element (str): The target element
    """
    logger.debug('replace atoms')
    positions = batoms.positions
    species = batoms.arrays['species']
    bondlists = batoms.bonds.bondlists
    removeH = []
    for i in indices:
        # copy mol
        mol = molecules[element].copy()
        # number of H atoms to be added
        nh = len(mol) - 1
        # how many atoms (nbond) are connected to atoms i
        bai0 = np.where((bondlists[:, 0] == i))[0]
        bai1 = np.where((bondlists[:, 1] == i))[0]
        bai = np.append(bondlists[bai0, 1], bondlists[bai1, 0]).astype(int)
        nbond = len(bai)
        sps = species[bai]
        indh = bai[np.where(sps == 'H')[0]]
        # difference between nbond and nh
        dbond = nbond - nh
        if dbond >= 0:
            # nbond >= nh,
            # 1. no more new H will be add
            # 2. dh old H will be remove
            dh = min(len(indh), dbond)
            removeH.extend(indh[:dh].tolist())
            mol = mol[0:1]
            mol[0].position = positions[i]
        else:
            # nbond < nh,
            # 1. dh new H will be add
            # 2. no more old H will be remove
            dh = min(nh, -dbond)
            if nbond == 1:
                # second nearest neighbour need be search
                j = bai[0]
                a1 = mol[0].position - mol[1].position
                a2 = positions[i] - positions[j]
                bondlength = covalent_radii[chemical_symbols.index(element)] + \
                    covalent_radii[chemical_symbols.index(species[j])]
                mol.positions += positions[j] + \
                    a2/np.linalg.norm(a2)*bondlength
                # how many atoms (nbond) are connected to atoms i
                baj0 = np.where((bondlists[:, 0] == j))[0]
                baj1 = np.where((bondlists[:, 1] == j))[0]
                baj = np.append(bondlists[baj0, 1],
                                bondlists[baj1, 0]).astype(int)
                baj = baj[np.where(baj != i)[0]]
                b2 = np.cross([0, 0, 1], a2)
                rotate(mol, a1, a2, [1, 0, 0], b2, center=mol[0].position)
                # find another verctor
                maxd = 0
                maxi = 0
                for i in range(20):
                    angle = i*9
                    mol1 = mol.copy()
                    mol1.rotate(angle, a2, center=mol[0].position)
                    p1 = mol1.positions[2:]
                    p2 = positions[baj]
                    D, D_len = get_distances(p1, p2)
                    maxd1 = np.min(D_len)
                    if maxd1 > maxd:
                        maxd = maxd1
                        maxi = i
                mol.rotate(maxi*9, a2, center=mol[0].position)
                del mol[[1]]
                # find minimized overlap
            else:
                baj = np.insert(bai, 0, i, axis=0)
                # minimize_rotation_and_translation from ASE
                p = mol.positions[:(nbond + 1), :]
                p0 = positions[baj]
                # centeroids to origin
                c = np.mean(p, axis=0)
                p -= c
                c0 = np.mean(p0, axis=0)
                p0 -= c0
                # Compute rotation matrix
                R = rotation_matrix_from_points(p.T, p0.T)
                mol.set_positions(np.dot(mol.positions, R.T) + c0)
                del mol[range(1, (nbond + 1))]
        # rotate mol
        batoms.add_atoms({'species': mol.get_chemical_symbols(),
                          'positions': mol.positions})
    indices.extend(removeH)
    batoms.delete(indices)
    batoms.model_style = 1


class MolecueEditElement(Operator):
    bl_idname = "batoms.molecule_edit_atom"
    bl_label = "Edit Atoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Edit Atoms")

    element: StringProperty(
        name="Element", default='C',
        description="element")
    bond_order: IntProperty(
        name="Bond order", default=1,
        description="bond order")

    @classmethod
    def poll(cls, context):
        return context.mode in {'EDIT_MESH'} and \
            context.object.batoms.type != 'OTHER'

    def execute(self, context):
        obj = context.object
        data = obj.data
        if data.total_vert_sel > 0:
            bm = bmesh.from_edit_mesh(data)
            indices = [s.index for s in bm.select_history
                       if isinstance(s, bmesh.types.BMVert)]
            self.report({'INFO'}, '%s atoms were replaced' % len(indices))
            batoms = Batoms(label=obj.batoms.label)
            edit_atom(batoms, indices, self.element)
            bpy.context.view_layer.objects.active = batoms.obj
            bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}
