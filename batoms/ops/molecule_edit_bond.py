import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       BoolProperty,
                       )
import bmesh
from batoms import Batoms
from ase.data import covalent_radii, chemical_symbols
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def edit_bond(batoms, indices, order):
    """Change bond order

    Args:
        order (int): The target order
    """
    logger.debug('edit bond')
    positions = batoms.positions
    species = batoms.arrays['species']
    bondlists = batoms.bonds.bondlists
    removeH = []
    for i in indices:
        # bond order
        order0 = batoms.bonds[i].order
        # bond atoms idnex
        ai = bondlists[i, 0]
        aj = bondlists[i, 1]
        #
        # how many atoms (nbond) are connected to atoms i and j
        bai0 = np.where((bondlists[:, 0] == ai))[0]
        bai1 = np.where((bondlists[:, 1] == ai))[0]
        baj0 = np.where((bondlists[:, 0] == aj))[0]
        baj1 = np.where((bondlists[:, 1] == aj))[0]
        bai = np.append(bondlists[bai0, 1], bondlists[bai1, 0]).astype(int)
        baj = np.append(bondlists[baj0, 1], bondlists[baj1, 0]).astype(int)
        nbondi = len(bai)
        nbondj = len(baj)
        spsi = species[bai]
        spsj = species[baj]
        indhi = bai[np.where(spsi == 'H')[0]]
        indhj = bai[np.where(spsj == 'H')[0]]
        # difference between order and target order
        dorder = order - order0
        if dorder >= 0:
            # nbond >= nh,
            # 1. no more new H will be add
            # 2. dh old H will be remove
            dhi = min(len(indhi), dorder)
            dhj = min(len(indhj), dorder)
            removeH.extend(indhi[:dhi].tolist())
            removeH.extend(indhj[:dhj].tolist())
        else:
            # nbond < nh,
            # 1. dh new H will be add
            # 2. no more old H will be remove
            dhi = dorder
            dhj = dorder

    indices.extend(removeH)
    batoms.delete(indices)
    batoms.model_style = 1


class MolecueEditBond(Operator):
    bl_idname = "batoms.molecule_edit_bond"
    bl_label = "Edit bond"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Edit bond")

    order: IntProperty(
        name="order", default=1,
        description="order")
    # TODO rotate around bond
    rotate: IntProperty(
        name="rotate", default=1,
        description="rotate")
    hydrogen: BoolProperty(
        name="Ajust hydrogen", default=1,
        description="Ajust hydrogen")

    @classmethod
    def poll(cls, context):
        return context.mode in {'EDIT_MESH'}

    def execute(self, context):
        obj = context.object
        if not obj.batoms.bbond.flag:
            logger.critical('Please select a Batom.')
            return {'FINISHED'}
        data = obj.data
        if data.total_vert_sel > 0:
            bm = bmesh.from_edit_mesh(data)
            indices = [s.index for s in bm.select_history
                       if isinstance(s, bmesh.types.BMVert)]
            self.report({'INFO'}, '%s atoms were replaced' % len(indices))
            batoms = Batoms(label=obj.batoms.label)
            edit_bond(batoms, indices, self.order)
            bpy.context.view_layer.objects.active = batoms.obj
            bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}
