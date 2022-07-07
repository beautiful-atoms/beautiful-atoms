import bpy
from batoms.ribbon.profile import ellipse, rectangle
from batoms.base.collection import Setting
from time import time
import numpy as np
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class Protein():
    def __init__(self, batoms=None,
                 chains={}, residues=[],
                 sheets={},
                 helixs={},
                 turns={},
                 ) -> None:
        self.batoms = batoms
        self.sheetsetting = SheetSetting(self.batoms.label, batoms)
        self.helixsetting = HelixSetting(self.batoms.label, batoms)
        self.chains = chains
        self.residues = residues
        self.sheets = sheets
        self.helixs = helixs
        self.turns = turns

    def setting(self, datas):
        if datas is None:
            return
        if 'sheet' in datas:
            self.sheetsetting.from_dict(datas['sheet'])
        if 'helix' in datas:
            self.helixsetting.from_dict(datas['helix'])

    def update(self):
        tstart = time()
        arrays = self.batoms.arrays
        if 'types' not in arrays:
            return
        # build Chains
        chains = {}
        chainIDs = np.unique(arrays['chainids'])
        for chainID in chainIDs:
            chain = Chain(chainID)
            indices = np.where(arrays['chainids'] == chainID)[0]
            chain.indices = indices
            chains[chainID] = chain
        # build residues
        residues = []
        # get Sheet
        names = np.core.defchararray.add(
            arrays['chainids'], arrays['residuenumbers'].astype('U20'))
        residuenames, indices = np.unique(names, return_index=True)
        residuenames = residuenames[np.argsort(indices)]
        indices = indices[np.argsort(indices)]
        lastResi = np.where(names == residuenames[-1])[0][-1]
        indices = np.append(indices, lastResi + 1)
        # find real residue and its Ca, C, O
        nresi = len(residuenames)
        for i in range(nresi):
            resi = Residue(names[indices[i]],
                           arrays['residuenames'][indices[i]],
                           arrays['chainids'][indices[i]],
                           arrays['residuenumbers'][indices[i]],
                           )
            resi.indices = list(range(indices[i], indices[i + 1]))
            if 'ATOM' in arrays['types']:
                for j in range(indices[i], indices[i + 1]):
                    if arrays['atomtypes'][j] == 'CA':
                        resi.Ca = j
                    elif arrays['atomtypes'][j] == 'C':
                        resi.C = j
                    elif arrays['atomtypes'][j] == 'O':
                        resi.O = j
            else:
                resi.type = -1
            residues.append(resi)
        # calculate the plane
        for i in range(nresi - 1):
            if residues[i].type == -1 or residues[i + 1].type == -1:
                continue
            plane = GetPeptidePlane(
                residues[i], residues[i + 1], arrays['positions'])
            residues[i].plane = plane
        # init sheets
        sheets = {}
        for sheet in self.sheetsetting.bpy_setting:
            sheet1 = Sheet(sheet.name, 1, sheet.startChain,
                           sheet.endChain,
                           sheet.startResi,
                           sheet.endResi)
            sheet1.color = sheet.color
            sheets[sheet.name] = sheet1
        # init helixs
        helixs = {}
        for helix in self.helixsetting.bpy_setting:
            helix1 = Helix(helix.name, 2, helix.startChain,
                           helix.endChain,
                           helix.startResi,
                           helix.endResi)
            helix1.color = helix.color
            helixs[helix.name] = helix1
        for i in range(nresi - 1):
            if residues[i].type == -1 or residues[i + 1].type == -1:
                continue
            # add sheet
            for sheet in self.sheetsetting:
                if sheet.startChain == residues[i].chainID:
                    if residues[i].resSeq == sheet.startResi and i > 0:
                        residues[i].type = 1
                        sheets[sheet.name].append(residues[i - 1])
                        sheets[sheet.name].append(residues[i])
                        break
                    elif residues[i].resSeq > sheet.startResi and \
                            residues[i].resSeq < sheet.endResi:
                        residues[i].type = 1
                        sheets[sheet.name].append(residues[i])
                        break
                    elif residues[i].resSeq == sheet.endResi and i < nresi - 1:
                        residues[i].type = 1
                        sheets[sheet.name].append(residues[i])
                        sheets[sheet.name].append(residues[i + 1])
                        break
            # add helix
            for helix in self.helixsetting:
                if helix.startChain == residues[i].chainID:
                    if residues[i].resSeq == helix.startResi and i > 0:
                        residues[i].type = 2
                        helixs[helix.name].append(residues[i - 1])
                        helixs[helix.name].append(residues[i])
                        break
                    elif residues[i].resSeq > helix.startResi and \
                            residues[i].resSeq < helix.endResi:
                        residues[i].type = 2
                        helixs[helix.name].append(residues[i])
                        break
                    elif residues[i].resSeq == helix.endResi and i < nresi - 1:
                        residues[i].type = 2
                        helixs[helix.name].append(residues[i])
                        helixs[helix.name].append(residues[i + 1])
                        break
        # init turns
        turns = {}
        if residues[0].type == 0:
            turn = Turn(residues[0].name, 0, residues[0].chainID,
                        residues[0].resSeq,
                        )
            turns[turn.name] = turn
            turns[turn.name].append(residues[0])
        # add turn
        for i in range(1, nresi - 1):
            if residues[i-1].chainID != turn.startChain or  \
                    (residues[i-1].type != 0 and residues[i].type == 0):
                turn = Turn(residues[i].name, 0, residues[i].chainID,
                            residues[i].resSeq,
                            )
                turns[turn.name] = turn
                if i > 1:
                    turns[turn.name].append(residues[i - 2])
                turns[turn.name].append(residues[i - 1])
                turns[turn.name].append(residues[i])
            elif residues[i].type == 0 and residues[i + 1].type == 0:
                turns[turn.name].append(residues[i])
            elif residues[i].type == 0 and residues[i + 1].type != 0:
                turns[turn.name].append(residues[i])
                if i < nresi - 1:
                    turns[turn.name].append(residues[i + 1])
                    if i < nresi - 2:
                        turns[turn.name].append(residues[i + 2])

        self.residues = residues
        self.chains = chains
        self.sheets = sheets
        self.helixs = helixs
        self.turns = turns
        self.batoms.selects['all'].show = False
        # HETATM
        indices = np.where(arrays['types'] == 'HETATM')[0]
        self.batoms.selects.add('heta', indices)
        logger.debug('update ribbon: %s' % (time() - tstart))


class Chain():
    def __init__(self, chainID, residues=[]) -> None:
        self.chainID = chainID
        self.residues = residues
        self.indices = []


class Residue():
    """
    Type:
        0:
        1: Sheet
        2: Helix
        -1: water or ...
    """

    def __init__(self, name, resName, chainID, resSeq, type=0) -> None:
        self.name = name
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.Ca = 0
        self.C = 0
        self.O = 0
        self.plane = None
        self.type = type
        self.indices = []


class secondaryStructure():
    """
    Type:
        0:
        1: Sheet
        2: Helix
        -1: water or ...
    """

    def __init__(self, name, type,
                 startChain, endChain,
                 startResi, endResi,
                 ) -> None:
        self.name = name
        self.type = type
        self.startChain = startChain
        self.endChain = endChain
        self.startResi = startResi
        self.endResi = endResi
        self.residues = []
        self.indices = []
        self.Ca = []
        self.C = []
        self.O = []
        self.extrude = 1.0
        self.depth = 0.25
        self.resolution = 20

    def append(self, residue):
        self.residues.append(residue)
        self.indices.extend(residue.indices)
        self.Ca.append(residue.Ca)
        self.C.append(residue.C)
        self.O.append(residue.O)

    @property
    def positions(self):
        n = len(self.residues)
        positions = np.ones((n, 3))
        for i in range(n):
            positions[i] = self.residues[i].plane['position']
        return positions

    @property
    def profiles(self):
        profiles = rectangle(self.extrude, self.depth)
        return profiles

    @property
    def sides(self):
        n = len(self.residues)
        sides = np.ones((n, 3))
        for i in range(n):
            sides[i] = self.residues[i].plane['side']
        return sides

    @property
    def normals(self):
        n = len(self.residues)
        normals = np.ones((n, 3))
        for i in range(n):
            normals[i] = self.residues[i].plane['normal']
        return normals

    @property
    def tilts(self):
        n = len(self.residues)
        tilts = np.zeros(n)
        for i in range(n):
            tilts[i] = np.arccos(self.residues[i].plane['side'][2])
        return tilts

    @property
    def scales(self):
        n = len(self.residues)
        scales = np.ones((self.resolution*(n-1), 3))
        scales[-self.resolution:, 0] *= np.linspace(1.5, 0.0, self.resolution)
        return scales

    def as_dict(self):
        """
         A complication arises when the direction of the carbonyl
        oxygen flips, as is always the case between adjacent
        residues of B sheets.
        """
        n = len(self.residues)
        for i in range(1, n):
            if np.dot(self.residues[i - 1].plane['side'],
                      self.residues[i].plane['side']) < 0:
                self.residues[i].plane['flipped'] = True
                self.residues[i].plane['normal'] = - \
                    self.residues[i].plane['normal']
                self.residues[i].plane['side'] = - \
                    self.residues[i].plane['side']
        data = {
            'vertices': self.positions,
            'color': self.color,
            'tilts': self.tilts,
            'sides': self.sides,
            'normals': self.normals,
            'profiles': self.profiles,
            'scales': self.scales,
            'depth': self.depth,
            'resolution': self.resolution,
        }
        return data


class Sheet(secondaryStructure):
    def __init__(self, name, type,
                 startChain, endChain,
                 startResi, endResi,
                 ) -> None:
        secondaryStructure.__init__(self, name, type,
                                    startChain, endChain,
                                    startResi, endResi,)
        self.type = 1


class Helix(secondaryStructure):
    """
    """

    def __init__(self, name, type,
                 startChain, endChain,
                 startResi, endResi,
                 ) -> None:
        secondaryStructure.__init__(self, name, type,
                                    startChain, endChain,
                                    startResi, endResi,)
        self.type = 2

    @property
    def positions(self):
        n = len(self.residues)
        positions = np.ones((n, 3))
        for i in range(n):
            positions[i] = self.residues[i].plane['position'] + \
                self.residues[i].plane['normal']*1.5
        return positions

    @property
    def profiles(self):
        profiles = ellipse(32, self.extrude, self.depth)
        return profiles

    @property
    def scales(self):
        n = len(self.residues)
        scales = np.ones((self.resolution*(n-1), 3))
        return scales


class Turn(secondaryStructure):
    """
    """

    def __init__(self, name, type,
                 startChain=None, endChain=None,
                 startResi=None, endResi=None,
                 ) -> None:
        secondaryStructure.__init__(self, name, type,
                                    startChain, endChain,
                                    startResi, endResi,)
        self.type = 2
        self.extrude = 1.0
        self.depth = 0.3
        self.color = [0, 1, 0, 1]

    @property
    def positions(self):
        n = len(self.residues)
        positions = np.ones((n, 3))
        for i in range(n):
            if self.residues[i].type == 2:
                positions[i] = self.residues[i].plane['position'] + \
                    self.residues[i].plane['normal']*1.5
            else:
                positions[i] = self.residues[i].plane['position']
        return positions

    def as_dict(self):
        data = {'vertices': self.positions,
                'color': self.color,
                'tilts': self.tilts,
                'extrudes': self.extrude,
                'radius': self.depth}
        return data


class SheetSetting(Setting):
    """
    sheet = {
    'name': '%s-%s'%(chainId, sheetId),
    'startChain': startChain,
    'startResi': startResi,
    'endChain': endChain,
    'endResi': endResi,
    }
    """

    def __init__(self, label, batoms=None,
                 ) -> None:
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Bprotein'
        self.batoms = batoms
        self.sheet_datas = {}

    def get_bpy_setting(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            data = getattr(coll, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.settings_sheet

    @property
    def show(self):
        return self.get_show()

    @show.setter
    def show(self, state):
        self.set_show(state)

    def get_show(self):
        return self.batoms.attributes['show'][self.indices]

    def set_show(self, show, only_atoms=True):
        show0 = self.batoms.show
        show0[self.indices] = show
        self.batoms.set_attributes({'show': show0})


class HelixSetting(Setting):
    """

    """

    def __init__(self, label, batoms=None,
                 ) -> None:
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.name = 'Bprotein'
        self.batoms = batoms
        self.helix_datas = {}

    def get_bpy_setting(self):
        if self.coll_name:
            coll = bpy.data.collections.get(self.coll_name)
            data = getattr(coll, self.name)
        else:
            raise KeyError("The collection property {} not exist!".format(self.name))
        return data.settings_helix

def GetPeptidePlane(resi1, resi2, positions):
    """
    Peptide plane: The atoms of the group, O=C-N-H,
    are fixed on the same plane, known as the peptide plane.
    The whole plane may rotate around the N-Cα bond (φ angle)
    or C-Cα bond (ψ angle).

    Cα is the carbon atom connected to the R group.

    In our case, we use the -Cα-O-Cα- plane as PetidePlane
    """
    # get GetPeptidePlane
    v1 = positions[resi2.Ca] - positions[resi1.Ca]
    v2 = positions[resi1.O] - positions[resi1.Ca]
    forward = v1/np.linalg.norm(v1)
    position = (positions[resi1.Ca] + positions[resi2.Ca])/2
    # position = positions[resi1.Ca]
    normal = np.cross(v1, v2)
    normal = normal/np.linalg.norm(normal)
    side = np.cross(normal, forward)
    # tilt = np.arccos(side[2])
    flipped = False
    plane = {
        'position': position,
        'forward': forward,
        'normal': normal,
        'side': side,
        'flipped': flipped,
    }
    return plane


def GetBackbone():
    """
    the repeating -Cα-C-N-Cα- is the backbone of the peptide chain.
    """
