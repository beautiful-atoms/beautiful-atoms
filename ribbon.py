"""
https://en.wikipedia.org/wiki/Ribbon_diagram

https://behreajj.medium.com/scripting-curves-in-blender-with-python-c487097efd13





"""

from numpy.core.records import array
import bpy
import numpy as np
from time import time
from batoms.base import Setting

def GetPeptidePlane():
    """
    Peptide plane: The atoms of the group, O=C-N-H, are fixed on the same plane, 
    known as the peptide plane. The whole plane may rotate 
    around the N-Cα bond (φ angle) or C-Cα bond (ψ angle). 

    Cα is the carbon atom connected to the R group.

    In our case, we use the -Cα-O-Cα- plane as PetidePlane
    """
    

def GetBackbone():
    """
    the repeating -Cα-C-N-Cα- is the backbone of the peptide chain.
    """

def draw_rope_from_vertices_nurbs(name, data, coll,
                node_type = 'Principled BSDF', 
                use_smooth = True,
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 30
    crv.fill_mode = 'FULL'
    spline = crv.splines.new(type='NURBS')
    nvert = len(vertices)
    spline.points.add(nvert-1)
    # vertices = np.append(vertices, np.zeros((nvert, 1)), axis = 1)
    vertices = np.append(vertices, np.ones((len(vertices), 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    spline.points.foreach_set('co', vertices)
    # Create bevel control curve.
    bpy.ops.curve.primitive_bezier_circle_add(radius=0.25, enter_editmode=True)
    bevel_control = bpy.context.active_object
    bevel_control.data.name = bevel_control.name = '%s_bevel'%name
    # Set the main curve's bevel control to the bevel control curve.
    crv.bevel_object = bevel_control
    crv.bevel_mode = 'OBJECT'
    bpy.ops.object.mode_set(mode='OBJECT')
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name, 
                    data['color'], 
                    node_type = node_type, 
                    node_inputs = node_inputs, 
                    material_style = material_style, 
                    backface_culling = backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)


def draw_rope_from_vertices(name, data, coll,
                node_type = 'Principled BSDF', 
                use_smooth = True,
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 30
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
    #
    # Create bevel control curve.
    bpy.ops.curve.primitive_bezier_circle_add(radius=data['radius'], enter_editmode=True)
    bevel_control = bpy.context.active_object
    bevel_control.data.name = bevel_control.name = '%s_bevel'%name
    # Set the main curve's bevel control to the bevel control curve.
    crv.bevel_object = bevel_control
    crv.bevel_mode = 'OBJECT'
    bpy.ops.object.mode_set(mode='OBJECT')
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name, 
                    data['color'], 
                    node_type = node_type, 
                    node_inputs = node_inputs, 
                    material_style = material_style, 
                    backface_culling = backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)

def draw_sheet_from_vertices(name, data, coll,
                node_type = 'Principled BSDF', 
                use_smooth = True,
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 30
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP'
    spline = crv.splines.new(type='BEZIER')
    nvert = len(vertices)
    spline.bezier_points.add(nvert-1)
    # vertices = np.append(vertices, np.zeros((nvert, 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    tilts = data['tilt']
    spline.bezier_points.foreach_set('co', vertices)
    spline.bezier_points.foreach_set('tilt', tilts)
    for i in range(nvert):
        spline.bezier_points[i].handle_right_type = 'AUTO'
        spline.bezier_points[i].handle_left_type = 'AUTO'
    crv.extrude = data['extrude']
    crv.bevel_depth = data['depth']
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name, 
                    data['color'], 
                    node_type = node_type, 
                    node_inputs = node_inputs, 
                    material_style = material_style, 
                    backface_culling = backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)

def draw_sheet_from_vertices_nurbs(name, data, coll,
                node_type = 'Principled BSDF', 
                use_smooth = True,
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 30
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP'
    spline = crv.splines.new(type='NURBS')
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis = 1)
    vertices[0, 3] = 1000
    vertices[-1, 3] = 1000
    vertices = vertices.reshape(-1, 1)
    tilts = data['tilt']
    spline.points.foreach_set('co', vertices)
    spline.points.foreach_set('tilt', tilts)
    crv.extrude = data['extrude']
    crv.bevel_depth = data['depth']
    obj = bpy.data.objects.new(name, crv)
    #
    material = create_material(name, 
                    data['color'], 
                    node_type = node_type, 
                    node_inputs = node_inputs, 
                    material_style = material_style, 
                    backface_culling = backface_culling)
    obj.data.materials.append(material)
    coll.objects.link(obj)

class Chain():
    def __init__(self, ChainID, Residues) -> None:
        self.ChainID = ChainID
        self.Residues = Residues

class Residue():
    """
    Type:
        0:
        1: Sheet
        2: Helix
    """
    def __init__(self, ResName, ChainID, ResSeq, Atoms, Type = 0) -> None:
        self.ResName = ResName
        self.ChainID = ChainID
        self.ResSeq = ResSeq
        self.Atoms = Atoms
        self.Type = Type  

class Sheet(Setting):
    """
    sheet = {
    'name': '%s-%s'%(chainId, sheetId),
    'startChain': startChain,
    'startResi': startResi,
    'endChain': endChain,
    'endResi': endResi,
    }
    """
    def __init__(self, label, batoms = None,
                ) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bsheet'
        self.batoms = batoms
        self.sheet_name = '%s_sas'%self.label
        self.sheet_datas = {}
        
    @property
    def show(self):
        return self.get_show()
    
    @show.setter
    def show(self, state):
        self.set_show(state)
    
    def get_show(self):
        return self.batoms.attributes['show'][self.indices]
    
    def set_show(self, show, only_atoms = True):
        show0 = self.batoms.show
        show0[self.indices] = show
        self.batoms.set_attributes({'show': show0})

class Helix(Setting):
    """
    helix = {
    'name': '%s-%s'%(chainId, helixId),
    'startChain': startChain,
    'startResi': startResi,
    'endChain': endChain,
    'endResi': endResi,
    }
    """
    def __init__(self, label, batoms = None,
                ) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bhelix'
        self.batoms = batoms
        self.helix_name = '%s_sas'%self.label
        self.helix_datas = {}

class Turn(Setting):
    """
    turn = {
    'name': '%s-%s'%(chainId, turnId),
    'startChain': startChain,
    'startResi': startResi,
    'endChain': endChain,
    'endResi': endResi,
    }
    """
    def __init__(self, label, batoms = None,
                ) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bturn'
        self.batoms = batoms
        self.turn_name = '%s_sas'%self.label
        self.turn_datas = {}

class Ribbon():
    """
    """
    def __init__(self, label, batoms = None, datas = {}, update = False) -> None:
        self.label = label
        self.batoms = batoms
        self.sheet = Sheet(label, batoms)
        self.helix = Helix(label, batoms)
        self.turn = Turn(label, batoms)
        self.datas = datas
        self.import_data(datas)
        if update:
            self.update()
    
    def update(self):
        """
        """
        arrays = self.batoms.arrays
        if 'types' not in arrays:
            return
        # get Ca
        self.calpha = np.where((arrays['types'] == 'ATOM') & 
                    ((arrays['atomtypes'] == 'C1') | 
                    (arrays['atomtypes'] == 'CA')))[0]
        # get chain
        Chains = {}
        chainIDs = np.unique(arrays['chainids'])
        for chainID in chainIDs:
            Chains[chainID] = {}
            indices = np.where(arrays['chainids'] == chainID)[0]
            # get residue
            residueid = np.unique(arrays['residuenumbers'][indices])
            # Ca
            Chains[chainID]['indices'] = indices
            Chains[chainID]['residueid'] = residueid
            Chains[chainID]['Ca'] = self.calpha[np.where((self.calpha > indices[0]) & (self.calpha < indices[-1]))[0]]
            Chains[chainID]['C'] = Chains[chainID]['Ca'] + 1
            Chains[chainID]['O'] = Chains[chainID]['Ca'] + 2
            # get GetPeptidePlane
            v1 = arrays['positions'][Chains[chainID]['C']] - arrays['positions'][Chains[chainID]['Ca']]
            v2 = arrays['positions'][Chains[chainID]['O']] - arrays['positions'][Chains[chainID]['Ca']]
            x = v2/np.linalg.norm(v2)
            z = np.cross(v1, v2)
            z = z/np.linalg.norm(v2)
            y = np.cross(z, x)
            tilt = np.arccos(y[:, 2])
            Chains[chainID]['positions'] = arrays['positions'][Chains[chainID]['Ca']]
            Chains[chainID]['normal'] = z
            Chains[chainID]['tilt'] = tilt
        self.Chains = Chains
        self.arrays = arrays
        # get Sheet
        sheet_datas = {}
        for data in self.sheet.collection:
            sheet_datas[data.name] = data.as_dict()
            istart = np.where(self.Chains[data.startChain]['residueid'] == data.startResi)[0][0]
            iend = np.where(self.Chains[data.endChain]['residueid'] == data.endResi)[0][0]
            

    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        return bpy.data.collections.get('%s_ribbon'%self.label)

    def import_data(self, datas):
        if 'sheet' in datas:
            self.sheet.from_dict(datas['sheet'])
        if 'helix' in datas:
            self.helix.from_dict(datas['helix'])
    
    def build_sheet(self):
        """
        
        """
        tstart = time()
        arrays = self.arrays
        sheet_datas = {}
        for data in self.sheet.collection:
            sheet_datas[data.name] = data.as_dict()
            istart = np.where(self.Chains[data.startChain]['residueid'] == data.startResi)[0][0]
            iend = np.where(self.Chains[data.endChain]['residueid'] == data.endResi)[0][0]
            indices = np.where((self.calpha >= istart) & (self.calpha <= iend))[0]
            # print(data.name, istart, iend, indices)
            sheet_datas[data.name]['vertices'] = self.Chains[data.startChain]['positions'][istart:iend+1]
            sheet_datas[data.name]['tilt'] = self.Chains[data.startChain]['tilt'][istart:iend+1]
        # calc normal, tilt

        # print(sheet_datas)
        self.sheet_datas = sheet_datas
        print('build_sheet: %s'%(time() - tstart))

    def draw_sheet(self):
        self.build_sheet()
        for name, data in self.sheet_datas.items():
            # draw_sheet_from_vertices('sheet-%s'%name, data, self.coll)
            draw_sheet_from_vertices_nurbs('sheet-%s'%name, data, self.coll)
        self.batoms.selects['sel0'].show = False
    
    def build_helix(self):
        """
        
        """
        tstart = time()
        arrays = self.arrays
        helix_datas = {}
        # indices_C = np.where((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA'))[0]
        indices_C = np.where((arrays['types'] == 'ATOM') & ((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA')))[0]
        # print(indices_C)
        for data in self.helix.collection:
            helix_datas[data.name] = data.as_dict()
            istart = np.where(self.Chains[data.startChain]['residueid'] == data.startResi)[0][0]
            iend = np.where(self.Chains[data.endChain]['residueid'] == data.endResi)[0][0]
            indices = np.where((self.calpha >= istart) & (self.calpha <= iend))[0]
            # print(data.name, istart, iend, indices)
            helix_datas[data.name]['vertices'] = self.Chains[data.startChain]['positions'][istart:iend+1]
            helix_datas[data.name]['tilt'] = self.Chains[data.startChain]['tilt'][istart:iend+1]
            # print(istart, iend, indices_C[indices])
        # print(helix_datas)
        self.helix_datas = helix_datas
        print('build_helix: %s'%(time() - tstart))

    def draw_helix(self):
        self.build_helix()
        for name, data in self.helix_datas.items():
            draw_sheet_from_vertices('helix-%s'%name, data, self.coll)
            # draw_sheet_from_vertices_nurbs('helix-%s'%name, data, self.coll)
        self.batoms.selects['sel0'].show = False
    
    def build_turn_dict(self):
        """
        
        """
        self.turn.collection.clear()
        arrays = self.arrays
        mask = {}
        for data in self.sheet.collection:
            if data.startChain not in mask:
                mask[data.startChain] = []
            mask[data.startChain].append([data.startResi, data.endResi])
        for data in self.helix.collection:
            if data.startChain not in mask:
                mask[data.startChain] = []
            mask[data.startChain].append([data.startResi, data.endResi])
        # print(mask)
        for chainId in mask:
            chain_indices = np.where(arrays['chainids'] == chainId)
            imin = min(arrays['residuenumbers'][chain_indices])
            imax = max(arrays['residuenumbers'][chain_indices])
            mask[chainId].append([imin, imin])
            mask[chainId].append([imax, imax])
            mask[chainId].sort()
            n = len(mask[chainId])
            for i in range(1, n):
                if mask[chainId][i][0] > mask[chainId][i - 1][1]:
                    data = {'name': '%s-%s-%s-%s'%(chainId, mask[chainId][i - 1][1], chainId, mask[chainId][i][0]),
                        'startChain': chainId,
                        'startResi': mask[chainId][i - 1][1],
                        'endChain': chainId,
                        'endResi': mask[chainId][i][0],
                        }
                    self.turn.from_dict(data)
                    # print(data)
    
    def build_turn(self):
        """
        
        """
        tstart = time()
        self.build_turn_dict()
        arrays = self.arrays
        turn_datas = {}
        indices_C = np.where((arrays['types'] == 'ATOM') & ((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA')))[0]
        imin = min(indices_C)
        imax = max(indices_C)
        # print(indices_C)
        for data in self.turn.collection:
            turn_datas[data.name] = data.as_dict()
            istart = np.where((arrays['chainids'] == data.startChain) & (arrays['residuenumbers'] == data.startResi))[0][0]
            iend = np.where((arrays['chainids'] == data.endChain) & (arrays['residuenumbers'] == data.endResi))[0][-1]
            indices = np.where((indices_C >= istart) & (indices_C <= iend))[0]
            # print(data.name, istart, iend, indices)
            positions = arrays['positions'][indices_C[indices]]
            turn_datas[data.name]['vertices'] = positions
            # print(istart, iend, indices_C[indices])
        # print(turn_datas)
        self.turn_datas = turn_datas
        print('build_turn: %s'%(time() - tstart))

    def draw_turn(self):
        self.build_turn()
        for name, data in self.turn_datas.items():
            draw_rope_from_vertices('turn-%s'%name, data, self.coll)
        self.batoms.selects['sel0'].show = False
    
    def draw(self):
        self.draw_sheet()
        self.draw_helix()
        self.draw_turn()