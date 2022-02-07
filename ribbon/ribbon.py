"""
https://en.wikipedia.org/wiki/Ribbon_diagram

https://behreajj.medium.com/scripting-curves-in-blender-with-python-c487097efd13


Carson, M.; Bugg, C. E. (1986), "Algorithm for Ribbon Models of Proteins", Journal of Molecular Graphics, 4 (2): 121-122, doi:10.1016/0263-7855(86)80010-8.


"""

import profile
from time import time

import bpy
import numpy as np
from batoms.base import Setting
from batoms.ribbon.profile import ellipse, rectangle, build_mesh

def GetPeptidePlane(resi1, resi2, positions):
    """
    Peptide plane: The atoms of the group, O=C-N-H, are fixed on the same plane, 
    known as the peptide plane. The whole plane may rotate 
    around the N-Cα bond (φ angle) or C-Cα bond (ψ angle). 

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

def draw_rope_from_vertices(name, data, coll, bevel_control,
                node_type = 'Principled BSDF', 
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
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
    #

    # Create bevel control curve.
    # bpy.ops.curve.primitive_bezier_circle_add(radius=data['radius'], enter_editmode=True)
    # bevel_control = bpy.context.active_object
    # bevel_control.data.name = bevel_control.name = '%s_bevel'%name
    # Set the main curve's bevel control to the bevel control curve.
    crv.bevel_object = bevel_control
    crv.bevel_mode = 'OBJECT'
    # bpy.ops.object.mode_set(mode='OBJECT')
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


def draw_rope_from_vertices_nurbs(name, data, coll, bevel_control,
                node_type = 'Principled BSDF', 
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 10
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP' #'TANGENT' #'Z_UP'
    crv.twist_smooth = 1
    spline = crv.splines.new(type='NURBS')
    # spline.use_endpoint_u = True
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    tilts = data['tilts']
    spline.points.foreach_set('co', vertices)
    # spline.points.foreach_set('tilt', tilts)
    # spline.order_u = 4
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
    crv.resolution_u = 10
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP' # 'Z_UP'
    crv.twist_smooth = 1
    spline = crv.splines.new(type='BEZIER')
    nvert = len(vertices)
    spline.bezier_points.add(nvert-1)
    vertices = vertices.reshape(-1, 1)
    tilts = data['tilts']
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
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True):
    """
    """
    from batoms.material import create_material
    vertices = data['vertices']
    crv = bpy.data.curves.new(name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = 10
    crv.fill_mode = 'FULL'
    crv.use_fill_caps = True
    crv.twist_mode = 'Z_UP' #''
    crv.twist_smooth = 1
    spline = crv.splines.new(type='NURBS')
    # spline.use_endpoint_u = True
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    tilts = data['tilts']
    spline.points.foreach_set('co', vertices)
    # spline.points.foreach_set('tilt', tilts)
    # spline.order_u = 4
    # print('order_u: ', spline.order_u)
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


def curve2mesh(name, vertices, resolution_u = 20, coll = None):
    """
    Use Blender's curve to get b-spline
    resolution: number of vertices = (n(control points) - 1)*preview_U
    default preview_U  is 12
    """
    tstart = time()
    crv = bpy.data.curves.new('%s-curve'%name, 'CURVE')
    crv.dimensions = '3D'
    crv.resolution_u = resolution_u
    crv.fill_mode = 'FULL'
    spline = crv.splines.new(type='NURBS')
    # spline.use_endpoint_u = True
    nvert = len(vertices)
    spline.points.add(nvert-1)
    vertices = np.append(vertices, np.ones((nvert, 1)), axis = 1)
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
                node_type = 'Principled BSDF', 
                node_inputs = None, 
                material_style = 'plastic', 
                backface_culling = True,
                shade_smooth = True,):
    """
    """
    from batoms.material import create_material
    from batoms.ribbon.profile import build_mesh
    vertices = curve2mesh('%s-vertices'%name, 
                data['vertices'], data['resolution'], coll)
    sides = curve2mesh('%s-sides'%name, 
                data['sides'], data['resolution'])
    normals = curve2mesh('%s-normals'%name, 
                data['normals'], data['resolution'])
    profiles = data['profiles']
    scales = data['scales']
    vertices, faces = build_mesh(vertices, normals, sides, profiles, scales)
    me = bpy.data.meshes.new(name)
    # edges = np.zeros((n - 1, 2), dtype=int)
    # edges[:, 0] = np.arange(n - 1)
    # edges[:, 1] = np.arange(1, n)
    # print(edges)
    # print(faces)
    me.from_pydata(vertices, [], faces)
    if shade_smooth:
        me.polygons.foreach_set('use_smooth', [True]*len(me.polygons))
        # bpy.context.view_layer.objects.active = obj
        # bpy.ops.object.shade_smooth()
    # crv.bevel_depth = data['depth']
    obj = bpy.data.objects.new(name, me)
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
        -1: water or ...
    """
    def __init__(self, name, resName, chainID, resSeq, type = 0) -> None:
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

class Sheet():
    """
    """
    def __init__(self, name, startChain, endChain,
                startResi, endResi,
                ) -> None:
        self.name = name
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
        scales[-self.resolution:, 0] *= np.linspace(1.5, 0.2, self.resolution)
        return scales
    
    def as_dict(self):
        """
         A complication arises when the direction of the carbonyl 
        oxygen flips, as is always the case between adjacent residues of B sheets.
        """
        n = len(self.residues)
        for i in range(1, n):
            if np.dot(self.residues[i - 1].plane['side'], self.residues[i].plane['side']) < 0:
                self.residues[i].plane['flipped'] = True
                self.residues[i].plane['normal'] = -self.residues[i].plane['normal']
                self.residues[i].plane['side'] = -self.residues[i].plane['side']
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

class Helix():
    """
    """
    def __init__(self, name, startChain, endChain,
                startResi, endResi,
                ) -> None:
        self.name = name
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
            positions[i] = self.residues[i].plane['position'] + \
                self.residues[i].plane['normal']*1.5
        return positions
    
    @property
    def profiles(self):
        profiles = ellipse(32, self.extrude, self.depth)
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
            tilts[i] = -np.arccos(self.residues[i].plane['side'][2])
        return tilts
    
    @property
    def scales(self):
        n = len(self.residues)
        scales = np.ones((self.resolution*(n-1), 3))
        return scales
    
    def as_dict(self):
        data = {
                'vertices': self.positions,
                'color': self.color,
                'tilts': self.tilts,
                'sides': self.sides,
                'normals': self.normals,
                'profiles': self.profiles,
                'scales': None,
                'depth': self.depth,
                'resolution': self.resolution,
                }
        return data


class Turn():
    """
    """
    def __init__(self, name, startChain,
                startResi, 
                ) -> None:
        self.name = name
        self.startChain = startChain
        self.startResi = startResi
        self.residues = []
        self.indices = []
        self.Ca = []
        self.C = []
        self.O = []
        self.extrude = 1.0
        self.depth = 0.3
        self.color = [0, 1, 0, 1]
    
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
            if self.residues[i].type == 2:
                positions[i] = self.residues[i].plane['position'] + \
                self.residues[i].plane['normal']*1.5
            else:
                positions[i] = self.residues[i].plane['position']
        return positions
    
    @property
    def tilts(self):
        n = len(self.residues)
        tilts = np.zeros(n)
        for i in range(n):
            tilts[i] = -np.arccos(self.residues[i].plane['side'][2])
        return tilts
    
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

class HelixSetting(Setting):
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

class TurnSetting(Setting):
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
        self.sheetsetting = SheetSetting(label, batoms)
        self.helixsetting = HelixSetting(label, batoms)
        self.turnsetting = TurnSetting(label, batoms)
        self.datas = datas
        self.import_data(datas)
        if update:
            self.update()
    
    def update(self):
        tstart = time()
        arrays = self.batoms.arrays
        if 'types' not in arrays:
            return
        # build Chains
        Chains = {}
        chainIDs = np.unique(arrays['chainids'])
        for chainID in chainIDs:
            Chains[chainID] = {}
            indices = np.where(arrays['chainids'] == chainID)[0]
            Chains[chainID]['indices'] = indices
        # build residues
        residues = []
        # get Sheet
        names = np.core.defchararray.add(arrays['chainids'], arrays['residuenumbers'].astype('U20'))
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
            plane = GetPeptidePlane(residues[i], residues[i + 1], arrays['positions'])
            residues[i].plane = plane
        # init sheets
        sheets = {}
        for sheet in self.sheetsetting.collection:
            sheet1 = Sheet(sheet.name, sheet.startChain,
                                sheet.endChain,
                                sheet.startResi,
                                sheet.endResi)
            sheet1.color = sheet.color
            sheets[sheet.name] = sheet1
        # init helixs
        helixs = {}
        for helix in self.helixsetting.collection:
            helix1 = Helix(helix.name, helix.startChain,
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
            turn = Turn(residues[0].name, residues[0].chainID,
                                residues[0].resSeq,
                                )
            turns[turn.name] = turn
            turns[turn.name].append(residues[0])
        # add turn
        for i in range(1, nresi - 1):
            if residues[i-1].chainID != turn.startChain or  \
                (residues[i-1].type !=0 and residues[i].type == 0):
                turn = Turn(residues[i].name, residues[i].chainID,
                            residues[i].resSeq,
                            )
                turns[turn.name] = turn
                if i > 1:
                    turns[turn.name].append(residues[i - 2])
                turns[turn.name].append(residues[i - 1])
                turns[turn.name].append(residues[i])
            elif residues[i].type == 0 and residues[i + 1].type ==0:
                turns[turn.name].append(residues[i])
            elif residues[i].type == 0 and residues[i + 1].type !=0:
                turns[turn.name].append(residues[i])
                if i < nresi -1:
                    turns[turn.name].append(residues[i + 1])
                    if i < nresi -2:
                        turns[turn.name].append(residues[i + 2])
            
        self.residues = residues
        self.Chains = Chains
        self.sheets = sheets
        self.helixs = helixs
        self.turns = turns
        self.batoms.selects['sel0'].show = False
        # HETATM    
        indices = np.where(arrays['types'] == 'HETATM')[0]
        self.batoms.selects.add('heta', indices)

        print('update ribbon: %s'%(time() - tstart))

    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        return bpy.data.collections.get('%s_ribbon'%self.label)

    def import_data(self, datas):
        if 'sheet' in datas:
            self.sheetsetting.from_dict(datas['sheet'])
        if 'helix' in datas:
            self.helixsetting.from_dict(datas['helix'])
    
    def build_sheet(self):
        """
        
        """
        tstart = time()
        print('build_sheet: %s'%(time() - tstart))

    def draw_sheet(self):
        tstart = time()
        for name, sheet in self.sheets.items():
            # draw_sheet_from_vertices('sheet-%s'%name, sheet.as_dict(), self.coll)
            # draw_sheet_from_vertices_nurbs('sheet-%s'%name, sheet.as_dict(), self.coll)
            draw_sheet_from_vertices_spline('sheet-%s'%name, 
                    sheet.as_dict(), self.coll, shade_smooth = False)
        print('draw sheet: %s'%(time() - tstart))
        self.batoms.selects['sel0'].show = False
    
    def draw_helix(self):
        tstart = time()
        for name, helix in self.helixs.items():
            # draw_sheet_from_vertices('helix-%s'%name, helix.as_dict(), self.coll)
            # draw_sheet_from_vertices_nurbs('helix-%s'%name, helix.as_dict(), self.coll)
            draw_sheet_from_vertices_spline('helix-%s'%name, 
                        helix.as_dict(), self.coll, shade_smooth = True)
        print('draw helix: %s'%(time() - tstart))
        self.batoms.selects['sel0'].show = False
    
    def draw_turn(self):
        # self.build_turn()
        # Create bevel control curve.
        tstart = time()
        bpy.ops.curve.primitive_bezier_circle_add(radius=0.25, enter_editmode=False)
        bevel_control = bpy.context.active_object
        bevel_control.data.name = bevel_control.name = '%s_turn_bevel'%self.label
        for name, turn in self.turns.items():
            # draw_rope_from_vertices('turn-%s'%name, turn.as_dict(), self.coll, bevel_control)
            draw_rope_from_vertices_nurbs('turn-%s'%name, turn.as_dict(), self.coll, bevel_control)
        print('draw turn: %s'%(time() - tstart))
        self.batoms.selects['sel0'].show = False
    
    def draw(self):
        tstart = time()
        self.draw_sheet()
        self.draw_helix()
        self.draw_turn()
        print('draw ribbon: %s'%(time() - tstart))
