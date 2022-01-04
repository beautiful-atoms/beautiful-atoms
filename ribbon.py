"""
https://en.wikipedia.org/wiki/Ribbon_diagram

https://behreajj.medium.com/scripting-curves-in-blender-with-python-c487097efd13

"""

import bpy
import numpy as np
from time import time
from batoms.base import Setting


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
    spline = crv.splines.new(type='BEZIER')
    nvert = len(vertices)
    spline.bezier_points.add(nvert-1)
    # vertices = np.append(vertices, np.zeros((nvert, 1)), axis = 1)
    vertices = vertices.reshape(-1, 1)
    spline.bezier_points.foreach_set('co', vertices)
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
    def __init__(self, label, batoms = None, datas = {}) -> None:
        self.label = label
        self.batoms = batoms
        self.sheet = Sheet(label, batoms)
        self.helix = Helix(label, batoms)
        self.turn = Turn(label, batoms)
        self.datas = datas
        self.import_data(datas)
    
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
        arrays = self.batoms.arrays
        sheet_datas = {}
        indices_C = np.where((arrays['types'] == 'ATOM') & ((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA')))[0]
        for data in self.sheet.collection:
            sheet_datas[data.name] = data.as_dict()
            istart = np.where((arrays['chainids'] == data.startChain) & (arrays['residuenumbers'] == data.startResi))[0][0]
            iend = np.where((arrays['chainids'] == data.endChain) & (arrays['residuenumbers'] == data.endResi))[0][-1]
            indices = np.where((indices_C >= istart) & (indices_C <= iend))[0]
            # print(data.name, istart, iend, indices)
            positions = arrays['positions'][indices_C[indices]]
            sheet_datas[data.name]['vertices'] = positions
        # print(sheet_datas)
        self.sheet_datas = sheet_datas

    def draw_sheet(self):
        self.build_sheet()
        for name, data in self.sheet_datas.items():
            draw_sheet_from_vertices('sheet-%s'%name, data, self.coll)
        self.batoms.set_hide(True, only_atoms = True)
    
    def build_helix(self):
        """
        
        """
        arrays = self.batoms.arrays
        helix_datas = {}
        # indices_C = np.where((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA'))[0]
        indices_C = np.where((arrays['types'] == 'ATOM') & ((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA')))[0]
        # print(indices_C)
        for data in self.helix.collection:
            helix_datas[data.name] = data.as_dict()
            istart = np.where((arrays['chainids'] == data.startChain) & (arrays['residuenumbers'] == data.startResi))[0][0]
            iend = np.where((arrays['chainids'] == data.endChain) & (arrays['residuenumbers'] == data.endResi))[0][-1]
            indices = np.where((indices_C >= istart) & (indices_C <= iend))[0]
            # print(data.name, istart, iend, indices)
            positions = arrays['positions'][indices_C[indices]]
            helix_datas[data.name]['vertices'] = positions
            # print(istart, iend, indices_C[indices])
        # print(helix_datas)
        self.helix_datas = helix_datas

    def draw_helix(self):
        self.build_helix()
        for name, data in self.helix_datas.items():
            draw_sheet_from_vertices('helix-%s'%name, data, self.coll)
        self.batoms.set_hide(True, only_atoms = True)
    
    def build_turn_dict(self):
        """
        
        """
        self.turn.collection.clear()
        arrays = self.batoms.arrays
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
        arrays = self.batoms.arrays
        for chainId in mask:
            chain_indices = np.where(arrays['chainids'] == chainId)
            imin = min(self.batoms.arrays['residuenumbers'][chain_indices])
            imax = max(self.batoms.arrays['residuenumbers'][chain_indices])
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
        self.build_turn_dict()
        arrays = self.batoms.arrays
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

    def draw_turn(self):
        self.build_turn()
        for name, data in self.turn_datas.items():
            draw_rope_from_vertices('turn-%s'%name, data, self.coll)
        self.batoms.set_hide(True, only_atoms = True)
    
    def draw(self):
        self.draw_sheet()
        self.draw_helix()
        self.draw_turn()