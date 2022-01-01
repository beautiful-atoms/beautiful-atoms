"""
https://en.wikipedia.org/wiki/Ribbon_diagram

"""

import bpy
import numpy as np
from time import time
from batoms.base import Setting

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
    spline = crv.splines.new(type='BEZIER')
    crv.resolution_u = 30
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
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        return bpy.data.collections.get('%s_ribbon'%self.label)

    def build_sheet(self):
        """
        
        """
        arrays = self.batoms.arrays
        sheet_datas = {}
        indices_C = np.where((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA'))[0]
        for data in self.collection:
            sheet_datas[data.name] = data.as_dict()
            istart = np.where((arrays['chainids'] == data.startChain) & (arrays['residuenumbers'] == data.startResi))[0][0]
            iend = np.where((arrays['chainids'] == data.endChain) & (arrays['residuenumbers'] == data.endResi))[0][-1]
            indices = np.where((indices_C >= istart) & (indices_C <= iend))[0]
            # print(data.name, istart, iend, indices)
            positions = arrays['positions'][indices_C[indices]]
            sheet_datas[data.name]['vertices'] = positions
        # print(sheet_datas)
        self.sheet_datas = sheet_datas
        return sheet_datas

    def draw_sheet(self):
        self.build_sheet()
        for name, data in self.sheet_datas.items():
            draw_sheet_from_vertices(name, data, self.coll)
        self.batoms.set_hide(True, only_atoms = True)


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
    
    @property
    def coll(self):
        return self.get_coll()
    
    def get_coll(self):
        return bpy.data.collections.get('%s_ribbon'%self.label)

    def build_helix(self):
        """
        
        """
        arrays = self.batoms.arrays
        helix_datas = {}
        indices_C = np.where((arrays['atomtypes'] == 'C1') | (arrays['atomtypes'] == 'CA'))[0]
        # print(indices_C)
        for data in self.collection:
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
        return helix_datas

    def draw_helix(self):
        self.build_helix()
        for name, data in self.helix_datas.items():
            draw_sheet_from_vertices(name, data, self.coll)

class Ribbon():
    """
    """
    def __init__(self, label, batoms = None,) -> None:
        self.label = label
        self.sheet = Sheet(label, batoms)
        self.helix = Helix(label, batoms)