"""
"""
import bpy
import numpy as np
from time import time
from batoms.bondsetting import Setting

default_colors = [(1, 1, 0, 0.8), (0.0, 0.0, 1.0, 0.8)]

class SASAsetting(Setting):
    """
    SASA object

    The SASA object store the SASA information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    
    def __init__(self, label, probe = 1.4, SASA = None) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bSASA'
        self.probe = 1.4
        
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'probe = %s \n'%(self.probe)
        s += '-'*60 + '\n'
        return s
    
    def build_voronoi(self, positions, radii, cell):
        voronoi = {}
        probe = self.probe
        color = default_colors[0]
        verts, faces = calc_voronoi(positions, radii, cell)
        voronoi = {'vertices': verts, 
                            'edges': [], 
                            'faces': faces, 
                            'color': color,
                            'battr_inputs': {},
                        }
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
    # datas = datas[0]
    # print(datas)
    faces = []
    vertices = []
    for data in datas:
        # adjacency = data['adjacency']
        nvert = len(vertices)
        vertices.extend(data['vertices'])
        faces1 = [[x + nvert for x in face['vertices']] for face in data['faces']]
        faces.extend(faces1)
    print('calc_voronoi: %s'%(time() - tstart))
    return vertices, faces






