"""Definition of the isosurface class.

This module defines the isosurface object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.base.object import BaseObject
from .setting import IsosurfaceSettings
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class Isosurface(BaseObject):
    def __init__(self,
                 label=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 ):
        """Isosurface Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = 'isosurface'
        BaseObject.__init__(self, label, name)
        self.settings = IsosurfaceSettings(
            self.label, batoms=batoms, parent=self)

    def build_isosurface(self, cell):
        volume = self.batoms.volume
        isosurface = {}
        for iso in self.settings.bpy:
            verts, faces = calc_isosurface(volume, cell, iso.level)
            isosurface[iso.name] = {'vertices': verts,
                                    'edges': [],
                                    'faces': faces,
                                    'material_style': iso.material_style,
                                    'color': iso.color,
                                    'battr_inputs': {'Bisosurface': iso.as_dict()}
                                    }
        return isosurface

    def build_materials(self, name, color, node_inputs=None,
                        material_style='default'):
        """
        """
        from batoms.material import create_material
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                              color=color,
                              node_inputs=node_inputs,
                              material_style=material_style,
                              backface_culling=False)
        return mat

    def draw(self, isosurface_name='ALL'):
        """Draw isosurface.
        """
        from batoms.draw import draw_surface_from_vertices
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old isosurface
        clean_coll_object_by_type(self.batoms.coll, 'ISOSURFACE')
        isosurface = self.build_isosurface(self.batoms.cell)
        for name, isosurface_data in isosurface.items():
            if isosurface_name.upper() != "ALL" and name.name != isosurface_name:
                continue
            name = '%s_%s_%s' % (self.label, 'isosurface', name)
            self.delete_obj(name)
            obj = draw_surface_from_vertices(name,
                                             datas=isosurface_data,
                                             coll=self.batoms.coll,
                                             )
            obj.batoms.type = 'ISOSURFACE'
            obj.batoms.label = self.label
            obj.parent = self.batoms.obj
            # material
            mat = self.build_materials(name, isosurface_data['color'],
                                       material_style=isosurface_data['material_style'],
                                       )
            obj.data.materials.append(mat)

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

def calc_isosurface(volume, cell, level,
                    gradient_direction='descent',
                    step_size=1):
    """

    Computes an isosurface from a volume grid.

    Parameters:

    volume: np.array

    cell: np.array

    level: float

    """
    from skimage import measure

    cell_origin = cell.origin
    #
    spacing = tuple(1.0/np.array(volume.shape))
    mlevel = np.mean(volume)
    if not level:
        level = mlevel*10
    # print('iso level: {0:1.9f}, iso mean: {1:1.9f}'.format(level, mlevel))
    scaled_verts, faces, normals, values = \
        measure.marching_cubes(volume, level=level,
                               spacing=spacing,
                               gradient_direction=gradient_direction,
                               allow_degenerate=False,
                               step_size=step_size)
    scaled_verts = scaled_verts.dot(cell)
    scaled_verts -= cell_origin
    faces = list(faces)
    return scaled_verts, faces
