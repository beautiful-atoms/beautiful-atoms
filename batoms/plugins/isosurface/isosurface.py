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
        self.settings.bpy_data.active = True

    def build_isosurface(self, cell):
        isosurface = {}
        for iso in self.settings.bpy_setting:
            if iso.volumetric_data == '':
                volume = self.batoms.volumetric_data[self.batoms.volumetric_data.bpy_setting[0].name]
            else:
                volume = self.batoms.volumetric_data[iso.volumetric_data]
            scaled_verts, faces = calc_isosurface(volume, cell, iso.level)
            # color by another volumetric data
            if iso.color_by != "None":
                from batoms.utils import map_volumetric_data
                data = map_volumetric_data(
                                self.batoms.volumetric_data[iso.color_by],
                                scaled_verts
                                )
                attribute_data = (data - np.min(data))/(np.max(data) - np.min(data))
                color_by_attribute = {'attribute_name': '{}_data'.format(iso.color_by),
                                'ValToRGB':[iso.color1[:],
                                            iso.color2[:]]
                                            }
            else:
                color_by_attribute = None
                attribute_data = None

            # back to cell
            verts = scaled_verts.dot(cell)
            # verts -= self.batoms.cell.origin
            isosurface[iso.name] = {'vertices': verts,
                                    'edges': [],
                                    'faces': faces,
                                    'material_style': iso.material_style,
                                    'color': iso.color,
                                    'color_by': iso.color_by,
                                    'color_by_attribute': color_by_attribute,
                                    'battr_inputs': {'Bisosurface': iso.as_dict()},
                                    'attribute_data': attribute_data,
                                    }
        return isosurface

    def build_materials(self, name, color, node_inputs=None,
                        material_style='default',
                        color_by_attribute=None):
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
                              backface_culling=False,
                              color_by_attribute=color_by_attribute)
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
            if isosurface_data['attribute_data'] is not None:
                from batoms.utils.attribute import set_mesh_attribute
                obj.data.attributes.new(name=isosurface_data['color_by_attribute']['attribute_name'],
                                type='FLOAT', domain='POINT')
                set_mesh_attribute(obj, isosurface_data['color_by_attribute']['attribute_name'], isosurface_data['attribute_data'])
            # material
            mat = self.build_materials(name, isosurface_data['color'],
                                       material_style=isosurface_data['material_style'],
                                        color_by_attribute = isosurface_data['color_by_attribute'],
                                       )
            obj.data.materials.append(mat)

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated(
            '"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data


def calc_isosurface(volumetric_data, cell, level,
                    gradient_direction='descent',
                    step_size=1):
    """

    Computes an isosurface from a volumetric_data grid.

    Parameters:

    volumetric_data: np.array

    cell: np.array

    level: float

    """
    from skimage import measure

    spacing = tuple(1.0/np.array(volumetric_data.shape))
    mlevel = np.mean(volumetric_data)
    if not level:
        level = mlevel*10
    # print('iso level: {0:1.9f}, iso mean: {1:1.9f}'.format(level, mlevel))
    scaled_verts, faces, normals, values = \
        measure.marching_cubes(volumetric_data, level=level,
                               spacing=spacing,
                               gradient_direction=gradient_direction,
                               allow_degenerate=False,
                               step_size=step_size)
    faces = list(faces)
    return scaled_verts, faces
