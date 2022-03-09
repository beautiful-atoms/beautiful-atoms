"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting
from batoms.utils.butils import clean_coll_objects, object_mode
from batoms.draw import draw_surface_from_vertices
from time import time

default_colors = [(1, 1, 0, 0.8), (0.0, 0.0, 1.0, 0.8)]

class IsosurfaceSetting(Setting):
    def __init__(self, label, batoms=None, isosurfacesetting=None) -> None:
        """IsosurfaceSetting object

            The IsosurfaceSetting object store the isosurface information.

            Parameters:

            label: str
                The label define the batoms object that a Setting belong to.


        Args:
            label (_type_):
                _description_
            batoms (_type_, optional):
                _description_. Defaults to None.
            isosurfacesetting (_type_, optional):
                _description_. Defaults to None.
        """
        Setting.__init__(self, label, coll_name='%s_volume' % label)
        self.label = label
        self.batoms = batoms
        self.name = 'bisosurface'
        # add a default level
        if isosurfacesetting is not None:
            for key, data in isosurfacesetting.items():
                self[key] = data
        volume = self.batoms.volume
        if volume is not None:
            if len(self) == 0:
                self['1'] = {'level': volume.max()/8, 'color': [1, 1, 0, 0.8]}

    def set_default(self):
        """
        """
        for sp, data in self.species.items():
            self[sp] = {'color': np.append(
                data['color'][:3], 0.3), 'level': 0.005}

    def add(self, isosurfacepair):
        for key in isosurfacepair:
            self.set_default(key)

    def remove_isosurfaces(self, isosurfacepair):
        for key in isosurfacepair:
            name = '%s-%s' % (key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Center     level                color  \n"
        for iso in self.collection:
            s += "{0:10s}   {1:1.6f}  ".format(iso.name, iso.level)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}] \n".format(
                iso.color[0], iso.color[1], iso.color[2], iso.color[3])
        s += "-"*60 + "\n"
        return s

    def build_isosurface(self, cell):
        volume = self.batoms.volume
        isosurface = {}
        for iso in self.collection:
            name = iso.name
            level = iso.level
            color = iso.color
            verts, faces = calc_isosurface(volume, cell, level)
            isosurface[name] = {'vertices': verts,
                                'edges': [],
                                'faces': faces,
                                'color': color,
                                'battr_inputs': {'isosurface': iso.as_dict()}
                                }
        return isosurface

    def draw_isosurface(self):
        """Draw isosurface.
        """
        object_mode()
        clean_coll_objects(self.coll, 'isosurface')
        isosurface = self.build_isosurface(self.batoms.cell)
        for name, isosurface_data in isosurface.items():
            name = '%s_%s_%s' % (self.label, 'isosurface', name)
            draw_surface_from_vertices(name,
                                       datas=isosurface_data,
                                       coll=self.coll,
                                       )


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
