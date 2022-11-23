"""Definition of the Magres class.

This module defines the Magres object in the Batoms package.

"""

import bpy
import bmesh
from time import time
import numpy as np
from batoms.base.object import BaseObject
from .setting import MagresSettings
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


class Magres(BaseObject):
    def __init__(self,
                 label=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 tensors = None,
                 ):
        """Magres Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        self.tensors = tensors
        name = 'magres'
        BaseObject.__init__(self, label, name)
        self.settings = MagresSettings(
            self.label, parent=self)
        self.settings.bpy_data.active = True

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

    def draw(self, magres_name="ALL"):
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old plane
        clean_coll_object_by_type(self.batoms.coll, 'MAGRES')
        for magres in self.settings.bpy_setting:
            if magres_name.upper() != "ALL" and magres.name != magres_name:
                continue
            if magres.type == 'MS':
                self.draw_MS(magres)
            elif magres.type == 'CS':
                self.draw_CS(magres)

    @property
    def ms_objs(self):
        ms_objs = {}
        for magres in self.settings.bpy_setting:
            ms_name = '%s_%s_ms' % (self.label, magres.name)
            obj = bpy.data.objects.get(ms_name)
            ms_objs[magres.name] = obj
        return ms_objs

    @property
    def ms_mat(self):
        return bpy.data.materials.get(self.ms_name)

    @property
    def cs_objs(self):
        cs_objs = {}
        for magres in self.settings.bpy_setting:
            ms_name = '%s_%s_cs' % (self.label, magres.name)
            obj = bpy.data.objects.get(ms_name)
            cs_objs[magres.name] = obj
        return cs_objs

    def build_ellipsoids(self, positions, tensors, scale):
        """
        generate ellipsoid
        """
        from numpy.linalg import eig
        n = len(tensors)
        ellipsoids = []
        for i in range(n):
            tensor = tensors[i]
            # get symmetry part
            sym_tensor = (tensor + np.transpose(tensor))/2.0
            evals, evecs = eig(sym_tensor)
            # abs and normlize
            evals = np.abs(evals)
            evals = evals*scale
            # sort
            index = np.argsort(evals)
            evals = evals[index]
            evecs = evecs[:, index]
            # print(evecs)
            # evecs = evecs/np.linalg.norm(evecs, axis = 0)
            # rotation matrix
            rotation_matrix = np.linalg.inv(evecs)
            ellipsoids.append([positions[i], evals, rotation_matrix])
        return ellipsoids

    def draw_ellipsoids(self, name, datas, coll):
        objs = []
        for data in datas:
            # print(data)
            p = np.eye(4)
            bpy.ops.mesh.primitive_uv_sphere_add(location = data[0], scale = data[1])
            obj = bpy.context.view_layer.objects.active
            me = obj.data
            # rotatation
            p[0:3, 0:3] = data[2]
            me.transform(p)
            me.update()
            objs.append(obj)
        objs[0].name = name
        objs[0].data.name = name
        c = {}
        c["object"] = c["active_object"] = objs[0]
        c["selected_objects"] = c["selected_editable_objects"] = objs
        # bpy.context.view_layer.objects.active = obj
        bpy.ops.object.join(c)
        bpy.ops.object.shade_smooth()
        return bpy.data.objects[name]

    def draw_MS(self, magres):
        """
        """
        scale = magres.scale
        # print('Scale: {:1.3}'.format(scale))
        tstart = time()
        indices = self.batoms.selects[magres.select].indices
        if len(indices) == 0:
            return
        # select atoms
        positions = self.batoms.positions[indices]
        tensors = self.batoms.get_attribute('ms')[indices]
        ellipsoids = self.build_ellipsoids(positions, tensors, scale)
        ms_name = '%s_%s_ms' % (self.label, magres.name)
        self.delete_obj(ms_name)
        coll = self.batoms.coll.children['%s_surface' %
                                         self.batoms.coll_name]
        obj = self.draw_ellipsoids(ms_name,
                            datas=ellipsoids,
                            coll=coll,
                            )
        mat = self.build_materials(ms_name, color=magres.color,
                                   material_style=magres.material_style,
                                   )
        obj.data.materials.append(mat)
        obj.parent = self.batoms.obj
        obj.batoms.type = 'MS'
        obj.batoms.label = self.batoms.label
        logger.debug('Draw MS: %s' % (time() - tstart))

    def draw_CS(self, magres, parallel=1):
        """
        """
        pass

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data
