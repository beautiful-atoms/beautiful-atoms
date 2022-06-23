"""Definition of the plugin class.

This module defines the plugin object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.base.object import BaseObject
from batoms.plugin_template.templatesetting import TemplateSetting
import logging
logger = logging.getLogger(__name__)


class Template(BaseObject):
    def __init__(self,
                 label,
                 batoms=None,
                 tensors = None,
                 ):
        """Template Class

        Args:
            label (string): _description_.
            batoms (Batoms, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        self.tensors = tensors
        name = 'plugin'
        BaseObject.__init__(self, label, name)
        self.setting = TemplateSetting(
            self.label, parent=self)
    
    def build_materials(self, name, color,
                        node_inputs=None,
                        material_style='default',
                        backface_culling = False
                        ):
        """Create materials

        Args:
            name (string): _description_
            color (array, (1x4)): _description_
            node_inputs (dict, optional): _description_. Defaults to None.
            material_style (string, optional): _description_. Defaults to 'default'.

        Returns:
            bpy.type.material: _description_
        """
        from batoms.material import create_material
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                              color=color,
                              node_inputs=node_inputs,
                              material_style=material_style,
                              backface_culling=backface_culling)
        return mat

    def draw(self):
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old object
        clean_coll_object_by_type(self.batoms.coll, 'TEMPLATE')
        for setting in self.setting.collection:
            self.draw_plugin(setting)

    
    def draw_plugin(self, setting):
        """
        """
        tstart = time()
        indices = self.batoms.selects[setting.select].indices
        if len(indices) == 0:
            return
        # select atoms
        positions = self.batoms.positions[indices]
        logger.debug('Draw Template: %s' % (time() - tstart))

    @property
    def objs(self):
        objs = {}
        for setting in self.setting.collection:
            ms_name = '%s_%s' % (self.label, setting.name)
            obj = bpy.data.objects.get(ms_name)
            objs[setting.name] = obj
        return objs

    @property
    def mat(self):
        return bpy.data.materials.get(self.name)
   