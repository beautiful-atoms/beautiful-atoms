"""Definition of the parent plugin class.

"""

import bpy
from time import time
import numpy as np
from batoms.base.collection import Setting
from batoms.base.object import BaseObject
import logging
logger = logging.getLogger(__name__)


class PluginSettings(Setting):
    def __init__(self, label, parent) -> None:
        """Plugin Settings
        This object store the plugin setting.

        Parameters:

        label: str
            The label define the batoms object that
            a Setting belong to.
        parent: Plugin Class

        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.parent = parent
        self.name = 'B{}'.format(parent.name)



class PluginObject(BaseObject):
    def __init__(self,
                 label,
                 name='plugin',
                 batoms=None,
                 ):
        """Plugin Class

        Args:
            label (string): _description_.
            batoms (Batoms, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        self.name = name
        obj_name = '%s_%s' % (label, name)
        BaseObject.__init__(self, obj_name)
        self.settings = PluginSettings(
            self.label, parent=self)

    def build_materials(self, name, color,
                        node_inputs=None,
                        material_style='default',
                        backface_culling = False,
                        vertex_color = None,
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
                              backface_culling=backface_culling,
                              vertex_color=vertex_color)
        return mat

    def draw(self):
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old object
        clean_coll_object_by_type(self.batoms.coll, 'PLUGIN')
        for setting in self.settings.bpy_setting:
            self.draw_item(setting)


    def draw_item(self, setting):
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
        for setting in self.settings.bpy_setting:
            name = '{}_{}_{}'.format(self.label, self.name, setting.name)
            obj = bpy.data.objects.get(name)
            objs[name] = obj
        return objs

    @property
    def mat(self):
        return bpy.data.materials.get(self.name)

    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data

    @property
    def show(self):
        return self.get_show()

    @show.setter
    def show(self, state):
        self.set_show(state)

    def get_show(self):
        return self.settings.bpy_data.show

    def set_show(self, state):
        self.settings.bpy_data.show = state
        for name, obj in self.objs.items():
            obj.hide_render = not state
            obj.hide_set(not state)
