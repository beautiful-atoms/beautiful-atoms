"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

class HighlightSettings(Setting):
    def __init__(self, label, batoms=None, parent=None,
                 Highlightsetting=None) -> None:
        """HighlightSetting object

            The HighlightSetting object store the highlight information.

            Parameters:

            label: str
                The label define the batoms object that a Setting belong to.


        Args:
            label (_type_):
                _description_
            batoms (_type_, optional):
                _description_. Defaults to None.
            Highlightsetting (_type_, optional):
                _description_. Defaults to None.
        """
        Setting.__init__(self, label, coll_name=label)
        self.label = label
        self.batoms = batoms
        self.parent = parent
        self.name = 'Bhighlight'
        # add a default level
        if Highlightsetting is not None:
            for key, data in Highlightsetting.items():
                self[key] = data

    def add(self, highlight):
        if isinstance(highlight, str):
            self[highlight] = {}
        elif isinstance(highlight, dict):
            self[highlight['name']] = highlight

    def remove_highlights(self, highlight):
        for key in highlight:
            name = '%s-%s' % (key[0], key[1])
            i = self.bpy_setting.find(name)
            if i != -1:
                self.bpy_setting.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Name       min   max   scale                color  \n"
        for cav in self.bpy_setting:
            s += "{:4s}   {:1.3f}   {:1.3f}   {:1.3f}   " \
                .format(cav.name, cav.min, cav.max, cav.scale)
            s += "[{:1.2f}  {:1.2f}  {:1.2f}  {:1.2f}] \n".format(
                cav.color[0], cav.color[1], cav.color[2], cav.color[3])
        s += "-"*60 + "\n"
        return s

    def __setitem__(self, name, setdict):
        """
        Set properties
        """
        name = str(name)
        subset = self.find(name)
        if subset is None:
            subset = self.bpy_setting.add()
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_instancer(subset.as_dict())

    def build_instancer(self, hl):
        """Build object instancer for hl

        Returns:
            bpy Object: instancer
        """
        name = 'highlight_%s_%s' % (self.label, hl['name'])
        self.delete_obj(name)
        print(hl)
        if hl["style"] == "1":
            bpy.ops.mesh.primitive_cube_add(size=1)
        else:
            bpy.ops.mesh.primitive_uv_sphere_add(radius=1)
        obj = bpy.context.view_layer.objects.active
        obj.name = name
        obj.data.name = name
        obj.batoms.type = 'INSTANCER'
        #
        obj.users_collection[0].objects.unlink(obj)
        bpy.data.collections['%s_instancer' % self.label].objects.link(obj)
        bpy.ops.object.shade_smooth()
        obj.hide_set(True)
        obj.hide_render = True
        mat = self.build_materials(hl)
        obj.data.materials.append(mat)
        bpy.context.view_layer.update()
        self.parent.add_geometry_node(hl['name'], obj)
        return obj

    def build_materials(self, hl, node_inputs=None):
        """
        """
        from batoms.material import create_material
        name = 'highlight_%s_%s' % (
            self.label, hl['name'])
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                        color=hl['color'],
                        node_inputs=node_inputs,
                        material_style=hl['material_style'],
                        backface_culling=True)
        return mat
