"""
"""
import bpy
import numpy as np
from batoms.base.collection import Setting


class CavitySettings(Setting):
    def __init__(self, label, batoms=None, parent=None,
                 cavitysetting=None) -> None:
        """CavitySetting object

            The CavitySetting object store the cavity information.

            Parameters:

            label: str
                The label define the batoms object that a Setting belong to.


        Args:
            label (_type_):
                _description_
            batoms (_type_, optional):
                _description_. Defaults to None.
            cavitysetting (_type_, optional):
                _description_. Defaults to None.
        """
        Setting.__init__(self, label, coll_name='%s_cavity' % label)
        self.label = label
        self.batoms = batoms
        self.parent = parent
        self.name = 'bcavity'
        # add a default level
        if cavitysetting is not None:
            for key, data in cavitysetting.items():
                self[key] = data
        volume = self.batoms.volume
        if volume is not None:
            if len(self) == 0:
                self['1'] = {'level': volume.max()/8, 'color': [1, 1, 0, 0.8]}

    def add(self, cavity):
        if isinstance(cavity, str):
            self[cavity] = {}
        elif isinstance(cavity, dict):
            self[cavity['name']] = cavity

    def remove_cavitys(self, cavity):
        for key in cavity:
            name = '%s-%s' % (key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "Name     min   max              color  \n"
        for cav in self.collection:
            s += "{:10s}   {:1.6f}   {:1.6f}".format(cav.name, cav.min, cav.max)
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
            subset = self.collection.add()
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_instancer(subset.as_dict())

    def build_instancer(self, sp):
        """Build object instancer for species

        Returns:
            bpy Object: instancer
        """
        name = 'cavity_%s_%s' % (self.label, sp['name'])
        self.delete_obj(name)
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
        mat = self.build_materials(sp)
        obj.data.materials.append(mat)
        bpy.context.view_layer.update()
        self.parent.add_geometry_node(sp['name'], obj)
        return obj

    def build_materials(self, sp, node_inputs=None,
                        material_style='default'):
        """
        """
        from batoms.material import create_material
        name = 'cavity_%s_%s' % (
            self.label, sp['name'])
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                        color=sp['color'],
                        node_inputs=node_inputs,
                        material_style=material_style,
                        backface_culling=True)
        return mat