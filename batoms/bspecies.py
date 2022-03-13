import bpy
import numpy as np
from batoms.base.collection import Setting
from batoms.base.object import ObjectGN, BaseObject


class Species(BaseObject):
    """_summary_

    Args:
        BaseObject (_type_): _description_
    """

    def __init__(self, name, parent,
                 data=None,
                 ) -> None:
        self.name = name
        self.parent = parent
        self.label = parent.label
        self.coll_name = "%s_instancer" % (self.label)
        BaseObject.__init__(self, obj_name="%s_%s" % (self.label, name),
                            bobj_name="binstancer")
        if data is not None:
            self.update(data)

    def update(self, data):
        """Update species by new data

        Args:
            data (dict): properites of species
        """
        data['elements'] = self.check_occupancy(data['elements'])
        #
        sp = self.data
        for key, value in data.items():
            if key == 'elements':
                continue
            setattr(sp, key, value)
        # add elements
        sp.elements.clear()
        for name, eledata in data['elements'].items():
            ele = sp.elements.add()
            ele.name = name
            ele.occupancy = eledata['occupancy']
            ele.color = eledata['color']
        self.build_instancer()

    @staticmethod
    def check_occupancy(elements):
        """Check input elements

        Args:
            elements (str, dict): information of elements

        Raises:
            ValueError: Total occupancy should be 1

        Returns:
            _type_: dict of element with its occupancy
        """
        if isinstance(elements, str):
            elements = {elements: {
                "occupancy": 1.0, "color": (0, 0.2, 0.8, 1)}}
        elif isinstance(elements, dict):
            occupancy = [ele['occupancy'] for ele in elements.values()]
            total = sum(occupancy)
            # not fully occupied.
            if total < 1 - 1e-6:
                elements['X'] = {"occupancy": 1 - total,
                                 "color": (0.8, 0.8, 0.0, 1.0)}
            elif total > 1 + 1e-6:
                raise ValueError("Total occupancy should be smaller than 1!")
        return elements

    def build_instancer(self):
        """Build object instancer for species

        Returns:
            bpy Object: instancer
        """
        sp = self.data
        name = '%s_%s' % (self.label, sp.name)
        radius = sp.radius*sp.scale
        segments = sp.segments
        shape = 'UV_SPHERE'
        self.delete_obj(name)
        if shape.upper() == 'UV_SPHERE':
            bpy.ops.mesh.primitive_uv_sphere_add(segments=segments[0],
                                                 ring_count=segments[1],
                                                 radius=radius)
        if shape.upper() == 'METABALL':
            bpy.ops.object.metaball_add(type='BALL', location=[0, 0, 0])
        obj = bpy.context.view_layer.objects.active
        obj.name = name
        obj.data.name = name
        obj.batoms.atom.radius = radius
        obj.batoms.type = 'INSTANCER'
        #
        obj.users_collection[0].objects.unlink(obj)
        bpy.data.collections['%s_instancer' % self.label].objects.link(obj)
        bpy.ops.object.shade_smooth()
        if shape.upper() != 'METABALL':
            obj.hide_set(True)
            obj.hide_render = True
        #
        self.build_materials()
        self.assign_materials()
        self.color = sp.color
        bpy.context.view_layer.update()
        self.parent.batoms.add_geometry_node(sp.name, obj)
        return obj

    def build_materials(self, node_inputs=None, material_style='default'):
        """Create materials for instancer object

        Args:
            node_inputs (dict, optional):
                properties of materials. Defaults to None.
            material_style (str, optional):
                Materials style. Defaults to 'default'.
        """
        from batoms.material import create_material
        mesh = self.obj.data
        mesh.materials.clear()
        elements = self.elements
        for occ in self.sorted_occupancies:
            ele = elements[occ[0]]
            name = '%s_%s_%s' % (self.label, self.name, ele["name"])
            self.delete_material(name)
            mat = create_material(name,
                                  color=ele["color"],
                                  node_inputs=node_inputs,
                                  material_style=material_style,
                                  backface_culling=True)
            mesh.materials.append(mat)

    def assign_materials(self):
        """Assign materials for instancer with order
        """
        # find the face index for ele
        mesh = self.obj.data
        occs = self.sorted_occupancies
        nele = len(occs)
        # for occupancy
        if nele > 1:
            # calc the angles
            npoly = len(mesh.polygons)
            normals = np.zeros(npoly*3)
            material_indexs = np.zeros(npoly, dtype='int')
            mesh.polygons.foreach_get('normal', normals)
            mesh.polygons.foreach_get('material_index', material_indexs)
            normals = normals.reshape(-1, 3)
            xy = normals - np.dot(normals, [0, 0, 1])[:, None]*[0, 0, 1]
            angles = (np.arctan2(xy[:, 1], xy[:, 0]) + np.pi)/np.pi/2
            #
            tos = 0
            for i in range(1, nele):
                toe = tos + occs[i][1]
                index = np.where((angles > tos) & (angles < toe))[0]
                material_indexs[index] = i
                tos = toe
            mesh.polygons.foreach_set('material_index', material_indexs)

    @property
    def data(self):
        """Get data for a this species

        Returns:
            Bspecies: Class Bspecies
        """
        data = self.parent.find(self.name)
        return data

    @property
    def occupancies(self):
        """Get occupancies for all elements

        Returns:
            dict: Occupancies for all elements for this species
        """
        occupancies = {}
        data = self.data.elements
        for eledata in data:
            occupancies[eledata.name] = round(eledata.occupancy, 3)
        return occupancies

    @property
    def sorted_occupancies(self):
        """Sort occupancies

        Returns:
            tuple: Sorted occupancies
        """
        # sort
        occupancies = self.occupancies
        sorted_occupancies = sorted(occupancies.items(), key=lambda x: -x[1])
        return sorted_occupancies

    @property
    def elements(self):
        """Elements of this species

        Returns:
            dict: dict of Belement
        """
        elements = {}
        for ele in self.data.elements:
            elements[ele.name] = ele.as_dict()
        return elements

    @elements.setter
    def elements(self, elements):
        """Set elements

        Args:
            elements (dict): elements with occupancis

        >>> H2O['O'].elements = 'O'
        >>> H2O['O'].elements = {'O':1.0}
        >>> H2O['O'].elements = {'O':0.7, 'S': 0.3}
        """
        data = self.data.as_dict()
        data['elements'] = elements
        self.update(data)

    @property
    def main_element(self):
        """Main element

        Returns:
            str: name of the element with largest occupancies,
            but vacancy "X" is not included.
        """
        main_element = {}
        occs = self.sorted_occupancies
        if occs[0][0] == 'X' and len(occs) > 1:
            main_element = occs[1][0]
        else:
            main_element = occs[0][0]
        return main_element

    @property
    def materials(self):
        """Get materials

        Returns:
            dict: dict of bpy material object
        """
        materials = {}
        for ele in self.elements:
            name = '%s_%s_%s' % (self.label, self.name, ele)
            mat = bpy.data.materials.get(name)
            materials[ele] = mat
        return materials

    @materials.setter
    def materials(self, materials):
        """Set materials properties of the instancer

        Args:
            materials (dict): parameters for materials

        >>> h2o['H'].materials = {'Metallic': 0.9,
        ...     'Specular': 1.0, 'Roughness': 0.01}
        """
        self.build_materials(node_inputs=materials)
        self.assign_materials()

    @property
    def color(self):
        return self.get_color()

    @color.setter
    def color(self, color):
        """
        >>> h.color = [0.8, 0.1, 0.3, 1.0]
        """
        self.set_color(color)

    def get_color(self):
        """
        """
        # Viewpoint_color = self.materials[self.main_element].diffuse_color
        for node in self.materials[self.main_element].node_tree.nodes:
            if 'Base Color' in node.inputs:
                node_color = node.inputs['Base Color'].default_value[:]
            if 'Alpha' in node.inputs:
                Alpha = node.inputs['Alpha'].default_value
        color = [node_color[0], node_color[1], node_color[2], Alpha]
        return color

    def set_color(self, color):
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.materials[self.main_element].diffuse_color = color
        for node in self.materials[self.main_element].node_tree.nodes:
            if 'Base Color' in node.inputs:
                node.inputs['Base Color'].default_value = color
            if 'Alpha' in node.inputs:
                node.inputs['Alpha'].default_value = color[3]

    # @property
    # def radius_style(self):
    #     return self.get_radius_style()

    # @radius_style.setter
    # def radius_style(self, radius_style):
    #     self.set_radius_style(radius_style)

    # def get_radius_style(self):
    #     return self.data.radius_style

    # def set_radius_style(self, radius_style):
    #     # print(species_props)
    #     radius_style = str(radius_style)
    #     self.data.radius_style = radius_style
    #     data = self.data.as_dict()
    #     data['radius_style'] = radius_style
    #     self.update(data)

    # @property
    # def color_style(self):
    #     return self.get_color_style()

    # @color_style.setter
    # def color_style(self, color_style):
    #     self.set_color_style(color_style)

    # def get_color_style(self):
    #     return self.data.color_style

    # def set_color_style(self, color_style):
    #     # print(species_props)
    #     color_style = str(color_style)
    #     self.data.color_style = color_style
    #     data = self.data.as_dict()
    #     data['color_style'] = color_style
    #     self.update(data)

    @property
    def indices(self):
        indices = np.where(
            self.parent.batoms.attributes['species'] == self.name)[0]
        return indices

    @property
    def scale(self):
        return self.parent.batoms.get_attribute_with_indices('scale',
                                                             self.indices)

    @scale.setter
    def scale(self, scale):
        self.parent.batoms.set_attribute_with_indices(
            'scale', self.indices, scale)

    @property
    def segments(self):
        return self.get_segments()

    @segments.setter
    def segments(self, segments):
        self.set_segments(segments)

    def get_segments(self):
        nverts = len(self.instancer.data.vertices)
        return nverts

    def set_segments(self, segments):
        if not isinstance(segments, int):
            raise Exception('Segments should be int!')
        self._species.update(segments=segments)
        self.build_geometry_node()

    def as_dict(self) -> dict:
        data = {
            'species': self.name,
            'name': self.name,
            'elements': self.elements,
        }
        return data

    def __repr__(self):
        s = "Bspecies(species = '%s', " % self.name
        s += "elements = %s, " % str(self.sorted_occupancies)
        s += "color = %s)" % np.round(np.array(self.color), 2)
        # s += "radius_style = %s, " % self.radius_style
        # s += "color_style = %s)" % self.color_style
        return s


class Bspecies(Setting):
    """Bspecies Class
    """

    def __init__(self, label, coll_name, species=None,
                 batoms=None,
                 material_style='default',
                 segments=None,
                 subdivisions=2,
                 ) -> None:
        Setting.__init__(self, label, coll_name=coll_name)
        self.label = label
        self.name = 'bspecies'
        self.batoms = batoms
        #
        self.calc_segments(segments)
        self.subdivisions = subdivisions
        if species is not None:
            self.add(species)
        # for sp, data in species.items():
            # self[sp] = data

    def calc_segments(self, segments):
        if segments is not None:
            self.segments = segments
            return
        if self.batoms is None:
            self.segments = [32, 24]
            return
        natom = len(self.batoms)
        if natom <= 1e3:
            segments = [32, 24]
        elif natom <= 1e4:
            segments = [16, 16]
        elif natom <= 1e4:
            segments = [10, 10]
        elif natom <= 1e5:
            segments = [8, 8]
        else:
            segments = [6, 6]
        self.segments = segments

    def __getitem__(self, name):
        if isinstance(name, str):
            sp = Species(name, parent=self)
        else:
            raise ValueError("Please input a string")
        return sp

    def __setitem__(self, name, data):
        """
        Set species
        """
        self.add({name: data})

    def add(self, props, instancer=None):
        """
        """
        if props is None:
            return
        if isinstance(props, str):
            props = [props]
        if isinstance(props, (list, np.ndarray)):
            species_list = props
            props = {}
            for species in species_list:
                ele = species.split('_')[0]
                props[species] = {"elements": {ele: {"occupancy": 1.0,
                                                     "color": (0.8, 0.8, 0, 1)}}}
        for name, data in props.items():
            # if 'color_style' not in data:
            # data['color_style'] = '0'
            # if 'radius_style' not in data:
            # data['radius_style'] = '0'
            sp = self.find(name)
            if sp is None:
                sp = self.collection.add()
            sp.name = name
            sp.label = self.label
            sp.segments = self.segments
            sp = Species(sp.name, parent=self, data=data)

    def update_geometry_node(self):
        instancers = self.instancers
        for sp, obj in instancers.items():
            self.batoms.add_geometry_node(sp, obj)

    def keys(self):
        return self.species.keys()

    def items(self):
        return self.species.items()

    def __repr__(self):
        s = ''
        for name, sp in self.species.items():
            s += str(sp)
            s += '\n'
        return s

    @ property
    def instancers(self):
        return self.get_instancers()

    def get_instancers(self):
        instancers = {}
        for sp in self.species:
            name = '%s_%s' % (self.label, sp)
            instancers[sp] = bpy.data.objects.get(name)
        return instancers

    @ property
    def species(self):
        return self.get_species()

    def get_species(self):
        species = {}
        for sp in self.collection:
            species[sp.name] = Species(sp.name, parent=self)
        return species

    @ property
    def main_elements(self):
        main_elements = {}
        for name, sp in self.species.items():
            main_elements[name] = sp.main_element
        return main_elements

    @ property
    def species_props(self):
        return self.get_species_props()

    def get_species_props(self):
        species_props = {}
        instancers = self.instancers
        species_props = {}
        for sp in self.species:
            radius = instancers[sp].batoms.atom.radius
            matname = instancers[sp].material_slots[0].name
            mat = bpy.data.materials.get(matname)
            color = mat.node_tree.nodes[0].inputs[0].default_value[:]
            species_props[sp] = {'radius': radius, 'color': color}
        return species_props

    def __add__(self, other):
        self += other
        return self

    def __iadd__(self, other):
        self.extend(other)
        return self

    def copy(self, label):
        """Copy this Bspecies class

        Args:
            name (str): new label

        Returns:
            Class Bspecies: copy
        """
        for sp, instancer in self.instancers.items():
            instancer = instancer.copy()
            bpy.data.collections['Collection'].objects.link(instancer)
            instancer.hide_set(True)
            instancer.name = '%s_%s_instancer' % (label, sp)
        bspecies = self.__class__(label, label, species={},
                                  batoms=None)
        return bspecies

    def extend(self, other):
        """

        """
        species1 = self.species
        species2 = other.species
        species3 = {}
        for key, sp in species2.items():
            data = sp.as_dict()
            if key not in species1:
                instancer = other.instancers[key]
                name = '%s_%s' % (self.label, key)
                instancer.name = name
                # materials
                # for mat in instancer.materials_slots:
                # mat.name.replace()
                species3[key] = [data, instancer]
            elif data != species1[key].as_dict():
                instancer = other.instancers[key]
                spname = '%s_%s' % (key, other.label)
                name = '%s_%s' % (self.label, spname)
                instancer.name = name
                species3[spname] = [data, instancer]
        #
        for key, data in species3.items():
            self.add({key: data[0]}, instancer=data[1])
        return self

    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]

    def update(self, segments=None):
        if segments:
            self.segments = segments
        for name, sp in self.species.items():
            sp.build_instancer()
