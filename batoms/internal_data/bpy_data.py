import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       CollectionProperty,
                       )

def get_volumetric_data(self, context):
    keys = bpy.data.collections[context.object.batoms.label].batoms.settings_volumetric_data.keys()
    items = []
    i = 0
    for key in keys:
        items.append((key, key, key))
        i += 1
    return items

class Base(bpy.types.PropertyGroup):

    material_style: EnumProperty(
        name="material_style",
        description="Structural models",
        items=(('default', "default", "default"),
               ('metallic', "metallic", "metallic"),
               ('plastic', "plastic", "plastic"),
               ('ceramic', "ceramic", "ceramic"),
               ('mirror', "mirror", "mirror")),
        default='default')


class Belement(bpy.types.PropertyGroup):
    name: StringProperty(name="name", default='H')
    occupancy: FloatProperty(name="occupancy")
    radius: FloatProperty(name="radius", default=1)
    color: FloatVectorProperty(name="color", size=4, default=(0, 0.2, 0.8, 1))

    def as_dict(self) -> dict:
        setdict = {
            'occupancy': self.occupancy,
            'name': self.name,
            # 'radius': self.radius,
            'color': self.color[:],
        }
        return setdict


class Bspecies(Base):
    name: StringProperty(name="name", default='H')
    species: StringProperty(name="species", default='H')
    elements: CollectionProperty(name='Belements',
                                 type=Belement)
    scale: FloatProperty(name="scale", min=0, soft_max=2, default=1)
    radius: FloatProperty(name="radius", default=1.0)
    segments: IntVectorProperty(name="segments", size=2, default=(24, 16))
    color: FloatVectorProperty(
        name="color", size=4,
        subtype='COLOR',
        min=0, max=1,
        default=(0, 0.2, 0.8, 1))

    @property
    def element_dict(self):
        element_dict = {}
        for ele in self.elements:
            element_dict[ele.name] = ele.as_dict()
        return element_dict

    def as_dict(self) -> dict:
        setdict = {
            'species': self.species,
            'name': self.name,
            'material_style': self.material_style,
            'color': self.color[:],
            'scale': self.scale,
            'elements': self.element_dict,
        }
        return setdict


class Batom(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    # select: StringProperty(name="select", default = 'all')
    species: CollectionProperty(name='Bspecies',
                                type=Bspecies)
    radius: FloatProperty(name="radius")
    radii_style: EnumProperty(
        name="radii_style",
        description="Radii",
        items=(('0', "Covalent", "covalent"),
               ('1', "VDW", "van der Waals"),
               ('2', "Ionic", "ionic")),
        default='0')


class Battribute(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='batoms')
    data_type: EnumProperty(
        name="data_type",
        description="data type",
        items=(('FLOAT', "FLOAT", "Floating-point value"),
               ('INT', "INT", "32-bit integer"),
               ("FLOAT_VECTOR", "FLOAT_VECTOR", "Vector – 3D vector with floating-point values."),
               ("FLOAT_COLOR", "FLOAT_COLOR", "Color – RGBA color with floating-point values."),
               ("BYTE_COLOR", "BYTE_COLOR", "Byte Color – RGBA color with 8-bit values."),
               ("STRING", "STRING", "String – Text string."),
               ("BOOLEAN", "BOOLEAN", "Boolean – True or false."),
               ("FLOAT2", "FLOAT2", "2D Vector – 2D vector with floating-point values."),
               ),
        default='FLOAT')
    domain: EnumProperty(
        name="domain",
        description="Domain",
        items=(('POINT', "POINT", ""),
               ('EDGE', "EDGE", ""),
               ("FACE", "FACE", ""),
               ),
        default='POINT')
    dimension: IntProperty(name="index", default=1)
    shape_: IntVectorProperty(name="shape_", soft_min=0, size=32)
    delimiter: StringProperty(name="delimiter", default='@')

    @property
    def natt(self) -> int:
        import numpy as np
        return np.product(self.shape)

    @property
    def sub_name(self) -> int:
        sub_name = ["{}{}".format(self.name, i) for i in range(self.natt)]
        return sub_name

    @property
    def shape(self) -> int:
        import numpy as np
        return self.shape_[:self.dimension]

    @shape.setter
    def shape(self, data):
        tmp = list(self.shape_)
        tmp[0:self.dimension] = list(data)
        self.shape_ = tmp

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'data_type': self.data_type,
            'domain': self.domain,
            'dimension': self.dimension,
            'shape': self.shape,
        }
        return setdict

    def __repr__(self) -> str:
        s = "-"*60 + "\n"
        s = "{:20s}{:10s}{:10s}{:10s}{:10s}   {:20s}\n".format("Name", "Type", "Domain", "delimiter", "Dimension", "Shape")
        s += "{:20s}{:10s}{:10s}{:10s}{:10d}  [".format(
            self.name, self.data_type, self.domain, self.delimiter, self.dimension)
        for i in range(self.dimension):
            s += "  {}  ".format(self.shape[i])
        s += "] \n"
        s += "-"*80 + "\n"
        return s

class Bsite(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    index: IntProperty(name="index", default=0)


class Bcell(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    pbc: BoolVectorProperty(name="pbc", default=[False, False, False], size=3)
    width: FloatProperty(name="width", default=0.05, min=0)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0.2, 0.2, 0.2, 1])

    def as_dict(self) -> dict:
        setdict = {
            'pbc': self.pbc,
            'width': self.width,
            'color': self.color[:],
        }
        return setdict

class Bvolumetric_data(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    npoint: IntProperty(name="npoint")
    shape: IntVectorProperty(name="shape", size=3)

    def as_dict(self) -> dict:
        setdict = {
            'label': self.label,
            'shape': self.shape[:],
        }
        return setdict

class Bboundary(Base):
    """
    """
    name: StringProperty(name="name")
    flag: BoolProperty(name="flag", default=False)
    active: BoolProperty(name="active", default=False)
    label: StringProperty(name="label", default='batoms')
    boundary: FloatVectorProperty(name="boundary", default=[
                                  0.0, 1.0, 0.0, 1.0, 0.0, 1.0], size=6)
    scale: FloatProperty(name="scale", soft_min=0.01, soft_max=1, default=0.01)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 1, 1, 0.5])

    def as_dict(self) -> dict:
        setdict = {
            'name': self.name,
            'flag': self.flag,
            'label': self.label,
            'active': self.active,
            'color': self.color[:],
            'boundary': self.boundary[:],
            'scale': self.scale,
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        import numpy as np
        s = '-'*60 + '\n'
        s = 'Name    scale   color  \n'
        s += '{:6s}  {:1.3f}  [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.scale, self.color[0], self.color[1], self.color[2], self.color[3])
        boundary = np.array(self.boundary)
        s += ' %s \n '.format(boundary)
        s += '-'*60 + '\n'
        return s


class Bselect(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name", default='all')
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0', "Space-filling", "Use ball"),
               ('1', "Ball-and-stick", "Use ball and stick"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Wireframe", "Use wireframe")),
        default='0')
    radius_style: EnumProperty(
        name="radius_style",
        description="Radii",
        items=(('0', "Covalent", "covalent"),
               ('1', "VDW", "van der Waals"),
               ('2', "Ionic", "ionic")),
        default='0')
    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(('0', "atoms, bonds and polyhedra", "atoms, bonds and polyhedra"),
               ('1', "atoms, polyhedra", "atoms, polyhedra"),
               ('2', "central atoms, polyhedra", "central atoms, polyhedra"),
               ('3', "polyhedra", "polyhedra")),
        default='0')
    bspecies: CollectionProperty(name='Bspecies',
                                 type=Bspecies)

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   model_style   radius_style     show \n'
        s += '{0:10s} {1:10s}   {2:10s}   {3:10s} \n'.format(
            self.name, str(self.model_style), self.radius_style, str(self.show))
        s += '-'*60 + '\n'
        return s


items = [('OTHER', "Not Batoms", "Not Batoms"),
               ('BATOMS', "Batoms", "Batoms"),
               ('VOLUME', "Volume", "Volume"),
               ('CELL', "Cell", "Cell"),
               ('BOND', "Bond", "Bond"),
               ('BOUNDARY', "Boundary", "Boundary"),
               ('INSTANCER', "Instancer", "Instancer"),
               ('POLYHEDRA', "Polyhedra", "Polyhedra"),
               ('CRYSTALSHAPE', "Crystal Shape", "Crystal Shape"),
               ('LATTICEPLANE', "LatticePlane", "LatticePlane"),
               ('ISOSURFACE', "Isosurface", "Isosurface"),
               ('MS', "MS", "Molecular Surface"),
               ('CAVITY', "CAVITY", "Cavity"),
               ('MARGES', "Marges", "NMR tensors"),
               ('HIGHLIGHT', "HIGHLIGHT", "Highlight"),
               ]

class BatomsCollection(bpy.types.PropertyGroup):
    """
    """
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Structural models",
        items=items,
        default='OTHER')

    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0', "Space-filling", "Use ball"),
               ('1', "Ball-and-stick", "Use ball and stick"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Wireframe", "Use wireframe")),
        default='0')

    radius_style: EnumProperty(
        name="radius_style",
        description="Radii",
        items=(('0', "Covalent", "covalent"),
               ('1', "VDW", "van der Waals"),
               ('2', "Ionic", "ionic")),
        default='0')

    color_style: EnumProperty(
        name="color_style",
        description="Color",
        items=(('0', "JMOL", "JMOL"),
               ('1', "VESTA", "VESTA"),
               ('2', "CPK", "CPK")),
        default='0')

    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(('0', "atoms, bonds and polyhedra", "atoms, bonds and polyhedra"),
               ('1', "atoms, polyhedra", "atoms, polyhedra"),
               ('2', "central atoms, polyhedra", "central atoms, polyhedra"),
               ('3', "polyhedra", "polyhedra")),
        default='0')

    segments: IntVectorProperty(name="segments", size=2, default=(24, 16))
    scale: FloatProperty(name="scale", default=1)
    show: BoolProperty(name="show", default=True)
    show_unit_cell: BoolProperty(name="show_unit_cell", default=True)
    show_axes: BoolProperty(name="show_axes", default=True)
    show_label: StringProperty(name="show_label", default="")
    wrap: BoolVectorProperty(name="wrap", default=[
                             False, False, False], size=3)
    boundary: PointerProperty(name="Bboundary", type=Bboundary)
    cell: PointerProperty(name='Bcell', type=Bcell)
    crystal_view: BoolProperty(name="crystal_view", default=False)
    ui_list_index_species: IntProperty(name="ui_list_index_species",
                               default=0)
    ui_list_index_select: IntProperty(name="ui_list_index_select",
                               default=0)
    ui_list_index_volumetric_data: IntProperty(name="ui_list_index_volumetric_data",
                               default=0)
    # collection
    settings_select: CollectionProperty(name='settings_select',
                                type=Bselect)

    settings_species: CollectionProperty(name='settings_species',
                                 type=Bspecies)
    settings_volumetric_data: CollectionProperty(name='settings_volumetric_data',
                                 type=Bvolumetric_data)

    def as_dict(self) -> dict:
        setdict = {
            'model_style': self.model_style,
            'color_style': self.color_style,
            'radius_style': self.radius_style,
            'polyhedra_style': self.polyhedra_style,
            'segments': self.segments,
            'wrap': self.wrap,
            'show_axes': self.show_axes,
            'show_label': self.show_label,
            'crystal_view': self.crystal_view,
            'scale': self.scale,
        }
        return setdict

class BatomsObject(bpy.types.PropertyGroup):
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Structural models",
        items=items,
        default='OTHER')
    # obj
    atom: PointerProperty(name='Batom',
                          type=Batom)
    cell: PointerProperty(name='Bcell',
                          type=Bcell)
    volume: PointerProperty(name='Bvolumetric_data',
                            type=Bvolumetric_data)

    # collection
    settings_attribute: CollectionProperty(name='settings_attribute',
                                 type=Battribute)
