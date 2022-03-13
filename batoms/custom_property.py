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


class Bspecies(bpy.types.PropertyGroup):
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
    show_unit_cell: BoolProperty(name="show_unit_cell", default=True)


class Bbond(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='batoms')
    # select: StringProperty(name="select", default='all')
    species1: StringProperty(name="species1")
    species2: StringProperty(name="species2")
    species: StringProperty(name="species")
    min: FloatProperty(name="min", description="min", default=0.0)
    max: FloatProperty(name="max", description="max", default=3.0)
    search: IntProperty(name="search", default=0)
    polyhedra: BoolProperty(name="polyhedra", default=False)
    color1: FloatVectorProperty(
        name="color1", size=4,
        subtype='COLOR',
        default=(0, 0.2, 0.8, 1))
    color2: FloatVectorProperty(
        name="color2", size=4,
        subtype='COLOR',
        default=(0.6, 0.2, 0, 1))
    width: FloatProperty(name="width", default=0.10)
    order: IntProperty(name="order", default=1)
    segments: IntProperty(name="segments", default=16)
    order_offset: FloatProperty(name="order_offset", default=0.1)
    style: EnumProperty(
        name="style",
        description="bond style",
        items=(('0', "Unicolor cylinder", ""),
               ('1', "Bicolor cylinder", ""),
               ('2', "Dashed line", ""),
               ('3', "Spring", "Spring")),
        default='1')

    @property
    def name(self) -> str:
        return '%s-%s-%s' % (self.species1, self.species2)

    def as_dict(self, reversed=False) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'species1': self.species1,
            'species2': self.species2,
            'species': self.species,
            'name': self.name,
            'min': self.min,
            'max': self.max,
            'search': self.search,
            'polyhedra': self.polyhedra,
            'color1': self.color1[:],
            'color2': self.color2[:],
            'width': self.width,
            'order': self.order,
            'order_offset': self.order_offset,
            'style': self.style
        }
        if reversed:
            setdict['name'] = '%s-%s' % (self.species2, self.species1)
            setdict['species1'] = self.species2
            setdict['species2'] = self.species1
            setdict['color1'] = self.color2[:]
            setdict['color2'] = self.color1[:]
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair   select   min     max   Search_bond    Polyhedra \n'
        s += '{:10s} {:10s} {:4.3f}   {:4.3f}      {:10s}   {:10s} \n'.format(
            self.name, self.select, self.min, self.max, str(self.search), str(self.polyhedra))
        s += '-'*60 + '\n'
        return s


class Bpolyhedra(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    species: StringProperty(name="species")
    name: StringProperty(name="name")
    color: FloatVectorProperty(name="color",
                               subtype='COLOR',
                               size=4)
    width: FloatProperty(name="width", min=0, soft_max=1, default=0.01)
    show_edge: BoolProperty(name="show_edge", default=True)

    @property
    def name(self) -> str:
        return self.species

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'species': self.species,
            'name': self.name,
            'color': self.color[:],
            'width': self.width,
            'show_edge': self.show_edge,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Center                show_edge           width \n'
        s += '{0:10s}    {1:10s}   {2:1.3f} \n'.format(
            self.species, str(self.show_edge), self.width)
        s += '-'*60 + '\n'
        return s


class Bisosurface(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
    level: FloatProperty(name="level", default=0.10)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               default=[1, 1, 0, 0.5],
                               )

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'color': self.color[:],
            'level': self.level,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name        level        color            \n'
        s += '{0:10s}   {1:1.6f}  [{2:1.2f}  {3:1.2f}  {4:1.2f}   {5:1.2f}] \n'.format(
            self.name, self.level, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s


class Bplane(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    indices: IntVectorProperty(
        name="Miller indices", size=3, default=[1, 0, 0])
    distance: FloatProperty(name="distance",
                            description="Distance from origin",
                            default=1)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               default=[0, 0, 1, 0.5]
                               )
    crystal: BoolProperty(name="crystal", default=False)
    symmetry: BoolProperty(name="symmetry", default=False)
    slicing: BoolProperty(name="slicing", default=False)
    boundary: BoolProperty(name="boundary", default=False)
    scale: FloatProperty(name="scale", default=1)
    show_edge: BoolProperty(name="show_edge", default=True)
    width: FloatProperty(name="width", default=0.01)

    @property
    def name(self) -> str:
        return '%s-%s-%s' % (self.indices[0], self.indices[1], self.indices[2])

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'color': self.color[:],
            'indices': self.indices,
            'distance': self.distance,
            'crystal': self.crystal,
            'symmetry': self.symmetry,
            'slicing': self.slicing,
            'boundary': self.boundary,
            'show_edge': self.show_edge,
            'width': self.width,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name        distance  crystal symmetry slicing  show_edge  boundary   edgewidth        \n'
        s += '{0:10s}   {1:1.3f}  {2:8s}  {3:8s}  {4:8s} {5:8s} {6:8s} {7:1.3f}\n'.format(
            self.name, self.distance, str(self.crystal), str(self.symmetry),
            str(self.slicing), str(self.show_edge), str(self.boundary), self.width)
        s += '-'*60 + '\n'
        return s


class Bvolume(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    npoint: IntProperty(name="npoint")
    shape: IntVectorProperty(name="shape", size=3)


class Blight(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    name: StringProperty(name="name", default='X')
    lock_to_camera: BoolProperty(name="lock_to_camera", default=False)
    direction: FloatVectorProperty(name="direction", default=[0, 0, 1], size=3)
    look_at: FloatVectorProperty(name="look_at", default=[0, 0, 0], size=3)


class Bcamera(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    name: StringProperty(name="name", default='X')
    direction: FloatVectorProperty(name="direction", default=[0, 0, 1], size=3)
    look_at: FloatVectorProperty(name="look_at", default=[0, 0, 0], size=3)


class Brender(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='X')
    engine: StringProperty(name="engine", default='BLENDER_EEVEE')
    compute_device_type: StringProperty(
        name="compute_device_type", default='CUDA')
    animation: BoolProperty(name="animation", default=False)
    run_render: BoolProperty(name="run_render", default=True)
    gpu: BoolProperty(name="gpu", default=False)
    viewport: FloatVectorProperty(name="viewport", default=[0, 0, 1], size=3)
    center: FloatVectorProperty(name="center", default=[0, 0, 1], size=3)
    distance: FloatProperty(name="distance",
                            description="Distance from camera",
                            default=-1)
    padding: FloatVectorProperty(name="padding", default=[1, 1, 1, 1], size=4)


class Bsheet(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    name: StringProperty(name="name")
    # sheetId: IntProperty(name="sheetId", default=0)
    # chainId: StringProperty(name="chainId", default='A')
    startChain: StringProperty(name="startChain")
    startResi: IntProperty(name="startResi", default=0)
    endChain: StringProperty(name="endChain")
    endResi: IntProperty(name="endResi", default=0)
    color: FloatVectorProperty(
        name="color1", size=4, default=(0.0, 0.0, 1.0, 1.0))
    extrude: FloatProperty(name="extrude", default=0.8)
    depth: FloatProperty(name="depth", default=0.1)

    @property
    def name(self) -> str:
        return '%s-%s-%s-%s' % (self.startChain, self.startResi, self.endChain, self.endResi)

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'startChain': self.startChain,
            'startResi': self.startResi,
            'endChain': self.endChain,
            'endResi': self.endResi,
            'name': self.name,
            'color': self.color[:],
            'extrude': self.extrude,
            'depth': self.depth,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s


class Bhelix(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    name: StringProperty(name="name")
    # helixId: IntProperty(name="helixId", default=0)
    # chainId: StringProperty(name="chainId", default='A')
    startChain: StringProperty(name="startChain")
    startResi: IntProperty(name="startResi", default=0)
    endChain: StringProperty(name="endChain")
    endResi: IntProperty(name="endResi", default=0)
    color: FloatVectorProperty(
        name="color1", size=4, default=(1.0, 0.0, 0.0, 1))
    extrude: FloatProperty(name="extrude", default=1.0)
    depth: FloatProperty(name="depth", default=0.1)

    @property
    def name(self) -> str:
        return '%s-%s-%s-%s' % (self.startChain, self.startResi, self.endChain, self.endResi)

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'startChain': self.startChain,
            'startResi': self.startResi,
            'endChain': self.endChain,
            'endResi': self.endResi,
            'name': self.name,
            'color': self.color[:],
            'extrude': self.extrude,
            'depth': self.depth,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s


class Bturn(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    name: StringProperty(name="name")
    # turnId: IntProperty(name="turnId", default=0)
    # chainId: StringProperty(name="chainId", default='A')
    startChain: StringProperty(name="startChain")
    startResi: IntProperty(name="startResi", default=0)
    endChain: StringProperty(name="endChain")
    endResi: IntProperty(name="endResi", default=0)
    color: FloatVectorProperty(
        name="color1", size=4, default=(0.0, 1.0, 0.0, 1))
    radius: FloatProperty(name="radius", default=0.2)

    @property
    def name(self) -> str:
        return '%s-%s-%s-%s' % (self.startChain, self.startResi, self.endChain, self.endResi)

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'startChain': self.startChain,
            'startResi': self.startResi,
            'endChain': self.endChain,
            'endResi': self.endResi,
            'name': self.name,
            'color': self.color[:],
            'radius': self.radius,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s


class Bms(bpy.types.PropertyGroup):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Surface type",
        items=(('SAS', "SAS", "Solvent Accessible Surface"),
               ('SES', "SES", "solvent-excluded surfaces"),
               ),
        default='SAS')
    probe: FloatProperty(name="probe", soft_min=0.4, soft_max=2, default=1.4)
    resolution: FloatProperty(name="resolution", soft_min=0.2, soft_max=2,
                              default=0.5)
    select: StringProperty(name="select", default='all')
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               default=[0, 1, 1, 1.0])

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'color': self.color[:],
            'probe': self.probe,
            'resolution': self.resolution,
            'select': self.select,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   select     probe   resolution    color  \n'
        s += '{:6s}  {:6s}  {:1.3f}  {:1.3f}  [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.select, self.probe, self.resolution, self.color[0], self.color[1], self.color[2], self.color[3])
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


class BatomsCollection(bpy.types.PropertyGroup):
    """
    """
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Structural models",
        items=(('OTHER', "Not Batoms", "Not Batoms"),
               ('BATOMS', "Batoms", "Batoms"),
               ('VOLUME', "Volume", "Volume"),
               ('BOND', "Bond", "Bond"),
               ('POLYHEDRA', "Polyhedra", "Polyhedra"),
               ('CRYSTALSHAPE', "Crystal Shape", "Crystal Shape"),
               ('LATTICEPLANE', "LatticePlane", "LatticePlane"),
               ('ISOSURFACE', "Isosurface", "Isosurface"),
               ('MS', "MS", "Molecular Surface"),
               ),
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

    show_unit_cell: BoolProperty(name="show_unit_cell", default=True)

    show_search: BoolProperty(name="show_search", default=False)

    wrap: BoolVectorProperty(name="wrap", default=[
                             False, False, False], size=3)

    boundary: FloatVectorProperty(name="boundary", default=[
                                  0.0, 1.0, 0.0, 1.0, 0.0, 1.0], size=6)

    # collection
    bbond: CollectionProperty(name='Bbond',
                              type=Bbond)

    bond_index: IntProperty(name="bond_index",
                            default=0)

    blatticeplane: CollectionProperty(name='Blatticeplane',
                                      type=Bplane)

    latticeplane_index: IntProperty(name="latticeplane_index",
                                    default=0)

    bcrystalshape: CollectionProperty(name='Bcrystalshape',
                                      type=Bplane)

    crystalshape_index: IntProperty(name="crystalshape_index",
                                    default=0)

    bisosurface: CollectionProperty(name='Bisosurface',
                                    type=Bisosurface)

    isosurface_index: IntProperty(name="isosurface_index",
                                  default=1)

    bpolyhedra: CollectionProperty(name='Bpolyhedra',
                                   type=Bpolyhedra)

    polyhedra_index: IntProperty(name="polyhedra_index",
                                 default=0)

    bsheet: CollectionProperty(name='Bsheet',
                               type=Bsheet)

    bhelix: CollectionProperty(name='Bhelix',
                               type=Bhelix)

    bturn: CollectionProperty(name='Bturn',
                              type=Bturn)

    bselect: CollectionProperty(name='Bselect',
                                type=Bselect)

    brender: PointerProperty(name='Brender',
                             type=Brender)

    bspecies: CollectionProperty(name='Bspecies',
                                 type=Bspecies)

    species_index: IntProperty(name="species_index",
                               default=1)

    bms: CollectionProperty(name='Bms',
                            type=Bms)

    ms_index: IntProperty(name="ms_index",
                          default=1)


class BatomsObject(bpy.types.PropertyGroup):
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Structural models",
        items=(('OTHER', "Not Batoms", "Not Batoms"),
               ('BATOMS', "Batoms", "Batoms"),
               ('VOLUME', "Volume", "Volume"),
               ('CELL', "Cell", "Cell"),
               ('BOND', "Bond", "Bond"),
               ('INSTANCER', "Instancer", "Instancer"),
               ('POLYHEDRA', "Polyhedra", "Polyhedra"),
               ('CRYSTALSHAPE', "Crystal Shape", "Crystal Shape"),
               ('LATTICEPLANE', "LatticePlane", "LatticePlane"),
               ('ISOSURFACE', "Isosurface", "Isosurface"),
               ('MS', "MS", "Molecular Surface"),
               ),
        default='OTHER')
    # obj
    atom: PointerProperty(name='Batom',
                          type=Batom)
    cell: PointerProperty(name='Bcell',
                          type=Bcell)
    bond: PointerProperty(name='Bbond',
                          type=Bbond)
    plane: PointerProperty(name='Bplane',
                           type=Bplane)
    volume: PointerProperty(name='Bvolume',
                            type=Bvolume)
    isosurface: PointerProperty(name='Bisosurface',
                                type=Bisosurface)
    polyhedra: PointerProperty(name='Bpolyhedra',
                               type=Bpolyhedra)
    sheet: PointerProperty(name='Bsheet',
                           type=Bsheet)
    light: PointerProperty(name='Blight',
                           type=Blight)
    camera: PointerProperty(name='Bcamera',
                            type=Bcamera)
