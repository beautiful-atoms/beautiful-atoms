import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       FloatVectorProperty,
                       CollectionProperty,
                       )

from batoms.internal_data import Base


class BondSetting(Base):
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
        min=0, max=1,
        default=(0, 0.2, 0.8, 1))
    color2: FloatVectorProperty(
        name="color2", size=4,
        subtype='COLOR',
        min=0, max=1,
        default=(0.6, 0.2, 0, 1))
    width: FloatProperty(name="width", default=0.10)
    scale: FloatProperty(name="scale", default=1.0)
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
    type: IntProperty(name="type", default=0)

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
            'material_style': self.material_style,
            'color1': self.color1[:],
            'color2': self.color2[:],
            'width': self.width,
            'order': self.order,
            'order_offset': self.order_offset,
            'style': self.style,
            'type': self.type,
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
        s = 'Bondpair   min     max   Search_bond    Polyhedra \n'
        s += '{:10s} {:4.3f}   {:4.3f}      {:10s}   {:10s} \n'.format(
            self.name, self.min, self.max, str(self.search), str(self.polyhedra))
        s += '-'*60 + '\n'
        return s


class Bond(Base):
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='batoms')
    show_search: BoolProperty(name="show_search", default=False)
    show_hydrogen_bond: BoolProperty(name="show_hydrogen_bond", default=False)
    ui_list_index: IntProperty(name="ui_list_index",
                            default=0)
    # collection
    settings: CollectionProperty(name='bondsetting',
                              type=BondSetting)
