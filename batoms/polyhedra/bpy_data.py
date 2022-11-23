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



class PolyhedraSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    species: StringProperty(name="species")
    name: StringProperty(name="name")
    color: FloatVectorProperty(name="color",
                               subtype='COLOR',
                               min=0, max=1,
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
            'material_style': self.material_style,
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


class Polyhedra(Base):
    label: StringProperty(name="label", default='batoms')
    # polyhedra_style
    ui_list_index: IntProperty(name="ui_list_index",
                            default=0)
    # collection
    settings: CollectionProperty(name='polyhedrasetting',
                              type=PolyhedraSetting)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
