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



class MagresSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='batoms')
    type: EnumProperty(
        name="type",
        description="Surface type",
        items=(('MS', "MS", ""),
               ('CS', "CS", ""),
               ),
        default='MS')
    scale: FloatProperty(name="scale", soft_min=0.01, soft_max=1, default=0.01)
    resolution: FloatProperty(name="resolution", soft_min=0.2, soft_max=2,
                              default=0.5)
    select: StringProperty(name="select", default='all')
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 1, 1, 0.5])

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'color': self.color[:],
            'scale': self.scale,
            'resolution': self.resolution,
            'select': self.select,
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   select     scale   resolution    color  \n'
        s += '{:6s}  {:6s}  {:1.3f}  {:1.3f}  [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.select, self.scale, self.resolution, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s



class Magres(bpy.types.PropertyGroup):
    """This module defines the Magres properties to extend
    Blender’s internal data.

    """
    active: BoolProperty(name="active", default=False)
    settings: CollectionProperty(name='MagresSetting',
                                type=MagresSetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
