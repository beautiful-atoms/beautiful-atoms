import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       CollectionProperty,
                       )

from batoms.custom_property import Base


class IsosurfaceSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
    level: FloatProperty(name="level", default=0.10)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[1, 1, 0, 0.5],
                               )

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'material_style': self.material_style,
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


class Isosurface(bpy.types.PropertyGroup):
    """This module defines the Isosurface properties to extend 
    Blender’s internal data.

    """
    settings: CollectionProperty(name='IsosurfaceSetting',
                                type=IsosurfaceSetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)
