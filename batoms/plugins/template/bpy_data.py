import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       CollectionProperty,
                       )

from batoms.custom_property import Base


class TemplateSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
   

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'name': self.name,
            'label': self.label,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name      label\n'
        s += '{:10s}    {:10s} \n'.format(self.name, self.label)
        s += '-'*60 + '\n'
        return s

class Template(bpy.types.PropertyGroup):
    """This module defines the Template properties to extend 
    Blender’s internal data.

    """
    show: BoolProperty(name="show", default=True)
    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)
    setting: CollectionProperty(name='TemplateSetting',
                                type=TemplateSetting)
