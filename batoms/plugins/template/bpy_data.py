import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       CollectionProperty,
                       )

from batoms.internal_data import Base



class TemplateSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
    prop1: FloatProperty(name="prop1", default=1)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 0, 1, 0.5]
                               )

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'name': self.name,
            'label': self.label,
            'color': self.color[:],
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
    active: BoolProperty(name="active", default=False)
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0', "Surface", "Surface"),
               ('1', "Dot", "Dot surface"),
               ('2', "Wireframe", "Use wireframe")),
        default='0')
    show: BoolProperty(name="show", default=True)
    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)
    settings: CollectionProperty(name='TemplateSetting',
                                type=TemplateSetting)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
