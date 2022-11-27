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



class HighlightSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
    scale: FloatProperty(name="scale", default=1.0)
    select: StringProperty(name="select", default='all')
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[1, 1, 0, 0.4],
                               )
    style: EnumProperty(
        name="style",
        description="bond style",
        items=(('0', "Sphere", ""),
               ('1', "Cubic", "")),
        default='0')

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'scale': self.scale,
            'name': self.name,
            'select': self.select,
            'material_style': self.material_style,
            'color': self.color[:],
            'style': self.style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name      select    scale        color            \n'
        s += '{:10s}   {:1.2f} {:1.2f} [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.select, self.scale, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s

class Highlight(bpy.types.PropertyGroup):
    """This module defines the highlight properties to extend
    Blender’s internal data.

    """
    active: BoolProperty(name="active", default=False)
    show: BoolProperty(name="show", default=True)
    settings: CollectionProperty(name='HighlightSetting',
                                type=HighlightSetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
