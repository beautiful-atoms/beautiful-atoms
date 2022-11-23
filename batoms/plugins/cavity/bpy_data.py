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



class CavitySetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    name: StringProperty(name="name")
    label: StringProperty(name="label", default='')
    min: FloatProperty(name="min", default=5)
    max: FloatProperty(name="max", default=6)
    scale: FloatProperty(name="scale", default=1.0)
    resolution: FloatProperty(name="resolution", soft_min=0.5, soft_max=2,
                              default=1.0)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[1, 1, 0, 1.0],
                               )

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'min': self.min,
            'max': self.max,
            'scale': self.scale,
            'name': self.name,
            'material_style': self.material_style,
            'color': self.color[:],
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name      min    max        color            \n'
        s += '{:10s}   {:1.2f} {:1.2f} [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.min, self.max, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s

class Cavity(bpy.types.PropertyGroup):
    """This module defines the cavity properties to extend
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
    atomRadius: FloatProperty(name="atomRadius", default=0.5)
    minCave: FloatProperty(name="minCave", default=3)
    resolution: FloatProperty(name="resolution", default=1)
    settings: CollectionProperty(name='CavitySetting',
                                type=CavitySetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
