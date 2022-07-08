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



class MolecularSurfaceSetting(Base):
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
                               min=0, max=1,
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
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   select     probe   resolution    color  \n'
        s += '{:6s}  {:6s}  {:1.3f}  {:1.3f}  [{:1.2f}  {:1.2f}  {:1.2f}   {:1.2f}] \n'.format(
            self.name, self.select, self.probe, self.resolution, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s


class MolecularSurface(bpy.types.PropertyGroup):
    """This module defines the MolecularSurface properties to extend 
    Blender’s internal data.

    """
    active: BoolProperty(name="active", default=False)
    settings: CollectionProperty(name='MolecularSurfaceSetting',
                                type=MolecularSurfaceSetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)

    def as_dict(self) -> dict:
        setdict = {
            
        }
        return setdict