import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       CollectionProperty,
                       )

from batoms.internal_data import Base


class SheetSetting(Base):
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
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s


class HelixSetting(Base):
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
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s


class TurnSetting(Base):
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
            'material_style': self.material_style,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   startChain   startResi   endChain   endResi\n'
        s += '{0:10s} {1:10s}   {2:10s}      {3:10s}   {4:10s} \n'.format(
            self.name, self.startChain, str(self.startResi), self.endChain, str(self.endResi))
        s += '-'*60 + '\n'
        return s

class Protein(bpy.types.PropertyGroup):
    """This module defines the Protein properties to extend
    Blender’s internal data.

    """
    turn_ui_list_index: IntProperty(name="turn_ui_list_index",
                              default=0)
    sheet_ui_list_index: IntProperty(name="sheet_ui_list_index",
                              default=0)
    helix_ui_list_index: IntProperty(name="helix_ui_list_index",
                              default=0)

    # colleciton properties
    settings_turn: CollectionProperty(name='settings_turn',
                                type=TurnSetting)
    settings_sheet: CollectionProperty(name='settings_sheet',
                                type=SheetSetting)
    settings_helix: CollectionProperty(name='settings_helix',
                                type=HelixSetting)
