import bpy
from bpy.props import (
                       IntVectorProperty,
                       StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       FloatVectorProperty,
                       CollectionProperty,
                       )

from batoms.internal_data import Base



class LatticePlaneSetting(Base):
    """
    """
    flag: BoolProperty(name="flag", default=False)
    label: StringProperty(name="label", default='')
    indices: IntVectorProperty(
        name="Miller indices", size=3, default=[1, 0, 0])
    distance: FloatProperty(name="distance",
                            description="Distance from origin",
                            default=1)
    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 0, 1, 0.5]
                               )
    crystal: BoolProperty(name="crystal", default=False)
    symmetry: BoolProperty(name="symmetry", default=False)
    slicing: BoolProperty(name="slicing", default=False)
    boundary: BoolProperty(name="boundary", default=False)
    scale: FloatProperty(name="scale", default=1)
    show_edge: BoolProperty(name="show_edge", default=False)
    width: FloatProperty(name="width", default=0.01)

    color: FloatVectorProperty(name="color", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 1, 1, 1.0])
    color1: FloatVectorProperty(name="color1", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[1, 0, 0, 1.0])
    color2: FloatVectorProperty(name="color2", size=4,
                               subtype='COLOR',
                               min=0, max=1,
                               default=[0, 0, 1, 1.0])

    def get_color_by(self, context):
        keys = bpy.data.collections[self.label].batoms.settings_volumetric_data.keys()
        items = [("None", "None", "None"),
                ]
        i = 1
        for key in keys:
            items.append((key, key, key))
            i += 1
        return items

    color_by : EnumProperty(
        items=get_color_by,
        name="color_by",
        description="color by",
        default=None,
        update=None,
        )

    @property
    def name(self) -> str:
        return '%s-%s-%s' % (self.indices[0], self.indices[1], self.indices[2])

    def as_dict(self) -> dict:
        setdict = {
            'flag': self.flag,
            'label': self.label,
            'name': self.name,
            'material_style': self.material_style,
            'color_by': self.color_by,
            'color': self.color[:],
            'color1': self.color1[:],
            'color2': self.color2[:],
            'indices': list(self.indices),
            'distance': self.distance,
            'crystal': self.crystal,
            'symmetry': self.symmetry,
            'slicing': self.slicing,
            'boundary': self.boundary,
            'show_edge': self.show_edge,
            'width': self.width,
        }
        return setdict

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name        distance  crystal symmetry slicing  show_edge  boundary   edgewidth        \n'
        s += '{0:10s}   {1:1.3f}  {2:8s}  {3:8s}  {4:8s} {5:8s} {6:8s} {7:1.3f}\n'.format(
            self.name, self.distance, str(self.crystal), str(self.symmetry),
            str(self.slicing), str(self.show_edge), str(self.boundary), self.width)
        s += '-'*60 + '\n'
        return s


class LatticePlane(bpy.types.PropertyGroup):
    """This module defines the LatticePlane properties to extend
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
    settings: CollectionProperty(name='LatticePlaneSetting',
                                type=LatticePlaneSetting)

    ui_list_index: IntProperty(name="ui_list_index",
                              default=0)

    def as_dict(self) -> dict:
        setdict = {

        }
        return setdict
