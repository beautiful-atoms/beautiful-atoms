import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

class Batoms(bpy.types.PropertyGroup):
    """
    """
    is_batoms: BoolProperty(name="is_batoms", default=False)
    model_type: EnumProperty(
        name="model_type",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball"),
               ('1',"Ball-and-stick", "Use ball and stick"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Wireframe", "Use wireframe")),
               default='0')
    polyhedra_type: EnumProperty(
        name="polyhedra_type",
        description="Polhhedra models",
        items=(('0',"atoms, bonds and polyhedra", "atoms, bonds and polyhedra"),
               ('1',"atoms, polyhedra", "atoms, polyhedra"),
               ('2',"central atoms, polyhedra","central atoms, polyhedra"),
               ('3',"polyhedra", "polyhedra")),
               default='0')
    pbc: BoolVectorProperty(name="pbc", default = [False, False, False], size = 3)
    cell: FloatVectorProperty(name="cell", default = [0, 0, 0, 0, 0, 0, 0, 0, 0], size = 9)
    show_unit_cell: BoolProperty(name="show_unit_cell", default = True)
    boundary: FloatVectorProperty(name="boundary", default = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0], size = 6)

class Batom(bpy.types.PropertyGroup):
    """
    """
    is_batom: BoolProperty(name="is_batom", default=False)
    label: StringProperty(name="species", default = 'X')
    species: StringProperty(name="species", default = 'X')
    element: StringProperty(name="element", default = '')
    radius: FloatProperty(name="radius")
class Bcell(bpy.types.PropertyGroup):
    """
    """
    is_bcell: BoolProperty(name="is_bcell", default=False)
    label: StringProperty(name="label", default = '')

class BBond(bpy.types.PropertyGroup):
    """
    """
    symbol1: StringProperty(name="symbol1")
    symbol2: StringProperty(name="symbol2")
    name:StringProperty(name = "name")
    min: FloatProperty(name="min", description = "min", default = 0.0)
    max: FloatProperty(name="max", description = "max", default = 2.0)
    search: IntProperty(name="search", default=0)
    polyhedra: BoolProperty(name="polyhedra", default=False)
    color1: FloatVectorProperty(name="color1", size = 4)
    color2: FloatVectorProperty(name="color1", size = 4)
    width: FloatProperty(name="width", default = 0.10)
    style: EnumProperty(
        name="style",
        description="bond style",
        items=(('0',"Unicolor cylinder", ""),
               ('1',"Bicolor cylinder", ""),
               ('2',"Dashed line", ""),
               ('3',"Dotted line", "")),
               default='1')
    def as_list(self) -> list:
        return [self.min, self.max, self.search, self.polyhedra, self.color1, self.color2, self.width, self.style]
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair      min     max   Search_bond    Polyhedra \n'
        s += '{0:10s} {1:4.3f}   {2:4.3f}      {3:10s}   {4:10s} \n'.format(\
                self.name, self.min, self.max, str(self.search), str(self.polyhedra))
        s += '-'*60 + '\n'
        return s

class BPolyhedra(bpy.types.PropertyGroup):
    """
    """
    symbol: StringProperty(name="symbol")
    name:StringProperty(name = "name")
    color: FloatVectorProperty(name="color", size = 4)
    edgewidth: FloatProperty(name="edgewidth", default = 0.10)
    def as_list(self) -> list:
        return [self.color, self.edgewidth]
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Center                color           edgewidth \n'
        s += '{0:10s}   [{1:1.2f}  {2:1.2f}  {3:1.2f}  {4:1.2f}]   {5:1.3f} \n'.format(\
                self.symbol, self.color[0], self.color[1], self.color[2], self.color[3], self.edgewidth)
        s += '-'*60 + '\n'
        return s


class BIsosurface(bpy.types.PropertyGroup):
    """
    """
    name: StringProperty(name = "name")
    npoint: IntProperty(name="npoint")
    level: FloatProperty(name="level", default = 0.10)
    color: FloatVectorProperty(name="color", size = 4, default = [1, 1, 0, 0.5])

    def as_list(self) -> list:
        return [self.level, self.color]
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name        level        color            \n'
        s += '{0:10s}   {1:1.6f}  [{2:1.2f}  {3:1.2f}  {4:1.2f}   {5:1.2f}] \n'.format(\
                self.name, self.level, self.color[0], self.color[1], self.color[2], self.color[3])
        s += '-'*60 + '\n'
        return s

class BVolume(bpy.types.PropertyGroup):
    """
    """
    is_bvolume: BoolProperty(name="is_bvolume", default=False)
    npoint: IntProperty(name="npoint")
    shape: IntVectorProperty(name="shape", size = 3)