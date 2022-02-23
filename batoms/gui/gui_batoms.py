from ase.build import molecule, bulk
from ase import Atom, Atoms
import bpy
from bpy.types import (Panel,
                       Operator,
                       )
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
from batoms.utils.butils import get_selected_batoms, get_selected_vertices, \
                        get_selected_objects
from batoms import Batoms, Batom

# The panel.
class Batoms_PT_prepare(Panel):
    bl_label       = "Batoms"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    def draw(self, context):
        layout = self.layout
        bapanel = context.scene.bapanel

        layout.label(text="Model style")
        col = layout.column()
        col.prop(bapanel, "model_style", expand  = True)
        layout.label(text="Radius style")
        layout.prop(bapanel, "radius_style", expand  = True)
        
        layout.prop(bapanel, "show", expand  = True)
        layout.prop(bapanel, "scale")

        layout.operator("batoms.replace")
        layout.prop(bapanel, "species")

        layout.operator("batoms.fragmentate")
        layout.prop(bapanel, "suffix", expand  = True)

        layout.label(text="Measurement")
        layout.operator("batoms.record_selection")
        layout.prop(bapanel, "measurement")

        layout.label(text="Export atoms")
        layout.prop(bapanel, "filetype")

class BatomsProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_vertices(self):
        return get_selected_vertices()
    def Callback_model_style(self, context):
        bapanel = bpy.context.scene.bapanel
        model_style = list(bapanel.model_style)[0]
        modify_batoms_attr(self.selected_batoms, 'model_style', model_style)
    def Callback_radius_style(self, context):
        bapanel = bpy.context.scene.bapanel
        radius_style = list(bapanel.radius_style)[0]
        modify_batoms_attr(self.selected_batoms, 'radius_style', radius_style)
    def Callback_modify_show(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'show', bapanel.show)
    def Callback_modify_scale(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'scale', bapanel.scale)
    def Callback_export_atoms(self, context):
        bapanel = bpy.context.scene.bapanel
        export_atoms(self.selected_batoms, bapanel.filetype)
    
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball and stick"),
               ('1',"Ball-and-stick", "Use ball"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
        default={'0'}, 
        update=Callback_model_style,
        options={'ENUM_FLAG'},
        )
    radius_style: EnumProperty(
        name="radius_style",
        description="Structural models",
        items=(('0',"Covalent", "covalent"),
               ('1',"VDW", "van der Waals"),
               ('2',"Ionic","ionic")),
        default={'0'}, 
        update=Callback_radius_style,
        options={'ENUM_FLAG'},
        )
    
    show: BoolProperty(name="show",
                default=False, 
                description = "show all object for view and rendering",
                update = Callback_modify_show)
    scale: FloatProperty(
        name="scale", default=1.0,
        description = "scale", update = Callback_modify_scale)
    measurement: StringProperty(
        name="", default='',
        description = "measurement in Angstrom, degree")
    species: StringProperty(
        name="species", default='O_1',
        description = "Replaced by this species")
    suffix: StringProperty(
        name="suffix", default='f',
        description = "Replaced by this suffix")
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Label", default='h2o',
        description = "Label")
    filetype: StringProperty(
        name = "File type", default='xyz',
        description = "save batoms to file", update = Callback_export_atoms)

def modify_batoms_attr(selected_batoms, key, value):
    """
    """
    batoms_list = []
    for name in selected_batoms:
        batoms = Batoms(name)
        setattr(batoms, key, value)
        batoms_list.append(batoms)
def export_atoms(selected_batoms, filetype = 'xyz'):
    batoms_list = []
    for name in selected_batoms:
        batoms = Batoms(name)
        batoms.write('%s.%s'%(batoms.label, filetype))
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True

def replace_atoms(species):
    selected_vertices = get_selected_vertices()
    for label, index in selected_vertices:
        batoms = Batoms(label = label)
        batoms.replace(index, species)

class ReplaceButton(Operator):
    bl_idname = "batoms.replace"
    bl_label = "Replace by"
    bl_description = "Replace selected atoms by new species"
    
    @classmethod
    def poll(cls, context):
        return context.mode in  {'EDIT_MESH'}

    def execute(self, context):
        bapanel = context.scene.bapanel
        replace_atoms(bapanel.species)
        return {'FINISHED'}
class MeasureButton(Operator):
    bl_idname = "batoms.measure"
    bl_label = "Measure by"
    bl_description = "Measure selected atoms by new species"

    def execute(self, context):
        bapanel = context.scene.bapanel
        replace_atoms(bapanel.species)
        return {'FINISHED'}

def fragmentate(suffix):
    batom_list = []
    selected_vertices = get_selected_vertices()
    for label, indices in selected_vertices:
        batoms = Batoms(label = label)
        batoms.fragmentate(species, indices, suffix)
        batom_list.append(name)
    for name in batom_list:
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
class FragmentateButton(Operator):
    bl_idname = "batoms.fragmentate"
    bl_label = "Fragmentate"
    bl_description = "Fragmentate selected atoms"

    @classmethod
    def poll(cls, context):
        return context.mode in  {'EDIT_MESH'}

    def execute(self, context):
        bapanel = context.scene.bapanel
        fragmentate(bapanel.suffix)
        return {'FINISHED'}
