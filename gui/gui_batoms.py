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
from batoms.butils import get_selected_batoms, get_selected_vertices, \
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

        box = layout.box()
        col = box.column()
        col.label(text="Model type")
        col = box.column()
        col.prop(bapanel, "model_type", expand  = True)
        
        row = box.row()
        row.prop(bapanel, "hide", expand  = True)
        
        row = box.row()
        row.prop(bapanel, "scale")

        box = layout.box()
        col = box.column(align=True)
        col.operator("batoms.replace")
        row = box.row(align=True)
        row.prop(bapanel, "species")

        col = box.column(align=True)
        col.operator("batoms.fragmentate")
        row = box.row(align=True)
        row.prop(bapanel, "suffix", expand  = True)

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Measurement")
        col.operator("batoms.record_selection")
        col.prop(bapanel, "measurement")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Export atoms")
        col.prop(bapanel, "filetype")

class BatomsProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_batom(self):
        return get_selected_objects('batom')
    @property
    def selected_vertices(self):
        return get_selected_vertices()
    def Callback_model_type(self, context):
        bapanel = bpy.context.scene.bapanel
        model_type = list(bapanel.model_type)[0]
        modify_batoms_attr(self.selected_batoms, 'model_type', model_type)
    def Callback_modify_hide(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'hide', bapanel.hide)
    def Callback_modify_scale(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'scale', bapanel.scale)
    def Callback_export_atoms(self, context):
        bapanel = bpy.context.scene.bapanel
        export_atoms(self.selected_batoms, bapanel.filetype)
    
    model_type: EnumProperty(
        name="model_type",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball and stick"),
               ('1',"Ball-and-stick", "Use ball"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
        default={'0'}, 
        update=Callback_model_type,
        options={'ENUM_FLAG'},
        )
    
    hide: BoolProperty(name="hide",
                default=False, 
                description = "Hide all object for view and rendering",
                update = Callback_modify_hide)
    single: BoolProperty(name="single", 
                default=False, 
                description = "Separate the species into single atom")
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
        batoms = Batoms(label = name)
        setattr(batoms, key, value)
        batoms_list.append(batoms)
    # Restor the selected state
    for batoms in batoms_list:
        batoms.select = True
def export_atoms(selected_batoms, filetype = 'xyz'):
    batoms_list = []
    for name in selected_batoms:
        batoms = Batoms(label = name)
        batoms.write('%s.%s'%(batoms.label, filetype))
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True

def replace_atoms(species1):
    batom_list = []
    selected_vertices = get_selected_vertices()
    for label, species, name, index in selected_vertices:
        batoms = Batoms(label = label)
        batoms.replace(species, species1, index)
        batom_list.append(name)
    for name in batom_list:
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')

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
    for label, species, name, index in selected_vertices:
        batoms = Batoms(label = label)
        batoms.fragmentate(species, index, suffix)
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
