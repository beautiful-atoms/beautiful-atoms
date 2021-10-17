import json
import bpy
from bpy.types import Operator, AddonPreferences
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.props import (
        StringProperty,
        BoolProperty,
        EnumProperty,
        IntProperty,
        FloatProperty,
        )
from ase.io.cube import read_cube_data
from batoms.batoms import Batoms

class IMPORT_OT_batoms(Operator, ImportHelper):
    bl_idname = "import_mesh.batoms"
    bl_label  = "Import batoms (*.batoms)"
    bl_options = {"PRESET", "REGISTER", "UNDO"}

    filename_ext = ".batoms"

    camera: BoolProperty(
        name="Camera", default=False,
        description="Do you need a camera?")
    light: BoolProperty(
        name="Light", default=False,
        description = "Do you need a light?")
    world: BoolProperty(
        name="World", default=False,
        description = "Do you need a world light?")
    model_type: EnumProperty(
        name="Type",
        description="Choose model",
        items=(('0',"Space-filling", "Use ball"),
               ('1',"Ball-and-stick", "Use ball and stick"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
               default='0',)
    label: StringProperty(
        name = "label", 
        description = "Label")
    show_unit_cell: BoolProperty(
        name = "Show_unit_cell", default=False,
        description = "Show unit cell")

    def draw(self, context):
        
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.label(text="Adding Structure")
        box = layout.box()
        row = box.row()
        row.prop(self, "label")
        #
        row = layout.row()
        row.prop(self, "camera")
        row.prop(self, "light")
        row.prop(self, "world")
        #
        box = layout.box()
        row = box.row()
        row.prop(self, "show_unit_cell")
        # Balls
        box = layout.box()
        row = box.row()
        row.label(text="Structure model")
        row = box.row()
        col = row.column()
        col.prop(self, "model_type")
        box = layout.box()
        row = box.row()
        row.active = (self.model_type == "1")
        #
        

    def execute(self, context):

        # This is to determine the path.
        self.inputfile = bpy.path.abspath(self.filepath)

        # Execute main routine
        import_batoms(self.inputfile, 
               self.model_type,
               self.camera,
               self.light,
               self.world,
               self.show_unit_cell,
               self.label,
               )

        return {'FINISHED'}


def import_batoms(inputfile, 
               model_type = '0',
               camera = 'True',
               light = 'True',
               world = 'False',
               show_unit_cell = False,
               label = None,
               ):
    #
    from ase.io import read
    if not label:
        label = inputfile
    if inputfile.split('.')[-1] == 'cube':
        images, volume = read_cube_data(inputfile)
        batoms = Batoms(label = label, atoms = images, model_type=model_type, volume=volume)
    else:
        images = read(inputfile)
        batoms = Batoms(label = label, atoms = images, model_type=model_type)
    
