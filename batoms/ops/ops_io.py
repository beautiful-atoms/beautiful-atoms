import bpy
from bpy.types import Operator
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.props import (
    BoolProperty,
    StringProperty,
    EnumProperty,
)
from batoms.batoms import Batoms
from batoms.bio.bio import read


class IMPORT_OT_batoms(Operator, ImportHelper):
    bl_idname = "batoms.import"
    bl_label = "Import"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".xyz"

    label: StringProperty(
        name="label", default='',
        description="Name of batoms.")


    render: BoolProperty(
        name="Add a default render", default=True,
        description="Do you need a render?")

    model_style: EnumProperty(
        name="Type",
        description="Choose model",
        items=(('0', "Space-filling", "Use ball"),
               ('1', "Ball-and-stick", "Use ball and stick"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Stick", "Use stick")),
        default='0',)

    def draw(self, context):

        layout = self.layout
        layout.label(text="Import a Structure")
        box = layout.box()
        col = box.column()
        col.prop(self, "label")
        col.prop(self, "render")
        #
        box = layout.box()
        row = box.row()
        row.label(text="Structure model")
        row = box.row()
        col = row.column()
        col.prop(self, "model_style")

    def execute(self, context):
        inputfile = bpy.path.abspath(self.filepath)
        if self.label == '':
            batoms = read(inputfile)
        else:
            batoms = read(inputfile, label = self.label)
        batoms.model_style = self.model_style
        return {'FINISHED'}


class EXPORT_OT_batoms(Operator, ExportHelper):
    bl_idname = "batoms.export"
    bl_label = "Export"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".xyz"

    def draw(self, context):

        layout = self.layout
        layout.label(text="Export a Structure")

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.write(self.filepath)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


def menu_func_import_batoms(self, context):
    lay = self.layout
    lay.operator(IMPORT_OT_batoms.bl_idname,
                 text="Batoms file (xyz, cif, pdb, ...)")


def menu_func_export_batoms(self, context):
    lay = self.layout
    lay.operator(EXPORT_OT_batoms.bl_idname,
                 text="Batoms file (xyz, cif, pdb, ...)")
