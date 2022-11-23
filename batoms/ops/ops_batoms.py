import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (
    StringProperty,
    BoolProperty,
    FloatVectorProperty,
    FloatProperty,
    IntProperty,
    IntVectorProperty,
    EnumProperty,
)
from batoms import Batoms
from batoms.ops.base import OperatorBatoms, OperatorBatomsEdit
from batoms.utils.butils import get_selected_vertices


class ApplyCell(OperatorBatoms):
    bl_idname = "batoms.apply_cell"
    bl_label = "Apply Cell"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new cell parameters")

    a: FloatVectorProperty(
        name="a", default=(1, 0, 0),
        # subtype = "XYZ",
        description="Cell in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1, 0),
        # subtype = "XYZ",
        description="Cell in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 0, 1),
        # subtype = "XYZ",
        description="Cell in c axis")

    def execute(self, context):
        cell = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms.cell = cell
        batoms.obj.select_set(True)
        return {'FINISHED'}


class ApplyTransform(OperatorBatoms):
    bl_idname = "batoms.apply_transform"
    bl_label = "Apply Transform"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new transform parameters")

    a: IntVectorProperty(
        name="a", default=(1, 0, 0, 10), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "MATRIX",
        description="Transform matrix in a axis")
    b: IntVectorProperty(
        name="b", default=(0, 1, 0, 10), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "XYZ",
        description="Transform matrix in b axis")
    c: IntVectorProperty(
        name="c", default=(0, 0, 1, 0), size=4,
        soft_min=-10, soft_max=10,
        # subtype = "XYZ",
        description="Transform matrix in c axis")

    def execute(self, context):
        transform = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms1 = batoms.transform(transform)
        # batoms.obj.select_set(True)
        # batoms1.obj.select_set(False)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class ApplyBoundary(OperatorBatoms):
    bl_idname = "batoms.apply_boundary"
    bl_label = "Apply Boundary"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new boundary parameters")

    a: FloatVectorProperty(
        name="a", default=(0, 1), size=2,
        description="boundary  in a axis")
    b: FloatVectorProperty(
        name="b", default=(0, 1), size=2,
        description="boundary  in b axis")
    c: FloatVectorProperty(
        name="c", default=(0, 1), size=2,
        description="boundary  in c axis")

    def execute(self, context):
        boundary = [self.a, self.b, self.c]
        batoms = Batoms(label=context.object.batoms.label)
        batoms.boundary = boundary
        batoms.obj.select_set(True)
        return {'FINISHED'}


class AddSurface(OperatorBatoms):
    bl_idname = "batoms.surface_add"
    bl_label = "Add surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add surface")

    label: StringProperty(
        name="Label", default='surf',
        description="Label")

    indices: IntVectorProperty(
        name="Miller indices", size=3, default=(1, 1, 1),
        soft_min=-10, soft_max=10,
        description="Miller indices for the plane")

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 1),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")

    nlayer: IntProperty(name="layers", min=1, soft_max=10, default=4)

    termination: StringProperty(
        name="termination", default='',
        description="termination")

    def execute(self, context):
        from ase.build import surface
        from ase.build.surfaces_with_termination import surfaces_with_termination
        batoms = Batoms(label=context.object.batoms.label)
        bulk = batoms.as_ase()
        if len(self.termination) == 0:
            atoms = surface(bulk, self.indices,
                            self.nlayer, self.vacuum)
        else:
            atoms = surfaces_with_termination(bulk, self.indices,
                                              self.nlayer, self.vacuum,
                                              termination=self.termination)
        batoms.hide = True
        label = batoms.label + ''.join(str(i) for i in self.indices)
        if self.label in bpy.data.collections:
            self.label = "%s.001" % self.label
        batoms = Batoms(label=label, from_ase=atoms, movie=True)
        batoms = batoms*self.size
        batoms.translate([2, 2, 2])
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class AddRootSurface(OperatorBatoms):
    bl_idname = "batoms.root_surface_add"
    bl_label = "Add root surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add root surface")

    label: StringProperty(
        name="Label", default='rootsurf',
        description="Label")

    root: IntProperty(name="root", min=1, soft_max=49, default=3)

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 1),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    eps: FloatProperty(name="eps", min=1e-8, soft_max=1e-3, default=1e-6)

    def execute(self, context):
        from ase.build import root_surface
        batoms = Batoms(label=context.object.batoms.label)
        surf = batoms.as_ase()
        atoms = root_surface(surf, self.root, eps = self.eps)
        batoms.hide = True
        label = batoms.label + '_root'
        atoms = atoms*self.size
        if self.label in bpy.data.collections:
            self.label = "%s.001" % self.label
        batoms = Batoms(label=label, from_ase=atoms, movie=True)
        batoms.translate([2, 2, 2])
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class BatomsJoin(OperatorBatoms):
    bl_idname = "batoms.join"
    bl_label = "Join batoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Join batoms")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        from batoms.utils.butils import get_selected_batoms
        batoms_list = get_selected_batoms()
        if len(batoms_list) < 2:
            self.report({"WARNING"}, "Please select at least two Batoms objects.")
            return {'FINISHED'}
        if self.label == '':
            self.label = batoms_list[0]
        if self.label not in batoms_list:
            self.report({"INFO"}, "Create a new Batoms object: {}.".format(self.label))
        batoms = Batoms(label=self.label)
        for label in batoms_list:
            if self.label == label: continue
            batoms += Batoms(label=label)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

class BatomsSeparate(OperatorBatoms):
    bl_idname = "batoms.separate"
    bl_label = "Separate batoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Separate batoms based on selects")

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.separate()
        # batoms.obj.select_set(True)
        self.report({"INFO"}, "Separate the Batoms object: {}.".format(self.label))
        return {'FINISHED'}

class deleteBatoms(Operator):
    bl_idname = "batoms.delete"
    bl_label = "Delete batoms"
    bl_description = ("Delete batoms")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        from batoms.utils.butils import remove_collection, read_batoms_list
        if self.label != '':
            coll = bpy.data.collections.get(self.label)
            remove_collection(self.label, keep_batom=False)
        else:
            batoms_list = read_batoms_list()
            for label in batoms_list:
                coll = bpy.data.collections.get(label)
                remove_collection(label, keep_batom=False)
                self.report({"INFO"}, "Delete {}".format(label))
        return {'FINISHED'}

class deleteSelectedBatoms(OperatorBatoms):
    bl_idname = "batoms.delete_selected_batoms"
    bl_label = "Delete selected batoms"
    bl_description = ("Delete selected batoms")

    def execute(self, context):
        from batoms.utils.butils import get_selected_batoms, remove_collection
        batoms_list = get_selected_batoms()
        for label in batoms_list:
            coll = bpy.data.collections.get(label)
            remove_collection(label, keep_batom=False)
            self.report({"INFO"}, "Delete {}".format(label))
        return {'FINISHED'}

class ApplyModelStyle(OperatorBatoms):
    bl_idname = "batoms.apply_model_style"
    bl_label = "Apply Model Style"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new model_style parameters")

    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=[('0', "Space-filling", "Use ball and stick"),
               ('1', "Ball-and-stick", "Use ball"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Stick", "Use stick")
            ],
        default='0',
    )

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.model_style = int(self.model_style)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Model style of {} is changed to {}"
            .format(context.object.batoms.label, self.model_style))
        return {'FINISHED'}

class ApplyRadiusStyle(OperatorBatoms):
    bl_idname = "batoms.apply_radius_style"
    bl_label = "Apply Radius Style"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new radius_style parameters")

    radius_style: EnumProperty(
        name="radius_style",
        description="Structural models",
        items=[('0', "Covalent", "Covalent"),
               ('1', "VDW", "van der Waals"),
               ('2', "Ionic", "Ionic")
            ],
        default='0',
    )

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.radius_style = int(self.radius_style)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Radius style of {} is changed to {}"
            .format(context.object.batoms.label, self.radius_style))
        return {'FINISHED'}


class ApplyColorStyle(OperatorBatoms):
    bl_idname = "batoms.apply_color_style"
    bl_label = "Apply Color Style"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new color_style parameters")

    color_style: EnumProperty(
        name="color_style",
        description="Structural models",
        items=[('0', "JMOL", "JMOL"),
               ('1', "VESTA", "VESTA"),
               ('2', "CPK", "CPK")
            ],
        default='0',
    )

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.color_style = int(self.color_style)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Color style of {} is changed to {}"
            .format(context.object.batoms.label, self.color_style))
        return {'FINISHED'}


class BatomsCopy(OperatorBatoms):
    bl_idname = "batoms.copy"
    bl_label = "Copy"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Copy a batoms")

    label: StringProperty(
        name="label", default='',
        description="The label of the copy.")


    def execute(self, context):
        if self.label == '':
            return {'FINISHED'}
        batoms = Batoms(label=context.object.batoms.label)
        batoms1 = batoms.copy(self.label)
        # batoms.obj.select_set(True)
        # batoms1.obj.select_set(False)
        bpy.context.view_layer.objects.active = batoms1.obj
        self.report({"INFO"}, "Copy {} to {}"
            .format(batoms.label, self.label))
        return {'FINISHED'}


class ApplyLabel(OperatorBatoms):
    bl_idname = "batoms.apply_label"
    bl_label = "Apply Label"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new label parameters")

    label: StringProperty(
        name="label", default='',
        description="Show label")

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.show_label = self.label
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Show label {} for {}"
            .format(self.label, context.object.batoms.label))
        return {'FINISHED'}

class BatomsAutoBuildSpecies(OperatorBatoms):
    bl_idname = "batoms.auto_build_species"
    bl_label = "Auto build species"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Auto build species")

    tol: FloatProperty(
        name="tol", default=1e-5,
        description="tol")

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        batoms.auto_build_species(self.tol)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Auto build species for {}"
            .format(context.object.batoms.label))
        return {'FINISHED'}
# ===========================================
# Below for Edit mode
# ===========================================


class BatomsReplace(OperatorBatomsEdit):
    bl_idname = "batoms.replace"
    bl_label = "Replace"
    bl_description = "Replace selected atoms by new species"

    species: StringProperty(
        name="species", default='',
        description="Replaced by this species")

    def execute(self, context):
        if self.species == '':
            return {'FINISHED'}
        obj = context.object
        v = get_selected_vertices(obj)
        batoms = Batoms(label=obj.batoms.label)
        batoms.replace(v, self.species)
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class BatomModify(OperatorBatomsEdit):
    bl_idname = "batoms.batom_modify"
    bl_label = "Modify batom"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify batom")

    key: StringProperty(
        name="key", default='scale',
        description="Replaced by this species")

    scale: FloatProperty(name="scale", default=0.6,
                         min=0, soft_max=2, description="scale",
                         )
    bond: BoolProperty(name="Bond", default=True,
                       )

    show: BoolProperty(name="Show", default=True,
                       )

    def execute(self, context):
        obj = context.object
        v = get_selected_vertices(obj)
        batoms = Batoms(label=obj.batoms.label)
        setattr(batoms[v], self.key, getattr(self, self.key))
        context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode="EDIT")
        return {'FINISHED'}


class ApplyModelStyleSelected(OperatorBatoms):
    bl_idname = "batoms.apply_model_style_selected"
    bl_label = "Apply Model Style Selected"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Apply new model_style parameters to selected atoms")

    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=[('0', "Space-filling", "Use ball and stick"),
               ('1', "Ball-and-stick", "Use ball"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Stick", "Use stick")
            ],
        default='0',
    )

    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        indices = get_selected_vertices(context.object)
        model_style_array = batoms.get_attribute("model_style")
        model_style_array[indices] = int(self.model_style)
        scale_array = batoms.get_attribute("scale")
        if int(self.model_style) == 0:
            scale = 1
        elif int(self.model_style) == 1:
            scale = 0.4
        elif int(self.model_style) == 2:
            scale = 0.4
        elif int(self.model_style) == 3:
            scale = 0.0001
        scale_array[indices] = scale
        batoms.set_model_style_array(model_style_array)
        batoms.set_attributes({'scale': scale_array})
        batoms.obj.select_set(True)
        self.report({"INFO"}, "Model style of {} atoms are changed to {}".format(len(indices), self.model_style))
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

class BatomsCopySelected(OperatorBatoms):
    bl_idname = "batoms.copy_selected"
    bl_label = "Copy selected atoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Copy selected atoms")

    select: StringProperty(
        name="select", default='',
        description="New select for the copied atoms.")

    label: StringProperty(
        name="label", default='',
        description="Label for the new batoms.")


    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        indices = get_selected_vertices(context.object)
        batoms.obj.select_set(True)
        self.report({"INFO"}, "Copy {} atoms".format(len(indices)))
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class BatomsDeleteSelectedAtoms(OperatorBatoms):
    """Delete selected vertices

    Args:
        OperatorBatoms (_type_): _description_

    Returns:
        _type_: _description_
    """
    bl_idname = "batoms.delete_selected_atoms"
    bl_label = "Delete selected atoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Delete selected atoms")


    def execute(self, context):
        batoms = Batoms(label=context.object.batoms.label)
        indices = get_selected_vertices(context.object)
        del batoms[indices]
        self.report({"INFO"}, "Delete {} atoms".format(len(indices)))
        bpy.context.view_layer.objects.active = batoms.obj
        bpy.ops.object.mode_set(mode="EDIT")
        return {'FINISHED'}
