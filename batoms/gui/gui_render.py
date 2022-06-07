import bpy
from bpy.types import Panel
from bpy.props import (
    FloatVectorProperty,
    IntVectorProperty,
    FloatProperty,
    EnumProperty,
)
from batoms.render.render import Render
from batoms import Batoms

class Render_PT_prepare(Panel):
    bl_label = "Render"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "RENDER_PT_Tools"

    def draw(self, context):
        layout = self.layout
        repanel = context.scene.repanel

        layout = self.layout
        # row = col.row(align=True)
        layout.operator("batoms.render_add")
        layout.label(text="Camera")
        row = layout.row()
        row.prop(repanel, "viewport_x")
        row.prop(repanel, "viewport_y")
        row.prop(repanel, "viewport_z")
        col = layout.column()
        col.label(text="Type")
        col.prop(repanel, "camera_type", expand=True)
        if repanel.camera_type == 'ORTHO':
            layout.prop(repanel, "scale")
        else:
            layout.prop(repanel, "camera_lens")
            layout.prop(repanel, "distance")
        # layout.prop(repanel, "resolution")
        layout.separator()
        layout.label(text="Light")
        row = layout.row()
        row.prop(repanel, "light_direction_x")
        row.prop(repanel, "light_direction_y")
        row.prop(repanel, "light_direction_z")
        layout.prop(repanel, "light_energy")

def get_viewport(i):
    """Helper function to easily get cell-property."""

    def getter(self):
        render = get_active_render_collection()
        if render is not None:
            return render.viewport[i]
        else:
            return 0
    return getter

def set_viewport(i):
    """Helper function to easily get cell-property."""

    def setter(self, value):
        render = get_active_render_collection()
        if render is not None:
            viewport = render.viewport
            viewport[i] = value
            render = get_active_render()
            render.viewport = viewport
    return setter



def get_camera_type(self):
    camera = get_default_camera()
    if camera is not None:
        if camera.data.type == 'ORTHO':
            return 0
        else:
            return 1
    else:
        return 0

def set_camera_type(self, value):
    items = self.bl_rna.properties["camera_type"].enum_items
    item = items[value]
    type = item.identifier
    # print(type)
    self["type"] = type
    set_camera_attr('type', type)

def get_lens(self):
    camera = get_default_camera()
    if camera is not None:
        return camera.data.lens
    else:
        return 50

def set_lens(self, value):
    self["lens"] = value
    set_camera_attr('lens', value)

def get_scale(self):
    camera = get_default_camera()
    if camera is not None:
        return camera.data.ortho_scale
    else:
        return 10

def set_scale(self, value):
    self["scale"] = value
    set_camera_attr('scale', value)

def get_distance(self):
    render = get_active_render_collection()
    if render is not None:
        return render.distance
    else:
        return 10

def set_distance(self, value):
    self["distance"] = value
    render = get_active_render()
    render.distance = value

def get_light_direction(i):
    """Helper function to easily get cell-property."""

    def getter(self):
        light = get_default_light()
        if light is not None:
            return light.batoms.light.direction[i]
        else:
            return 0
    return getter

def set_light_direction(i):
    """Helper function to easily set light-property."""

    def setter(self, value):
        light = get_default_light()
        if light is not None:
            light_direction = light.batoms.light.direction
            light_direction[i] = value
            set_light_attr('direction', light_direction)
    return setter

def get_energy(self):
    light = get_default_light()
    if light is not None:
        return light.data.energy
    else:
        return 10

def set_energy(self, value):
    self["energy"] = value
    set_light_attr('energy', value)

class RenderProperties(bpy.types.PropertyGroup):
    """
    """
    viewport_x: FloatProperty(
        name="viewport_x",default=0,
        soft_min = -5, soft_max = 5,
        description="Miller viewport for the render",
        get=get_viewport(0),
        set=set_viewport(0))
    
    viewport_y: FloatProperty(
        name="viewport_y",default=0,
        soft_min = -5, soft_max = 5,
        description="Miller viewport for the render",
        get=get_viewport(1),
        set=set_viewport(1))
    
    viewport_z: FloatProperty(
        name="viewport_z",default=1,
        soft_min = -5, soft_max = 5,
        description="Miller viewport for the render",
        get=get_viewport(2),
        set=set_viewport(2))
    
    camera_type: EnumProperty(
        name="type",
        description="camera type",
        items=(('ORTHO', "ORTHO", "", 0),
               ('PERSP', "PERSP", "", 1)
               ),
        default=0,
        get=get_camera_type,
        set=set_camera_type,
    )

    # resolution: IntVectorProperty(
        # name="resolution", size=2, default=(1000, 1000),
        # soft_min = 100, soft_max = 2000,
        # description="Miller resolution for the render", update=Callback_modify_resolution)
    
    scale: FloatProperty(
        name="scale", default=10,
        soft_min = 1, soft_max = 100,
        description="scale",
        get=get_scale,
        set=set_scale)

    light_direction_x: FloatProperty(
        name="light_direction_x",default=0,
        soft_min = -5, soft_max = 5,
        description="Light direction for the render",
        get=get_light_direction(0),
        set=set_light_direction(0))
    
    light_direction_y: FloatProperty(
        name="light_direction_y",default=0,
        soft_min = -5, soft_max = 5,
        description="Light direction for the render",
        get=get_light_direction(1),
        set=set_light_direction(1))
    
    light_direction_z: FloatProperty(
        name="light_direction_z",default=1,
        soft_min = -5, soft_max = 5,
        description="Light direction for the render",
        get=get_light_direction(2),
        set=set_light_direction(2))

    light_energy: FloatProperty(
        name="energy", default=10,
        soft_min = 1, soft_max = 100,
        description="light_energy",
        get=get_energy,
        set=set_energy)

    distance: FloatProperty(
        name="distance", default=3,
        description="distance from origin",
        get=get_distance,
        set=set_distance)


    camera_lens: FloatProperty(
        name="lens", default=100,
        soft_min = 1, soft_max = 100,
        get=get_lens,
        set=set_lens)



def modify_render_attr(context, key, value):
    from batoms.batoms import Batoms
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        setattr(batoms.render, key, value)
        if key not in ['distance']:
            batoms.render.init()
        context.space_data.region_3d.view_perspective = 'CAMERA'


def set_light_attr(key, value):
    from batoms.batoms import Batoms
    if bpy.context.object and bpy.context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=bpy.context.object.batoms.label)
        setattr(batoms.render.lights['Default'], key, value)
        batoms.render.init()
        bpy.context.space_data.region_3d.view_perspective = 'CAMERA'

def set_camera_attr(key, value):
    from batoms.batoms import Batoms
    if bpy.context.object and bpy.context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=bpy.context.object.batoms.label)
        if key == 'scale':
            batoms.render.camera.set_ortho_scale(value)
        elif key in ['lens']:
            setattr(batoms.render.camera, key, value)
        else:
            setattr(batoms.render.camera, key, value)
            batoms.render.init()
        bpy.context.space_data.region_3d.view_perspective = 'CAMERA'

def get_active_render_collection():
    """Get the collection of the active Batoms

    When get the attribute of Batoms object, 
    if the attribute if saved in the Batoms.coll.batoms,
    we only need to read data form the colleciton, 
    it is faster than get data from the Batoms itself.

    Returns:
        bpy.type.collection: _description_
    """
    if 'batoms_render' in bpy.data.collections:
        return bpy.data.collections['batoms_render'].batoms.brender
    return None

def get_active_render():
    context = bpy.context
    if 'batoms_render' in bpy.data.collections:
        render = Render(label='batoms')
        return render
    return None

def get_default_camera():
    if 'batoms_camera_Default' in bpy.data.objects:
        return bpy.data.objects['batoms_camera_Default']
    return None

def get_default_light():
    if 'batoms_light_Default' in bpy.data.objects:
        return bpy.data.objects['batoms_light_Default']
    return None