"""
Rendering ase atoms objects using blender.

"""
import bpy
from mathutils import Vector, Matrix
import os
import numpy as np

default_render_settings = {
        'transparent': True,  # transparent background
        'resolution_x': 1000,
        'resolution_y': None,  # 
        'lock_camera_to_view': True,
        'lock_light_to_camera': True,
        'camera_loc': [0, 0, 100],  # x, y is the image plane, z is *out* of the screen
        'camera_target': [0, 0, 0], #
        'camera_type': 'ORTHO',  #  ['PERSP', 'ORTHO']
        'ortho_scale': None, #
        'camera_lens': 50,  #
        'fstop': 0.5,
        'ratio': 1,
        'light_loc': [0, 0, 100],
        'light_type': 'SUN', # 'POINT', 'SUN', 'SPOT', 'AREA'
        'light_energy': 5.0,
        'studiolight': 'Default',  # "basic", "outdoor", "paint", "rim", "studio"
        'world': False,
        'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
        'use_motion_blur': False,
        'motion_blur_position': 'START', 
        'motion_blur_steps': 10,
        'motion_blur_shutter': 30.0,
        'frame': 0,
        'functions': [],
        'run_render': True,
        'output': None,
        'animation': False,
        'save_to_blend': False,
        'queue': None,
        'gpu': True,
        'num_samples': 32,
        }

class Render():
    """
    Object to render atomic structure.

    """
    #
    def __init__(self, label, batoms = None, **kwargs):
        self.label = label
        self.batoms = batoms
        self.camera_name = 'camera_%s'%self.label
        self.light_name = 'light_%s'%self.label
        default_render_settings.update(kwargs)
        self.set_parameters(default_render_settings)
        self.clean_default()
        self.set_coll()
        
        if self.world:
            self.set_world()
    def set_parameters(self, parameters):
        for k, v in parameters.items():
            setattr(self, k, v)
    def clean_default(self, camera = False, light = True):
        if 'Cube' in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects["Cube"], do_unlink=True)
        if camera and 'Camera' in bpy.data.cameras:
            bpy.data.cameras.remove(bpy.data.cameras['Camera'])
        if light and 'Light' in bpy.data.lights:
            bpy.data.lights.remove(bpy.data.lights['Light'])
    def set_coll(self):
        name = "%s_render"%self.label
        if name not in bpy.data.collections:
            coll = bpy.data.collections.new(name)
            self.scene.collection.children.link(coll)
    @property
    def scene(self):
        return self.get_scene()
    def get_scene(self):
        return bpy.data.scenes['Scene']
    @property
    def coll(self):
        return self.get_coll()
    def get_coll(self):
        name = "%s_render"%self.label
        return bpy.data.collections[name]
    @property
    def engine(self):
        return self.get_engine()
    @engine.setter
    def engine(self, engine):
        self.set_engine(engine)
    def get_engine(self):
        return self.scene.render.engine
    def set_engine(self, engine):
        if engine.upper() == 'EEVEE':
            engine = 'BLENDER_EEVEE'
        elif engine.upper() == 'WORKBENCH':
            engine = 'BLENDER_WORKBENCH'
        self.scene.render.engine = engine.upper()
    @property
    def camera(self):
        return self.get_camera()
    @camera.setter
    def camera(self, camera):
        self.set_camera(camera)
    def get_camera(self):
        return bpy.data.objects[self.camera_name]
    @property
    def light(self):
        return self.get_light()
    @light.setter
    def light(self, light):
        self.set_light(light)
    def get_light(self):
        return bpy.data.objects[self.light_name]
    def set_world(self, color = [0.9, 0.9, 0.9, 1.0]):
        world = self.scene.world
        world.use_nodes = True
        node_tree = world.node_tree
        rgb_node = node_tree.nodes.new(type="ShaderNodeRGB")
        rgb_node.outputs["Color"].default_value = color
        node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
        node_tree.links.new(rgb_node.outputs["Color"], node_tree.nodes["Background"].inputs["Color"])
    def set_camera(self, camera_type = None, camera_lens = None,):
        '''

        camera_type: str
            * PERSP
            * ORTHO
        camera_lens: float
            
        # todo lock camera to view
        
        '''
        if not camera_type:
            camera_type = self.camera_type
        if not camera_lens:
            camera_lens = self.camera_lens
        if self.camera_name in bpy.data.objects:
            camera = bpy.data.objects[self.camera_name]
        else:
            if self.camera_name in bpy.data.cameras:
                camera_data = bpy.data.cameras[self.camera_name]
            else:
                camera_data = bpy.data.cameras.new(self.camera_name)
            camera = bpy.data.objects.new(self.camera_name, camera_data)
            self.coll.objects.link(camera)
        camera.data.lens = camera_lens
        camera.data.dof.aperture_fstop = self.fstop
        camera.data.type = camera_type
        target = self.camera_target
        if self.ortho_scale and camera_type == 'ORTHO':
            camera.data.ortho_scale = self.ortho_scale
        camera.location = Vector(self.camera_loc)
        self.look_at(camera, target, roll = 0.0)
        bpy.context.scene.camera = camera
    def set_light(self):
        '''

        light_type: str
            * POINT, Omnidirectional point light source.
            * SUN, Constant direction parallel ray light source.
            * SPOT, Directional cone light source.
            * AREA, Directional area light source.
        light_energy: float

        # track to camera

        '''
        # check light exist or not
        if self.light_name in bpy.data.objects:
            light = bpy.data.objects[self.light_name]
        else:
            if self.light_name in bpy.data.lights:
                light_data = bpy.data.lights[self.light_name]
            else:
                light_data = bpy.data.lights.new(self.light_name, type = self.light_type)
            light = bpy.data.objects.new(self.light_name, light_data)
            self.coll.objects.link(light)
        light.data.type = self.light_type
        light.data.energy = self.light_energy
        light.location = Vector(self.camera_loc)
        light.data.use_nodes = True
        light.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
        if self.lock_light_to_camera:
            light.constraints.new(type = 'COPY_LOCATION')
            light.constraints["Copy Location"].target = self.camera
            light.constraints.new(type = 'COPY_ROTATION')
            light.constraints["Copy Rotation"].target = self.camera
        else:
            light.location = Vector(self.light_loc)
            self.look_at(light, self.camera_target)
    def look_at(self, obj, target, roll=0):
        """
        Rotate obj to look at target
        """
        if not isinstance(target, Vector):
            target = Vector(target)
        loc = obj.location
        direction = target - loc
        quat = direction.to_track_quat('-Z', 'Y')
        quat = quat.to_matrix().to_4x4()
        rollMatrix = Matrix.Rotation(roll, 4, 'Z')
        loc = loc.to_tuple()
        obj.matrix_world = quat @ rollMatrix
        obj.location = loc
    def motion_blur(self):
        bpy.context.scene.eevee.use_motion_blur = self.use_motion_blur
        bpy.context.scene.eevee.motion_blur_position = self.motion_blur_position
        bpy.context.scene.eevee.motion_blur_steps = self.motion_blur_steps
        bpy.context.scene.eevee.motion_blur_shutter = self.motion_blur_shutter
        bpy.context.scene.cycles.use_motion_blur = self.use_motion_blur
        bpy.context.scene.cycles.motion_blur_position = self.motion_blur_position
        bpy.context.scene.cycles.motion_blur_shutter = self.motion_blur_shutter

    def set_direction(self, direction = (0, 0, 1), margin = None, canvas = None):
        """
        """
        from batoms.tools import get_canvas
        batoms = self.batoms
        atoms = batoms.get_atoms_with_boundary()
        if not margin:
            sizes = [ba.size.max() for ba in batoms.batoms.values()]
            margin = max(sizes) + 0.5
        if canvas is None:
            canvas, canvas1 = get_canvas(atoms = atoms, direction = direction, margin = margin)
        else:
            if isinstance(canvas, (int, float)):
                canvas = np.array([[0, 0, 0], [canvas, canvas, canvas]])
            canvas1 = canvas
        camera_data = batoms.calc_camera_data(canvas, canvas1, direction = direction)
        self.set_parameters(camera_data)
    def run(self, direction = None, canvas = None, **kwargs):
        """
        """
        from batoms.butils import lock_camera_to_view
        if direction:
            self.set_direction(direction, canvas = canvas)
        if self.lock_camera_to_view:
            lock_camera_to_view(True)
        self.set_parameters(kwargs)
        self.set_camera()
        self.set_light()
        if self.use_motion_blur:
            self.motion_blur()
        self.prepare()
        if self.save_to_blend:
            print('saving to {0}.blend'.format(self.output))
            bpy.ops.wm.save_as_mainfile('EXEC_SCREEN', filepath = '{0}.blend'.format(self.output))
        elif self.run_render:
            bpy.ops.render.render(write_still = 1, animation = self.animation)
    def prepare(self, ):
        self.scene.frame_set(self.frame)
        if self.output is None:
            self.output = self.batoms.label
        self.directory = os.path.split(self.output)[0]
        if self.directory and not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.scene.render.image_settings.file_format = 'PNG'
        self.scene.render.engine = self.engine.upper()
        if self.engine.upper() == 'BLENDER_WORKBENCH':
            bpy.data.scenes['Scene'].display.shading.studio_light = '%s'%self.studiolight
        if self.engine.upper() == 'BLENDER_EEVEE':
            self.scene.eevee.taa_render_samples = self.num_samples
            self.scene.cycles.use_denoising = True
        if self.engine.upper() == 'CYCLES':
            self.scene.cycles.samples = self.num_samples
            self.scene.cycles.use_denoising = True
            if self.gpu:
                self.scene.cycles.device = 'GPU'
                prefs = bpy.context.preferences.addons['cycles'].preferences
                for device in prefs.devices:
                    device.use = True
                self.scene.render.tile_x = 256
                self.scene.render.tile_y = 256
        self.scene.render.film_transparent = self.transparent
        self.scene.render.resolution_x = self.resolution_x
        self.scene.render.resolution_y = int(self.resolution_x*self.ratio)
        self.scene.render.filepath = '{0}'.format(self.output)
    def export(self, filename = 'blender-ase.obj'):
        if filename.split('.')[-1] == 'obj':
            bpy.ops.export_scene.obj(filename)
        if filename.split('.')[-1] == 'x3d':
            bpy.ops.export_scene.x3d(filename)
        if filename.split('.')[-1] == 'xyz':
            self.export_xyz(filename)

    def render_move_camera(self, filename, loc1, loc2, n):
        from batoms.tools import getEquidistantPoints
        locs = getEquidistantPoints(loc1, loc2, n)
        i = 0
        for loc in locs:
            bpy.data.objects['Camera'].location = loc
            self.render(self, filename + '_{0:3d}'.format(i))
            i += 1
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s += 'Engine: %s \n'%(self.engine)
        s += 'Camera: %s %s %s \n'%(self.camera_type, self.camera_loc, self.ortho_scale)
        s += 'Light: %s %s \n'%(self.light_type, self.light_energy)
        s += '-'*60 + '\n'
        return s