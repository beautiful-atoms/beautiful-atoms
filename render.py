"""
Rendering Batoms object using blender.

"""
import bpy
from mathutils import Vector, Matrix
import os
import numpy as np
from batoms.light import Lights
from batoms.camera import Camera

default_render_settings = {
        'transparent': True,  # transparent background
        'ratio': 1,
        'world': True,
        # 'engine': 'BLENDER_EEVEE', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
        'frame': None,
        'run_render': True,
        'output': None,
        'animation': False,
        'save_to_blend': False,
        'gpu': True,
        'compute_device_type': 'CUDA',
        'num_samples': 32,
        }

subcollections = ['light']

class Render():
    """
    Object to render atomic structure.

    Has one camera, a list of light

    """
    #
    def __init__(self, label, batoms = None, from_coll = False, **kwargs):
        self.label = label
        self.name = '%s_render'%label
        self.batoms = batoms
        self.camera_name = 'camera_%s'%self.label
        self.set_parameters(default_render_settings)
        if from_coll:
            self.lights = Lights(label, from_coll = True)
            self.camera = Camera(label, from_coll = True)
        else:
            self.set_collection()
            default_render_settings.update(kwargs)
            self.clean_default()
            self.camera = Camera(label, coll = self.coll)
            self.lights = Lights(label)
            self.lights.add('Default', lock_to_camera = True)
            if self.world:
                self.set_world()
    def set_collection(self):
        """
        build main collection and its child collections.
        """
        if self.name not in bpy.data.collections:
            coll = bpy.data.collections.new(self.name)
            bpy.data.collections[self.label].children.link(coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s'%(self.label, sub_name))
            self.coll.children.link(subcoll)
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
    def set_world(self, color = [0.2, 0.2, 0.2, 1.0]):
        """
        """
        world = self.scene.world
        world.use_nodes = True
        node_tree = world.node_tree
        node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
        node_tree.nodes["Background"].inputs["Color"].default_value = color
    def motion_blur(self, use_motion_blur):
        """
        'use_motion_blur': False,
        'motion_blur_position': 'START', 
        'motion_blur_steps': 10,
        'motion_blur_shutter': 30.0,
        """
        bpy.context.scene.eevee.use_motion_blur = use_motion_blur
        bpy.context.scene.eevee.motion_blur_position = 'START'
        bpy.context.scene.eevee.motion_blur_steps = 10
        bpy.context.scene.eevee.motion_blur_shutter = 30.0
        bpy.context.scene.cycles.use_motion_blur = use_motion_blur
        bpy.context.scene.cycles.motion_blur_position = 'START'
        bpy.context.scene.cycles.motion_blur_shutter = 30.0
    def set_direction(self, direction = (0, 0, 1), 
                    distance = None, 
                    padding = None, canvas = None):
        """
        Calculate canvas and direction
        """
        from batoms.tools import rotate_frame
        batoms = self.batoms
        if padding is None:
            sizes = [0]
            sizes.extend([ba.size.max() for ba in batoms.batoms.values()])
            padding = max(sizes) + 0.5
        if canvas is None:
            canvas, canvas1 = batoms.get_canvas_box(direction = direction, 
                                        padding = padding)
        else:
            if isinstance(canvas, (int, float)):
                canvas = np.array([[0, 0, 0], [canvas, canvas, canvas]])
            canvas1 = canvas
        center = np.mean(canvas, axis=0)
        width = canvas1[1, 0] - canvas1[0, 0]
        height = canvas1[1, 1] - canvas1[0, 1]
        depth = canvas1[1, 2] - canvas1[0, 2]
        if distance is None:
            distance = max(10, depth*2)
        # camera
        direction = direction/np.linalg.norm(direction)
        self.camera.location = center + direction*distance
        self.camera.ortho_scale = max(width, height)
        self.camera.target = center
        self.ratio = height/width
        # plane
        frame = rotate_frame(direction)
        # light
        for name, light in self.lights.lights.items():
            if light.lock_to_camera:
                self.lights['Default'].lock_to(self.camera.obj)
            else:
                light.lock_to(None)
                # calculate direction in origial axis
                dirction = np.dot(light.direction, frame)
                dirction = dirction/np.linalg.norm(dirction)
                location = center + dirction*distance
                light.location = location
                light.target = center
    def run(self, direction = None, distance = None, canvas = None, 
                padding = None, **kwargs):
        """
        render the model and export result
        """
        if direction:
            self.set_direction(direction, distance = distance, 
                        canvas = canvas, padding = padding)
        self.set_parameters(kwargs)
        if kwargs.get('use_motion_blur'):
            self.motion_blur(kwargs.get('use_motion_blur'))
        self.prepare()
        if self.save_to_blend:
            print('saving to {0}.blend'.format(self.output))
            bpy.ops.wm.save_as_mainfile('EXEC_SCREEN', 
                                filepath = '{0}.blend'.format(self.output))
        elif self.run_render:
            bpy.data.scenes[0].frame_end = self.batoms.nframe
            bpy.ops.render.render(write_still = 1, animation = self.animation)
    def prepare(self, num_samples = 32, studiolight = 'Default',
                resolution_x = 1000, 
                ):
        """
        Set parameters for rendering
        """
        # frame_set(1) is wrong when mesh is modified.
        bpy.context.scene.camera = self.camera.obj
        if self.frame is not None:
            self.scene.frame_set(self.frame)
        if self.output is None:
            self.output = self.label
        self.directory = os.path.split(self.output)[0]
        if self.directory and not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.scene.render.image_settings.file_format = 'PNG'
        self.scene.render.engine = self.engine.upper()
        if self.engine.upper() == 'BLENDER_WORKBENCH':
            bpy.data.scenes['Scene'].display.shading.studio_light = '%s'%self.studiolight
        if self.engine.upper() == 'BLENDER_EEVEE':
            self.scene.eevee.taa_render_samples = num_samples
            self.scene.cycles.use_denoising = True
        if self.engine.upper() == 'CYCLES':
            self.scene.cycles.samples = num_samples
            self.scene.cycles.use_denoising = True
            if self.gpu:
                bpy.context.preferences.addons["cycles"].preferences.compute_device_type  \
                                = self.compute_device_type
                self.scene.cycles.device = 'GPU'
                bpy.context.preferences.addons["cycles"].preferences.get_devices()
                for device in bpy.context.preferences.addons["cycles"].preferences.devices:
                    device["use"] = 1 # Using all devices, include GPU and CPU
        self.scene.render.film_transparent = self.transparent
        self.scene.render.resolution_x = resolution_x
        self.scene.render.resolution_y = int(resolution_x*self.ratio)
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
        s += 'Light: \n'
        s += self.lights.__repr__()
        s += '-'*60 + '\n'
        return s
    
    