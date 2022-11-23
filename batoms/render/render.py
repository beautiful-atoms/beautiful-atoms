import bpy
import os
import numpy as np
from batoms.base.collection import BaseCollection
from batoms.render.light import Lights
from batoms.render.camera import Camera
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


subcollections = ['light']


class Render(BaseCollection):
    def __init__(self, label='batoms',
                 batoms=None,
                 viewport=[0, 0, 1],
                 engine='EEVEE',
                 output='batoms.png',
                 animation=False,
                 gpu=False,
                 run_render=True,
                 resolution=[1000, 758],
                 transparent=True,
                #  compute_device_type='CUDA',
                 studiolight='Default',
                 samples=64,
                 ):
        """Rendering Batoms object using blender.

        Object to render batoms object. A Render object has
        one camera, a list of lights.

        Parameters:

        label: str
            Name for the collection in Blender.
        batoms: None or Batoms object
        viewport: array
            The direction of the viewport windows.
        engine: str
            enum in ['BLENDER_WORKBENCH', 'BLENDER_EEVEE', 'CYCLES']
        output: str:
            filepath for the output image
        animation: bool
        gpu: bool
        run_render: bool
        transparent: bool
        resolution: list of 2 ints
        compute_device_type: str
            enum in ['NONE', 'CUDA', 'OPENGL', 'METAL', 'OPTIX'], depending
            on computer.
        studiolight: str
            enum in ['Default', 'basic.sl', 'outdoor.sl',
                    'paint.sl', 'rim.sl', 'studio.sl']
        samples: int
            Default 64
        """
        self.label = label
        self.coll_name = '%s_render' % label
        BaseCollection.__init__(self, coll_name=self.coll_name)
        self.batoms = batoms
        self.camera_name = '%s_camera' % self.label
        # self.compute_device_type = compute_device_type
        coll = bpy.data.collections.get(self.coll_name)
        # load from collection
        self.output = output
        if coll and coll.Brender.flag:
            self.lights = Lights(label)
            self.camera = Camera(label)
            # coll.Brender.viewport = viewport
            # coll.Brender.animation = animation
            # coll.Brender.run_render = run_render
        else:
            self.set_collection(viewport, animation=animation,
                                run_render=run_render)
            self.transparent = transparent
            self.resolution = resolution
            self.gpu = gpu
            self.engine = engine
            self.studiolight = studiolight
            self.samples = samples
            self.camera = Camera(label, coll=self.coll)
            self.lights = Lights(label)
            self.lights.add('Default', lock_to_camera=True)

    def set_collection(self, viewport, animation=False,
                       run_render=True):
        """
        build main collection and its child collections.
        """
        if bpy.data.collections.get(self.coll_name):
            raise Exception("Failed, the name %s already in use!" %
                            self.coll_name)
        coll = bpy.data.collections.new(self.coll_name)
        self.scene.collection.children.link(coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new(
                '%s_%s' % (self.label, sub_name))
            coll.children.link(subcoll)
        coll.Brender.flag = True
        coll.Brender.viewport = viewport
        coll.Brender.animation = animation
        coll.Brender.run_render = run_render

    def set(self, **kwargs):
        """
        Set parameters like set(key1=value1, key2=value2, ...).
        """
        for key, value in kwargs.items():
            setattr(self, key, value)

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
        elif engine.upper() == 'CYCLES':
            self.scene.cycles.use_denoising = True
        self.scene.render.engine = engine.upper()

    @property
    def compute_device_type(self):
        return self.get_compute_device_type()

    @compute_device_type.setter
    def compute_device_type(self, compute_device_type):
        self.set_compute_device_type(compute_device_type)

    def get_compute_device_type(self):
        return bpy.context.preferences.addons["cycles"].preferences.compute_device_type

    def set_compute_device_type(self, compute_device_type):
        bpy.context.preferences.addons["cycles"].preferences.compute_device_type = compute_device_type.upper()

    @property
    def viewport(self):
        return np.array(self.coll.Brender.viewport)

    @viewport.setter
    def viewport(self, viewport):
        if viewport is not None:
            # viewport = viewport/np.linalg.norm(viewport)
            self.coll.Brender.viewport = viewport
            self.update_camera()
            self.update_light()

    @property
    def distance(self):
        return self.coll.Brender.distance

    @distance.setter
    def distance(self, distance):
        if distance is not None:
            self.coll.Brender.distance = distance
            self.update_camera()
            self.update_light()

    @property
    def center(self):
        return np.array(self.coll.Brender.center)

    @center.setter
    def center(self, center):
        if center is not None:
            self.coll.Brender.center = center
            self.update_camera()
            self.update_light()

    @property
    def gpu(self):
        return self.get_gpu()

    @gpu.setter
    def gpu(self, gpu):
        self.set_gpu(gpu)

    def get_gpu(self):
        return self.coll.Brender.gpu

    def set_gpu(self, gpu):
        self.coll.Brender.gpu = gpu
        if gpu:
            self.scene.cycles.device = 'GPU'
            # bpy.context.preferences.addons["cycles"].preferences.get_devices()
            cprefs = bpy.context.preferences.addons["cycles"].preferences
            for device in cprefs.devices:
                device["use"] = 1  # Using all devices, include GPU and CPU
            # auto set compute_device_type, depend on the computer
            for compute_device_type in ('METAL', 'CUDA', 'OPENCL', 'NONE'):
                try:
                    cprefs.compute_device_type = compute_device_type
                    break
                except TypeError:
                    pass

    @property
    def run_render(self):
        return self.coll.Brender.run_render

    @run_render.setter
    def run_render(self, run_render):
        self.coll.Brender.run_render = run_render

    @property
    def animation(self):
        return self.coll.Brender.animation

    @animation.setter
    def animation(self, animation):
        self.coll.Brender.animation = animation

    @property
    def studiolight(self):
        return bpy.context.scene.display.shading.studio_light

    @studiolight.setter
    def studiolight(self, studiolight):
        bpy.context.scene.display.shading.studio_light = \
            '%s' % studiolight

    @property
    def frame(self):
        return bpy.context.scene.display.shading.studio_light

    @frame.setter
    def frame(self, frame):
        self.scene.frame_set(frame)

    @property
    def samples(self):
        return self.scene.cycles.samples

    @samples.setter
    def samples(self, samples):
        self.scene.cycles.samples = samples

    @property
    def transparent(self):
        return self.scene.render.film_transparent

    @transparent.setter
    def transparent(self, transparent):
        self.scene.render.film_transparent = transparent

    @property
    def resolution(self):
        return [self.scene.render.resolution_x,
                self.scene.render.resolution_y]

    @resolution.setter
    def resolution(self, resolution):
        self.scene.render.resolution_x = int(resolution[0])
        self.scene.render.resolution_y = int(resolution[1])

    @property
    def use_motion_blur(self):
        if self.engine == 'BLENDER_EEVEE':
            return bpy.context.scene.eevee.use_motion_blur
        elif self.engine == 'BLENDER_CYCLES':
            return bpy.context.scene.cycles.use_motion_blur

    @use_motion_blur.setter
    def use_motion_blur(self, use_motion_blur):
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

    def set_viewport_distance_center(self, center=None,
                                     padding=None, canvas=None):
        """
        Calculate canvas and direction
        """
        batoms = self.batoms
        if padding is None:
            padding = max(batoms.size) + 0.5
        if center is None:
            center = batoms.get_center_of_geometry()
        self.center = center
        if canvas is None:
            width, height, depth = \
                batoms.get_canvas_box(direction=self.viewport,
                                      padding=padding)
        else:
            width = canvas[0]
            height = canvas[1]
            depth = canvas[2]
        if self.distance < 0:
            self.distance = max(10, depth)
        self.update_camera(width, height, depth/2)
        self.update_light()

    def update_camera(self, width=None, height=None, depth=0):
        if width is not None:
            self.camera.ortho_scale = max(width, height)
            # sensor_with/lens = width/distance
            self.camera.lens = self.camera.obj.data.sensor_width *\
                (self.distance - depth)/width
            self.resolution = [self.resolution[0],
                               self.resolution[0]*height/width]
        self.camera.location = self.center  \
            + self.viewport*self.distance
        self.camera.look_at = self.center

    def update_light(self):
        """
        Calculate canvas and direction
        """
        from batoms.utils.butils import lock_to
        # plane
        # light
        for name, light in self.lights.lights.items():
            if light.lock_to_camera:
                lock_to(self.lights['Default'].obj, self.camera.obj)
            else:
                lock_to(light.obj, None)
                light.direction = light.direction
                # calculate direction in origial axis

    def init(self, batoms=None):
        if batoms is not None:
            self.batoms = batoms
        bpy.context.scene.camera = self.camera.obj
        self.distance = -1
        self.set_viewport_distance_center(
            center=None, padding=None, canvas=None)

    def run(self, batoms,
            center=None,
            padding=None,
            canvas=None,
            **kwargs):
        """
        render the model and export result
        """
        self.batoms = batoms
        self.set(**kwargs)
        self.set_viewport_distance_center(center=center,
                                          canvas=canvas, padding=padding)
        bpy.context.scene.camera = self.camera.obj
        #
        directory = os.path.split(self.output)[0]
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        self.scene.render.image_settings.file_format = 'PNG'
        if self.run_render:
            self.scene.frame_end = self.batoms.nframe
            self.scene.render.filepath = '{0}'.format(self.output)
            bpy.ops.render.render(
                write_still=1, animation=self.coll.Brender.animation)
        else:
            print('saving to {0}.blend'.format(self.output))
            bpy.ops.wm.save_as_mainfile('EXEC_SCREEN',
                                        filepath='{0}.blend'.format(self.output))

    def export(self, filename='blender-ase.obj'):
        if filename.split('.')[-1] == 'obj':
            bpy.ops.export_scene.obj(filename)
        if filename.split('.')[-1] == 'x3d':
            bpy.ops.export_scene.x3d(filename)
        if filename.split('.')[-1] == 'xyz':
            self.export_xyz(filename)

    def render_move_camera(self, filename, loc1, loc2, n):
        from batoms.utils import getEquidistantPoints
        locs = getEquidistantPoints(loc1, loc2, n)
        i = 0
        for loc in locs:
            bpy.data.objects['Camera'].location = loc
            self.render(self, filename + '_{0:3d}'.format(i))
            i += 1

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s += 'Label:   %s \n' % (self.label)
        if self.batoms is not None:
            s += 'Batoms:   %s \n' % (self.batoms.label)
        s += 'Engine:   %s \n' % (self.engine)
        s += 'Viewport: %s \n' % (self.viewport)
        s += 'Camera:   type = %s, location = %s, ortho_scale = %s \n' \
            % (self.camera.type, self.camera.location,
               self.camera.ortho_scale)
        s += 'Light: \n'
        s += self.lights.__repr__()
        s += '-'*60 + '\n'
        return s
