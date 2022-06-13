"""
"""
import bpy
import numpy as np
from ase.cell import Cell
from batoms.utils.butils import object_mode, clean_coll_objects
from batoms.base.object import ObjectGN
from batoms.draw import draw_cylinder
# from time import time


class Bcell(ObjectGN):
    def __init__(self, label,
                 array=None,
                 location=np.array([0, 0, 0]),
                 color=(0.0, 0.0, 0.0, 1.0),
                 batoms=None,
                 ) -> None:
        """
        Unit cell of three dimensions.
        ver: 3x3 verlike object
          The three cell vectors: cell[0], cell[1], and cell[2].
        """
        self.label = label
        name = 'cell'
        self.batoms = batoms
        ObjectGN.__init__(self, label, name)
        self.edges = [[0, 3], [0, 1], [4, 2], [4, 1],
                      [3, 5], [2, 6], [7, 5], [7, 6],
                      [0, 2], [3, 6], [1, 5], [4, 7]
                      ]
        self.width = 0.03
        self.color = color
        if array is not None:
            self.build_object(array, location)
        # if self.show_axes:
            # self.show_axes = True

    def build_object(self, array, location):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        array = np.array(array)
        if array is None:
            positions = np.zeros([8, 3])
            self._frames = {'positions': [positions]}
        elif array.shape in [(3, ), (3, 1), (6, ), (6, 1), (3, 3)]:
            cell = Cell.new(array)
            positions = self.array2verts(cell.array)
            self._frames = {'positions': [positions]}
        elif array.shape == (8, 3):
            positions = array - array[0]
            location = array[0]
            self._frames = {'positions': [positions]}
        elif len(array.shape) == 3:
            positions = array[0]
            self._frames = {'positions': [positions]}
        else:
            raise Exception('Shape of cell %s is not corrected!' % array.shape)

        self.delete_obj(self.obj_name)
        mesh = bpy.data.meshes.new(self.obj_name)
        mesh.from_pydata(positions, self.edges, [])
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj = bpy.data.objects.new(self.obj_name, mesh)
        obj.data = mesh
        obj.location = location
        obj.batoms.type = 'CELL'
        obj.batoms.label = self.label
        if self.batoms is not None:
            self.batoms.coll.objects.link(obj)
        else:
            bpy.data.collections['Collection'].objects.link(obj)
        self.build_geometry_node()
        self.set_frames(self._frames)
        bpy.context.view_layer.update()

    def build_geometry_node(self):
        """
        """
        from batoms.utils.butils import get_nodes_by_name, compareNodeType, \
                            build_modifier
        name = 'GeometryNodes_%s_cell' % self.label
        modifier = build_modifier(self.obj, name)
        # ------------------------------------------------------------------
        # select attributes
        gn = modifier
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = gn.node_group.nodes[1]
        # ------------------------------------------------------------------
        # transfer first 4 positions of cell
        PositionCell = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_PositionCell' % (self.label),
                                         'GeometryNodeInputPosition')
        TransferCells = []
        for i in range(4):
            TransferCell = get_nodes_by_name(gn.node_group.nodes,
                                             '%s_TransferCell_%s' % (
                                                 self.label, i),
                                             'GeometryNodeAttributeTransfer')
            TransferCell.mapping = 'INDEX'
            TransferCell.data_type = 'FLOAT_VECTOR'
            InputInt = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_InputInt_%s' % (self.label, i),
                                         'FunctionNodeInputInt')
            InputInt.integer = i
            gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                    TransferCell.inputs[0])
            gn.node_group.links.new(PositionCell.outputs['Position'],
                                    TransferCell.inputs['Attribute'])
            gn.node_group.links.new(InputInt.outputs[0],
                                    TransferCell.inputs['Index'])
            TransferCells.append(TransferCell)
        # ------------------------------------------------------------------
        VectorAdds = []
        for i in range(4):
            VectorAdd = get_nodes_by_name(gn.node_group.nodes,
                                          '%s_VectorAdd_%s' % (self.label, i),
                                          'ShaderNodeVectorMath')
            VectorAdd.operation = 'ADD'
            VectorAdds.append(VectorAdd)
        gn.node_group.links.new(
            TransferCells[1].outputs[0], VectorAdds[0].inputs[0])
        gn.node_group.links.new(
            TransferCells[2].outputs[0], VectorAdds[0].inputs[1])
        gn.node_group.links.new(
            TransferCells[1].outputs[0], VectorAdds[1].inputs[0])
        gn.node_group.links.new(
            TransferCells[3].outputs[0], VectorAdds[1].inputs[1])
        gn.node_group.links.new(
            TransferCells[2].outputs[0], VectorAdds[2].inputs[0])
        gn.node_group.links.new(
            TransferCells[3].outputs[0], VectorAdds[2].inputs[1])
        gn.node_group.links.new(
            TransferCells[3].outputs[0], VectorAdds[3].inputs[0])
        gn.node_group.links.new(
            VectorAdds[0].outputs[0], VectorAdds[3].inputs[1])
        # calculate unit cell vector
        VectorSubtracts = []
        for i in range(4):
            VectorSubtract = get_nodes_by_name(gn.node_group.nodes,
                                               '%s_VectorSubtract_%s' % (
                                                   self.label, i),
                                               'ShaderNodeVectorMath')
            VectorSubtract.operation = 'SUBTRACT'
            gn.node_group.links.new(
                VectorAdds[i].outputs[0], VectorSubtract.inputs[0])
            gn.node_group.links.new(
                TransferCells[0].outputs[0], VectorSubtract.inputs[1])
            VectorSubtracts.append(VectorSubtract)
        # for 7
        VectorSubtract = get_nodes_by_name(gn.node_group.nodes,
                                           '%s_VectorSubtract_%s' % (
                                               self.label, 4),
                                           'ShaderNodeVectorMath')
        VectorSubtract.operation = 'SUBTRACT'
        gn.node_group.links.new(
            VectorSubtracts[3].outputs[0], VectorSubtract.inputs[0])
        gn.node_group.links.new(
            TransferCells[0].outputs[0], VectorSubtract.inputs[1])
        VectorSubtracts.append(VectorSubtract)
        # set positions
        SetPositions = []
        IndexCell = get_nodes_by_name(gn.node_group.nodes,
                                      '%s_IndexCell' % (self.label),
                                      'GeometryNodeInputIndex')
        for i in range(4):
            SetPosition = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_SetPosition_%s' % (
                                                self.label, i),
                                            'GeometryNodeSetPosition')
            CompareSelect = get_nodes_by_name(gn.node_group.nodes,
                                              'select_%s_%s' % (
                                                  self.label, i),
                                              compareNodeType)

            CompareSelect.operation = 'EQUAL'
            CompareSelect.inputs[1].default_value = i + 4
            gn.node_group.links.new(
                IndexCell.outputs[0], CompareSelect.inputs[0])
            gn.node_group.links.new(
                CompareSelect.outputs[0], SetPosition.inputs['Selection'])
            SetPositions.append(SetPosition)
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                SetPositions[0].inputs['Geometry'])
        gn.node_group.links.new(SetPositions[0].outputs['Geometry'],
                                SetPositions[1].inputs['Geometry'])
        gn.node_group.links.new(SetPositions[1].outputs['Geometry'],
                                SetPositions[2].inputs['Geometry'])
        gn.node_group.links.new(SetPositions[2].outputs['Geometry'],
                                SetPositions[3].inputs['Geometry'])
        gn.node_group.links.new(SetPositions[3].outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])
        gn.node_group.links.new(VectorSubtracts[0].outputs[0],
                                SetPositions[0].inputs['Position'])
        gn.node_group.links.new(VectorSubtracts[1].outputs[0],
                                SetPositions[1].inputs['Position'])
        gn.node_group.links.new(VectorSubtracts[2].outputs[0],
                                SetPositions[2].inputs['Position'])
        gn.node_group.links.new(VectorSubtracts[4].outputs[0],
                                SetPositions[3].inputs['Position'])

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0:
            return
        name = self.label
        obj = self.obj
        self.set_obj_frames(name, obj, frames["positions"])

    def build_cell_cylinder(self):
        #
        cell_cylinder = {'lengths': [],
                         'centers': [],
                         'normals': [],
                         'vertices': 16,
                         'width': self.width,
                         'color': self.color,
                         'battr_inputs': {},
                         }
        if np.max(abs(self.positions)) < 1e-6:
            return cell_cylinder
        for e in self.edges:
            center = (self.positions[e[0]] + self.positions[e[1]])/2.0
            vec = self.positions[e[0]] - self.positions[e[1]]
            length = np.linalg.norm(vec)
            nvec = vec/length
            cell_cylinder['lengths'].append(length)
            cell_cylinder['centers'].append(center)
            cell_cylinder['normals'].append(nvec)
        return cell_cylinder

    def __repr__(self) -> str:
        numbers = np.round(self.array, 3).tolist()
        s = 'Cell({})'.format(numbers)
        return s

    def __getitem__(self, index):
        return self.local_array[index]

    def __setitem__(self, index, value):
        """Set unit cell vectors.

        Parameters:

        Examples:

        """
        array = self.array
        array[index] = value
        positions = self.array2verts(array)
        self.local_positions = positions

    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.local_array

    @property
    def local_array(self):
        return self.get_local_array()

    def get_local_array(self):
        """
        In this case, the origin is translated to vertices[0].
        While the orientation is reserved. 

        Returns:
            (3x3 local_array): The array of cell.
        """
        positions = self.local_positions
        local_array = np.array([positions[1] - positions[0],
                                positions[2] - positions[0],
                                positions[3] - positions[0]])
        return local_array

    @property
    def array(self):
        return self.get_array()

    def get_array(self):
        """
        In this case, the origin is translated to vertices[0].
        While the orientation is reserved. 

        Returns:
            (3x3 array): The array of cell.
        """
        positions = self.positions
        array = np.array([positions[1] - positions[0],
                         positions[2] - positions[0],
                         positions[3] - positions[0]])
        return array

    @property
    def origin(self):
        return self.positions[0]

    def array2verts(self, array):
        """
        """
        positions = np.array([[0, 0, 0],
                              [1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1],
                              [1, 1, 0],
                              [1, 0, 1],
                              [0, 1, 1],
                              [1, 1, 1],
                              ])
        positions = np.dot(positions, array)
        return positions

    def copy(self, label):
        object_mode()
        cell = Bcell(label, array=self.array, location=self.obj.location)
        return cell

    def repeat(self, m):
        self[:] = np.array([m[c] * self.local_array[c] for c in range(3)])

    @property
    def length(self):
        length = np.linalg.norm(self.local_array, axis=1)
        return length

    @property
    def reciprocal(self):
        from math import pi
        b1 = 2*pi/self.volume*np.cross(self[1], self[2])
        b2 = 2*pi/self.volume*np.cross(self[2], self[0])
        b3 = 2*pi/self.volume*np.cross(self[0], self[1])
        return np.array([b1, b2, b3])

    @property
    def volume(self):
        return np.dot(self[0], np.cross(self[1], self[2]))

    @property
    def center(self):
        """Center of unit cell.
        """
        return (self.array[0] + self.array[1] + self.array[2])/2.0

    def draw_curve_from_vertices_nurbs(self, label, vertices, coll,
                                       ):
        """
        """
        name = '%s_cell' % label
        crv = bpy.data.curves.new(name, 'CURVE')
        crv.dimensions = '3D'
        crv.resolution_u = 30
        crv.fill_mode = 'FULL'
        spline = crv.splines.new(type='POLY')
        nvert = len(vertices)
        spline.points.add(nvert-1)
        # vertices = np.append(vertices, np.zeros((nvert, 1)), axis = 1)
        vertices = np.append(vertices, np.ones((len(vertices), 1)), axis=1)
        vertices = vertices.reshape(-1, 1)
        spline.points.foreach_set('co', vertices)
        crv.bevel_mode = 'OBJECT'
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = bpy.data.objects.new(name, crv)
        #
        obj.data.materials.append(self.material)
        coll.objects.link(obj)

    def build_materials(self, label, color,
                        node_type='Principled BSDF',
                        use_smooth=True,
                        node_inputs=None,
                        material_style='plastic',
                        backface_culling=True):
        """
        """
        from batoms.material import create_material
        name = '%s_cell' % (label)
        if name not in bpy.data.materials:
            create_material(name,
                            color=color,
                            node_inputs=node_inputs,
                            material_style=material_style,
                            backface_culling=True)

    def build_bevel_object(self, label, radius):
        # Create bevel control curve.
        name = '%s_cell_bevel_object' % (label)
        bpy.ops.curve.primitive_bezier_circle_add(radius=radius)
        bevel_control = bpy.context.active_object
        bevel_control.data.name = bevel_control.name = '%s_bevel' % name
        # Set the main curve's bevel control to the bevel control curve.
        self.obj.data.bevel_object = self.bevel_control

    @property
    def bevel_object(self):
        name = '%s_cell_bevel_object' % (self.label)
        return bpy.data.objects.get(name)

    def draw(self):
        """Draw unit cell
        """
        object_mode()
        name = '%s_%s_%s' % (self.label, 'cell', 'cylinder')
        self.delete_obj(name)
        cell_cylinder = self.build_cell_cylinder()
        if self.batoms is not None:
            coll = self.batoms.coll
        else:
            coll = bpy.data.collections["Collection"]
        obj = draw_cylinder(name=name,
                            datas=cell_cylinder,
                            coll=coll
                            )
        if self.batoms is not None:
            obj.parent = self.batoms.obj

    @property
    def obj_cylinder(self):
        obj = bpy.data.objects.get('%s_cylinder' % self.obj_name)
        if obj is None:
            raise KeyError('%s object is not exist.' % self.obj_name)
        return obj


    @property
    def draw_crystal_axes(self):
        """draw_crystal_axes object.
        Shoud be a global variable.
        """
        from batoms.draw.draw_screen import DrawCrystalAxes
        dns = bpy.app.driver_namespace
        name = "{}_crystal_axes".format(self.label)
        if self.label in dns:
            return dns[self.label]
        context = bpy.context
        dns[name] = DrawCrystalAxes(context, name)
        return dns[name]

    @property
    def show_axes(self):
        if self.batoms is not None:
            return self.batoms.coll.batoms.show_axes
        else:
            return False

    @show_axes.setter
    def show_axes(self, show_axes):
        self.batoms.coll.batoms.show_axes = show_axes
        positions = self.array
        self.draw_crystal_axes.remove_handle()
        if show_axes:
            self.draw_crystal_axes.add_handle(positions)
