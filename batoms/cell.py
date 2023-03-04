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
                 batoms=None,
                 location=np.array([0, 0, 0]),
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
        if array is not None:
            self.build_object(array, location)
        # if self.show_axes:
            # self.show_axes = True

    def build_object(self, array, location):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        array = np.array(array)
        positions = np.zeros([4, 3])
        if array is None:
            self._frames = {'positions': [positions]}
        elif array.shape in [(3, ), (3, 1), (6, ), (6, 1), (3, 3)]:
            cell = Cell.new(array)
            positions[1:4, :] = cell.array
            self._frames = {'positions': [positions]}
        elif len(array.shape) == 3:
            positions = array[0]
            self._frames = {'positions': [positions]}
        else:
            raise Exception('Shape of cell %s is not corrected!' % array.shape)

        self.delete_obj(self.obj_name)
        mesh = bpy.data.meshes.new(self.obj_name)
        mesh.from_pydata(positions, [[0, 1], [0, 2], [0, 3]], [])
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
        # materials
        mat = self.build_materials(self.label, color=self.color)
        obj.data.materials.append(mat)
        self.build_geometry_node()
        self.set_frames(self._frames)
        bpy.context.view_layer.update()

    def build_geometry_node(self):
        """
        """
        from batoms.utils.butils import get_nodes_by_name, build_modifier
        name = 'GeometryNodes_%s_cell' % self.label
        modifier = build_modifier(self.obj, name)
        # ------------------------------------------------------------------
        # select attributes
        gn = modifier
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = gn.node_group.nodes[1]
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        gn.node_group.links.new(
            GroupInput.outputs['Geometry'], JoinGeometry.inputs['Geometry'])
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
                                            'GeometryNodeSampleIndex')
            TransferCell.data_type = 'FLOAT_VECTOR'
            InputInt = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_InputInt_%s' % (self.label, i),
                                         'FunctionNodeInputInt')
            InputInt.integer = i
            gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                    TransferCell.inputs[0])
            gn.node_group.links.new(PositionCell.outputs['Position'],
                                    TransferCell.inputs[3])
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
            TransferCells[1].outputs[2], VectorAdds[0].inputs[0])
        gn.node_group.links.new(
            TransferCells[2].outputs[2], VectorAdds[0].inputs[1])
        gn.node_group.links.new(
            TransferCells[1].outputs[2], VectorAdds[1].inputs[0])
        gn.node_group.links.new(
            TransferCells[3].outputs[2], VectorAdds[1].inputs[1])
        gn.node_group.links.new(
            TransferCells[2].outputs[2], VectorAdds[2].inputs[0])
        gn.node_group.links.new(
            TransferCells[3].outputs[2], VectorAdds[2].inputs[1])
        gn.node_group.links.new(
            TransferCells[3].outputs[2], VectorAdds[3].inputs[0])
        gn.node_group.links.new(
            VectorAdds[0].outputs[0], VectorAdds[3].inputs[1])
        #============================================================
        # In this implementation.
        # Use an object with four vertices to represent a cell, instead of eight vertices.
        # We use curve to mesh methods to draw the edges (cylinder)
        # One does not need to call cell.draw() anymore.
        # circle for profle of curve
        Circle = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_Circle' % self.label,
                                            'GeometryNodeCurvePrimitiveCircle')
        Circle.inputs[4].default_value = self.width
        # faces from vertices
        faces = [[0, 2, 6, 3],
                [1, 4, 7, 5],
                [0, 1, 4, 2],
                [3, 5, 7, 6],
                [0, 1, 5, 3],
                [2, 4, 7, 6],
                ]
        input_nodes = TransferCells + VectorAdds
        for i in range(6):
            # set quadrilaterial
            # maybe in the future, we have curve primitive: cube
            Quadrilateral = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_Quadrilateral_%s' % (
                                                self.label, i),
                                            'GeometryNodeCurvePrimitiveQuadrilateral')
            Quadrilateral.mode = "POINTS"
            CurveToMesh = get_nodes_by_name(gn.node_group.nodes,
                                            '%s_CurveToMesh_%s' % (
                                                self.label, i),
                                            'GeometryNodeCurveToMesh')
            gn.node_group.links.new(Quadrilateral.outputs[0],
                                CurveToMesh.inputs[0])
            gn.node_group.links.new(Circle.outputs[0],
                                CurveToMesh.inputs[1])
            gn.node_group.links.new(CurveToMesh.outputs[0],
                                JoinGeometry.inputs[0])
            for j in range(4):
                k = 2 if faces[i][j] < 4 else 0
                gn.node_group.links.new(input_nodes[faces[i][j]].outputs[k],
                                Quadrilateral.inputs[j + 7])
        setMaterial = get_nodes_by_name(self.gnodes.node_group.nodes,
                                                '%s_setMaterial' % (self.label),
                                                'GeometryNodeSetMaterial')
        setMaterial.inputs[2].default_value = self.material
        gn.node_group.links.new(JoinGeometry.outputs[0],
                                setMaterial.inputs[0])
        gn.node_group.links.new(setMaterial.outputs[0],
                                GroupOutput.inputs[0])

    @property
    def material(self):
        return self.get_material()

    def get_material(self):
        name = '%s_cell' % self.label
        mat = bpy.data.materials.get(name)
        return mat

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0:
            return
        name = self.label
        obj = self.obj
        self.set_obj_frames(name, obj, frames["positions"])

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
        local_positions = self.local_positions
        local_positions[1:4, :][index] = value
        self.local_positions = local_positions
        self.batoms.update_gn_cell()
        self.batoms.boundary.update_gn_cell()

    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.local_array

    @property
    def width(self):
        return self.batoms.coll.batoms.cell.width

    @width.setter
    def width(self, width):
        from batoms.utils.butils import get_nodes_by_name
        self.batoms.coll.batoms.cell.width = width
        Circle = get_nodes_by_name(self.gnodes.node_group.nodes,
                                            '%s_Circle' % (
                                                self.label),
                                            'GeometryNodeCurvePrimitiveCircle')
        Circle.inputs[4].default_value = width

    @property
    def color(self):
        # Viewpoint_color = self.materials[self.main_element].diffuse_color
        # for node in self.material.node_tree.nodes:
        #     if 'Base Color' in node.inputs:
        #         node_color = node.inputs['Base Color'].default_value[:]
        #     if 'Alpha' in node.inputs:
        #         Alpha = node.inputs['Alpha'].default_value
        # color = [node_color[0], node_color[1], node_color[2], Alpha]
        color = self.batoms.coll.batoms.cell.color
        return color

    @color.setter
    def color(self, color):
        from batoms.utils.butils import get_nodes_by_name
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.batoms.coll.batoms.cell.color = color
        self.material.diffuse_color = color
        for node in self.material.node_tree.nodes:
            if 'Base Color' in node.inputs:
                node.inputs['Base Color'].default_value = color
            if 'Alpha' in node.inputs:
                node.inputs['Alpha'].default_value = color[3]
        setMaterial = get_nodes_by_name(self.gnodes.node_group.nodes,
                                                '%s_setMaterial' % (self.label),
                                                'GeometryNodeSetMaterial')
        setMaterial.inputs[2].default_value = self.material

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

    @property
    def edges(self):
        """Edges of the cell
        """
        edge_indices = [[0, 3], [0, 1], [4, 2], [4, 1],
                      [3, 5], [2, 6], [7, 5], [7, 6],
                      [0, 2], [3, 6], [1, 5], [4, 7]
                      ]
        basis = np.array([[0, 0, 0],
                              [1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1],
                              [1, 1, 0],
                              [1, 0, 1],
                              [0, 1, 1],
                              [1, 1, 1],
                              ])
        positions = np.dot(basis, self.array)
        edges = []
        for indices in edge_indices:
            edges.append(positions[indices])
        return edges

    def copy(self, label):
        object_mode()
        cell = Bcell(label, batoms = self.batoms, array=self.array, location=self.obj.location)
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


    def build_materials(self, label, color=None,
                        node_type='Principled BSDF',
                        use_smooth=True,
                        node_inputs=None,
                        material_style='plastic',
                        backface_culling=True):
        """
        """
        from batoms.material import create_material
        name = '%s_cell' % (label)
        if color is None:
            color = [0.2, 0.2, 0.2, 1]
        self.delete_material(name)
        mat = create_material(name,
                        color=color,
                        node_inputs=node_inputs,
                        material_style=material_style,
                        backface_culling=True)
        return mat



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

    def draw(self):
        from batoms.utils import deprecated
        deprecated('"draw" will be deprecated in the furture. The cell is drawn automaticely now.')
        pass

    def as_dict(self):
        """
        """
        data = {
            'array': self.array,
        }
        data.update(self.batoms.coll.batoms.cell.as_dict())
        return data
