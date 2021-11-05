import bpy
from batoms.butils import get_selected_vertices
from batoms import Batom
import numpy as np
from ase.geometry.geometry import get_distances, get_angles, get_dihedrals

class EDIT_MESH_OT_record_selection(bpy.types.Operator):
    """Records the selection order while running and when finished with ESC """
    bl_idname = "batom.record_selection"
    bl_label = "Record Selection of batom"

    origSel: None
    selOrder: None

    @classmethod
    def poll(cls, context):
        return context.mode in  {'EDIT_MESH'}

    def modal(self, context, event):
        selected_vertices = get_selected_vertices()
        # Only add one more vertices
        nsel = []
        for label, species, name, index in selected_vertices:
            for i in index:
                if (name, np.array([i])) not in self.selOrder:
                    nsel.append((name, np.array([i])))
        if len(nsel) == 1:
            self.report({'INFO'}, 'New vertices was selected')
            self.selOrder.append(nsel[0])
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            cell = None
            pbc = None
            positions = np.array([]).reshape(-1, 3)
            print(self.selOrder)
            for name, index in self.selOrder:
                batom = Batom(label = name)
                positions = np.append(positions, batom.positions[index], axis = 0)
            bapanel = context.scene.bapanel
            if len(positions) == 2:
                results = get_distances([positions[0]], 
                            [positions[1]], 
                            cell=cell, pbc=pbc)[1]
            elif len(positions) == 3:
                v12 = positions[0] - positions[1]
                v32 = positions[2] - positions[1]
                results =  get_angles([v12], [v32], cell=cell, pbc=pbc)
            elif len(positions) == 4:
                v0 = positions[1] - positions[0]
                v1 = positions[2] - positions[1]
                v2 = positions[3] - positions[2]
                results =  get_dihedrals([v0], [v1], [v2], cell=cell, pbc=pbc)
            else:
                return {'CANCELLED'}
            # update measurement value
            results.shape = (-1,)
            results = [str(round(float(i), 2)) for i in results]
            results = ' '.join(results)
            bapanel.measurement = results
            return {'CANCELLED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.selOrder = []
        selected_vertices = get_selected_vertices()
        if len(selected_vertices) == 1:
            if selected_vertices[0][3] == 1:
                self.selOrder.append(selected_vertices[0])
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}



def add_edit_mesh_button(self, context):
    self.layout.operator(
        EDIT_MESH_OT_record_selection.bl_idname,
        text="Select atom",
        icon='RESTRICT_SELECT_OFF')

