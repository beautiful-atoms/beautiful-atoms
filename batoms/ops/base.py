"""
This module defines abstract Operator classes.
"""

import bpy
from bpy.types import Operator
from batoms import Batoms


class OperatorBatoms(Operator):
    """
    """
    bl_options = {'REGISTER', 'UNDO'}

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
        return {'FINISHED'}


class OperatorBatomsEdit(OperatorBatoms):
    """
    """
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER' and context.mode in {'EDIT_MESH'}
        else:
            return False
