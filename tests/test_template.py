import bpy
import pytest
from ase.build import bulk
from batoms import Batoms
from batoms.bio.bio import read
import numpy as np
import os
try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}

skip_test = bool(os.environ.get("NOTEST_CUBE", 0))


def test_template():
    bpy.ops.batoms.delete()
    au = bulk("Au", cubic=True)
    au = Batoms("au", from_ase=au)
    au.template.settings['test'] = {"prop1": 1.1}
    print(au.template.settings)
    au.template.draw()


def test_template_ops():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au")
    au = Batoms('Au')
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.template.template_add(name='test')
    au.template.settings['test'].scale = 3.0
    assert len(au.template.settings) == 1
    bpy.ops.template.template_draw()
    bpy.ops.template.template_remove(name="test")
    assert len(au.template.settings) == 0
    print(au.template.settings)
    bpy.ops.template.template_draw()

def test_gui():
    """latticeplane panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert bpy.context.scene.Btemplate.show == True
    bpy.context.scene.Btemplate.show = False
    assert au.template.show == False

def test_gui_uilist():
    """latticeplane panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert au.coll.Btemplate.ui_list_index==0
    bpy.ops.template.template_add(name='test1')
    bpy.ops.template.template_add(name='test2')
    assert au.coll.Btemplate.ui_list_index==1


if __name__ == "__main__":
    test_template()
    print("\n Template: All pass! \n")
