#!/bin/bash
# Special settings for testing on CI systems, such as github actions.
# Please run using
# PYTHONPATH=$(pwd) bash test_ci.sh
rm -rf *png
if ! [ -z "$CI" ]
then
    echo "Running inside a CI system, use headleass only."
    export DISPLAY=
    export NOTEST_CUBE=1
fi
blender -b -P test_animation.py
blender -b -P test_batom.py
blender -b -P test_batoms.py
blender -b -P test_bonds.py
blender -b -P test_boundary.py
blender -b -P test_camera.py
blender -b -P test_cell.py
blender -b -P test_light.py
blender -b -P test_ms.py
blender -b -P test_performance.py
blender -b -P test_plane.py
blender -b -P test_plugins.py
blender -b -P test_polyhedrasetting.py
blender -b -P test_real.py
blender -b -P test_render.py
blender -b -P test_ribbon.py
blender -b -P test_select.py
# skip the pytest if NOTEST_CUBE is set
blender -b -P test_isosurfacesetting.py
