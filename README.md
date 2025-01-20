## Beautiful Atoms
[![Build beautiful_atoms main image](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/build_main_image.yml/badge.svg)](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/build_main_image.yml)
[![Unit test](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/ci.yml/badge.svg)](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/unittests.yaml)
[![codecov](https://codecov.io/gh/beautiful-atoms/beautiful-atoms/branch/main/graph/badge.svg)](https://codecov.io/gh/beautiful-atoms/beautiful-atoms)

Batoms is a Python package for editing and rendering atoms and molecule objects using Blender. A Python interface that allows for automating workflows.

Features:

* Model: space-filling, ball-stick, polyhedral, cavity and so on.
* Supported File type: cif, xyz, cube, pdb, json, VASP-out and so on.
* Supported structure: ASE and Pymatgen
* Volumetric data (Isosurface)
* Ribbon diagram for protein
* Site occupancy
* Animation
* GUI
* Support periodic boundary conditions
* Support fetch structures from MaterialProject, Pubchem and RSCB
* ``Flexible``: Python script, run interactively or in the background.
* ``High quality rendering``:  3D models
* ``Free, Open Source``: Easy to download and install.
* ``Cross-platform``: (Linux, Windows, macOS)



## How to use
### **⚠️ Need Blender 4.2.0 or later.**

Please visit the [documentation](https://beautiful-atoms.readthedocs.io/en/latest/) for more information.

## View and edit structure in Jupyter Notebook
Another package, [weas-widget](https://github.com/superstar54/weas-widget), allows users to visualize and edit with atomistic structures in Jupyter Notebook.

## How to contribute


### Extension
One can build the extension locally by following the steps below.

Download wheels for ase and scikit-image:
```
pip download ase scikit-image --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=manylinux_2_17_x86_64
pip download ase scikit-image --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=win_amd64
pip download ase scikit-image --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=macosx_12_0_arm64
```

Update the wheels in the `batoms/blender_manifest.toml` file.

```console
python scripts/update_wheels.py
```

Build the extension.
```
blender --command extension build --source-dir batoms --output-dir .
```

### Editor
We recommend using [Visual Studio Code](https://code.visualstudio.com/) with the [Blender extension](https://github.com/JacquesLucke/blender_vscode).

In order to develop it inside VSCode, you need at least update the wheel file in the `batoms/blender_manifest.toml` file using the steps above.

### Test
To run the tests, run:

```console
# install pytest
~/apps/blender-4.3.2-linux-x64/4.3/python/bin/python3.11 -m pip install pytest
# run all tests
~/apps/blender-4.3.2-linux-x64/blender -b -P scripts/run_tests.py -- -vv tests
# run a specific test
~/apps/blender-4.3.2-linux-x64/blender -b -P scripts/run_tests.py -- -vv tests/test_batoms.py
```

### Pre-commit
To install the pre-commit hooks, run:

```console
$ pre-commit install
```

## Contact
* Xing Wang  <xingwang1991@gmail.com>
