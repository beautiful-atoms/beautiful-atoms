import bpy
from batoms.batoms import Batoms
from ase.build import molecule

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_polyhedra_species(tio2):
    """
    This is an example to show different polyhedra for the same element, but different species.
    """
    assert len(tio2.polyhedra.settings) == 1
    # We set the first one to a new species Ti_1.
    tio2.replace([0], "Ti_1")
    tio2.bond.settings
    tio2.model_style = 2
    assert len(tio2.polyhedra) == 12
    # remove Ti_1 to polyhedra.
    tio2.polyhedra.settings.remove("Ti_1")
    assert len(tio2.polyhedra.settings) == 1
    tio2.model_style = 2
    assert len(tio2.polyhedra) == 6


def test_polyhedra_molecule():
    bpy.ops.batoms.delete()
    ch4 = molecule("CH4")
    mh4 = ch4.copy()
    mh4.translate([2, 2, 0])
    mh4[0].symbol = "N"
    ch4 = ch4 + mh4
    ch4 = Batoms("ch4", from_ase=ch4)
    ch4.model_style = 2
    ch4.bond.settings[("C", "H")].polyhedra = True
    ch4.model_style = 1


def test_polyhedra_crystal(tio2):
    tio2.model_style = 2
    tio2 = tio2 * [3, 3, 3]
    tio2.pbc = False
    tio2.model_style = 2


def test_polyhedra_setting(ch4):
    ch4.bond.settings[("C", "H")].polyhedra = True
    ch4.model_style = 2
    ch4.pbc = True
    ch4.cell = [3, 3, 3]
    ch4 = ch4 * [2, 1, 1]
    sel1 = ch4.selects.add("sel1", [0, 1, 2, 3, 4])
    sel1.model_style = 1
    ch4.polyhedra.settings.remove("C")
    assert len(ch4.polyhedra.settings) == 1
    ch4.polyhedra.settings.add("C")
    assert len(ch4.polyhedra.settings) == 2
    ch4.polyhedra.settings["C"] = {"color": [0.8, 0.1, 0.3, 0.3]}
    ch4.model_style = 2
    ch4.render.engine = "workbench"
    if use_cycles:
        set_cycles_res(ch4)
    ch4.get_image([1, 1, 0], output="polyhedra.png", **extras)
