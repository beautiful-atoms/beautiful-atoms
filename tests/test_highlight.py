try:
    from _common_helpers import has_display

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_highlight(ch4):
    """ """
    ch4.selects.add("s1", [2])
    ch4.highlight.settings["s1"] = {
        "select": "s1",
        "scale": 0.5,
        "color": (0.5, 0.5, 0, 0.4),
    }
    ch4.highlight.draw()
    assert len(ch4.highlight.obj.data.vertices) == 1
    ch4.get_image(output="highlight_ch4.png")
    # area = ch4.molecular_surface.get_psasa()


def test_highlight_style(ch4):
    """ """
    ch4.selects.add("s1", [2])
    ch4.highlight.settings["s1"] = {
        "select": "s1",
        "scale": 0.8,
        "color": (0.5, 0.5, 0, 0.4),
        "style": "1",
    }
    ch4.highlight.draw()
    assert len(ch4.highlight.obj.data.vertices) == 1
    ch4.get_image(output="highlight_ch4_cube.png")
    # area = ch4.molecular_surface.get_psasa()


def test_highlight_multi(ch4):
    """ """
    ch4.selects.add("s1", [2])
    ch4.selects.add("s2", [3])
    ch4.highlight.settings["s1"] = {
        "select": "s1",
        "scale": 0.5,
        "color": (0.5, 0.5, 0, 0.4),
    }
    ch4.highlight.settings["s2"] = {
        "select": "s2",
        "scale": 0.8,
        "color": (0, 0.5, 0.5, 0.4),
    }
    ch4.highlight.draw()
    assert len(ch4.highlight.obj.data.vertices) == 2
    ch4.get_image(output="highlight_ch4_multi.png")
    # area = ch4.molecular_surface.get_psasa()
