import pytest
import sys
from pathlib import Path
curdir = Path(__file__).parent.resolve()
print(curdir)
sys.path.append(curdir.parent.as_posix())
print(sys.path)
import os
from tempfile import TemporaryDirectory

def test_empty_dir():
    from install import _is_empty_dir
    assert _is_empty_dir("/") is False
    with TemporaryDirectory() as tmpdir1:
        assert _is_empty_dir(tmpdir1) is True
        tmpdir2 = Path(tmpdir1) / "test"
        os.makedirs(tmpdir2, exist_ok=True)
        assert _is_empty_dir(tmpdir1) is False
        assert _is_empty_dir(tmpdir2) is True
    with pytest.raises(FileNotFoundError):
        _is_empty_dir(tmpdir1)
        _is_empty_dir(tmpdir2)
    

def test_default_location_linux(fs):
    """Test requires the `pyfakefs` package's fs fixture
    """
    from pathlib import Path
    from install import _get_default_locations
    with pytest.raises(NotImplementedError):
        _get_default_locations("linux")

def test_default_location_macos(fs, monkeypatch):
    """Test requires the `pyfakefs` package's fs fixture
    """
    monkeypatch.setattr('builtins.input', lambda _: "0")
    from pathlib import Path
    import os
    from install import _get_default_locations
    # macos case 1
    fdn1 = "/Applications/Blender.app/Contents/Resources/3.4"
    fdd1 = fs.create_dir(fdn1)
    assert _get_default_locations("macos") == Path(fdn1)
    

    # macos case 2: Blender 3.4 and 3.1 both exist, should return 3.4
    fdn2 = "/Applications/Blender.app/Contents/Resources/3.1"
    fdd2 = fs.create_dir(fdn2)
    assert _get_default_locations("macos").name == "3.4"

    fs.rmdir(fdn1)
    with pytest.raises(FileNotFoundError):
        _get_default_locations("macos")
    fs.rmdir(fdn2)

    # macos case 3: both user and application blender
    fdn3 = os.path.expandvars("$HOME/Applications/Blender.app/Contents/Resources/3.4")
    fdd1 = fs.create_dir(fdn1)
    fdd3 = fs.create_dir(fdn3)
    assert _get_default_locations("macos") == Path(fdn1)

def test_blender_bin(fs):
    # Windows
    from pathlib import Path
    from install import _get_blender_bin
    bundle_root = "/c/Program Files/Blender Foundation/Blender 3.4/3.4"
    fake_bin = "/c/Program Files/Blender Foundation/Blender 3.4/blender.exe"
    fs.create_dir(bundle_root)
    fs.create_file(fake_bin)
    assert _get_blender_bin("windows", bundle_root) == Path(fake_bin).resolve()
    # User should not point to the upper level
    with pytest.raises(FileNotFoundError):
        _get_blender_bin("windows", Path(bundle_root).parent)
    fs.rmdir(bundle_root)

    # Linux
    fs.create_dir("./blender/3.4")
    fs.create_file("./blender/blender")
    assert _get_blender_bin("linux", "./blender/3.4") == Path("./blender/blender").resolve()
    # fs.remove("./blender/blender")
    # fs.rmdir("./blender")

    #Macos
    fs.create_dir("/Applications/Blender.app/Contents/Resources/3.4")
    fs.create_file("/Applications/Blender.app/Contents/MacOS/Blender")
    assert _get_blender_bin("macos", "/Applications/Blender.app/Contents/Resources/3.4")


def test_blender_py_location(monkeypatch):
    def fake_output(*args, **argv):
        output = ("Blender 3.2.2 (hash bcfdb14560e7 built 2022-08-02 23:31:44)\n"
                  "Read prefs: /Users/ttian/Library/Application Support/Blender/3.2/config/userpref.blend\n"
                  "Python binary:  /Applications/Blender.app/Contents/Resources/3.2/python/bin/python3.10\n"
                  "\n"
                  "Blender quit\n")
        class MockClass():
            stdout = output.encode("utf8")
        return MockClass()
    from pathlib import Path
    import install
    monkeypatch.setattr("install._run_blender_multiline_expr", fake_output)
    from install import _get_blender_py
    assert _get_blender_py("test_bin") == Path("/Applications/Blender.app/Contents/Resources/3.2/python/bin/python3.10")


def test_factory_version(monkeypatch):
    def fake_output(*args, **argv):
        output = ("Blender 3.4.1 (hash 55485cb379f7 built 2022-12-20 00:46:45)\n"
                  "Read prefs: /home/pink/.config/blender/3.4/config/userpref.blend\n"
                  "Python Version:  3.10.8 (main, Oct 24 2022, 20:47:11) [GCC 9.3.1 20200408 (Red Hat 9.3.1-2)]\n"
                  "Numpy Version:  1.24.2\n"
                  "\n"
                  "Blender quit\n")
        class MockClass():
            stdout = output.encode("utf8")
        return MockClass()
    import install
    monkeypatch.setattr("install._run_blender_multiline_expr", fake_output)
    from install import _get_factory_versions
    py_version, numpy_version = _get_factory_versions("test_bin")
    assert py_version == "3.10.8"
    assert numpy_version == "1.24.2"

def test_blender_version(monkeypatch):
    def fake_output(*args, **argv):
        output = ("Blender 3.4.1 (hash 55485cb379f7 built 2022-12-20 00:46:45)\n"
                  "Read prefs: /home/pink/.config/blender/3.4/config/userpref.blend\n"
                  "\n"
                  "Blender quit")
        class MockClass():
            stdout = output.encode("utf8")
        return MockClass()
    import install
    monkeypatch.setattr("install._run_process", fake_output)
    from install import _get_blender_version
    assert _get_blender_version("test_bin") == "3.4.1"


def test_rename_dir(fs):
    from install import _rename_dir
    import os
    # case 1: no src
    with pytest.raises(FileNotFoundError):
        _rename_dir("_empty", "_empty")
    # case 2: has src but no dst
    fs.create_dir("source")
    _rename_dir("source", "target")
    assert os.path.isdir("target")
    # case 3: target exists
    fs.create_dir("source1/s1")
    _rename_dir("source1", "target")
    assert os.path.isdir("target")
    assert os.path.isdir("target/s1")
    assert not os.path.isdir("source1")
    # case 4: non-empty target, raise error
    fs.create_dir("source2/s2")
    with pytest.raises(OSError):
        _rename_dir("source2", "target")

def test_symlink(fs):
    from install import _symlink_dir
    import os
    # case 1: no src
    with pytest.raises(FileNotFoundError):
        _symlink_dir("_empty", "_empty")
    # case 2: has src but no dst
    fs.create_dir("source")
    _symlink_dir("source", "target")
    assert os.path.islink("target")
    # case 3: target exists
    fs.create_dir("source1/s1")
    _symlink_dir("source1", "target")
    assert os.path.islink("target")
    assert os.path.isdir("target/s1")
    assert os.path.isdir("source1")
    # case 4: non-empty linked target, will work
    fs.create_dir("source2/s2")
    _symlink_dir("source2", "target")
    assert os.path.isdir("target/s2")
    # case 5: non-empty standard target, will not work
    fs.create_dir("target1/t1")
    with pytest.raises(OSError):
        _symlink_dir("source2", "target1")



def test_binary_file(fs):
    from install import _is_binary_file
    import sys
    import struct
    ff = fs.create_file("test_plain")
    with open("test_plain", "w") as fd:
        fd.write("sdf")
    assert _is_binary_file("test_plain") is False
    fs.create_file("test_bin")
    with open('test_bin', 'wb') as f:
        packed_value = struct.pack('i', 1111)
        f.write(packed_value)
    assert _is_binary_file("test_bin")


def test_conda_name_recognition():
    from install import _is_conda_name_abbrev
    assert _is_conda_name_abbrev("base")
    assert _is_conda_name_abbrev("my-test")
    assert not _is_conda_name_abbrev("./local-env")
    assert not _is_conda_name_abbrev("C:\\local-env")
