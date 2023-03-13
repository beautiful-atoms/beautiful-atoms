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

# def test_git_process(monkeypatch):
#     import install
#     from install import _gitclone, _gitcheckout
    

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
    
    
