import pytest
import sys
from pathlib import Path
curdir = Path(__file__).parent.resolve()
print(curdir)
sys.path.append(curdir.parent.as_posix())
print(sys.path)
import install
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
    
