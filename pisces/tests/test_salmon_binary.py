import os
import sys
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from pisces import find_salmon_binary, find_data_directory

def test_find_salmon_binary_system(monkeypatch):
    # Simulate salmon in PATH
    monkeypatch.setenv('PATH', os.environ['PATH'])
    # If salmon is not installed, this will raise, so we skip if not present
    try:
        path = find_salmon_binary(find_data_directory())
        assert os.path.basename(path) == 'salmon'
    except FileNotFoundError:
        pytest.skip('No salmon binary in PATH or bundled')

def test_find_salmon_binary_bundled(tmp_path, monkeypatch):
    # Simulate no salmon in PATH, but bundled exists
    monkeypatch.setenv('PATH', '')
    data_dir = tmp_path
    salmon_dir = data_dir / 'redist' / 'salmon' / 'bin'
    salmon_dir.mkdir(parents=True)
    salmon_bin = salmon_dir / 'salmon'
    salmon_bin.write_text('#!/bin/bash\necho salmon')
    os.chmod(salmon_bin, 0o755)
    path = find_salmon_binary(str(data_dir))
    assert str(salmon_bin) == path

def test_find_salmon_binary_not_found(tmp_path, monkeypatch):
    # Simulate no salmon anywhere
    monkeypatch.setenv('PATH', '')
    with pytest.raises(FileNotFoundError):
        find_salmon_binary(str(tmp_path))
