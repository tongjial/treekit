import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TREEKIT_SCRIPT = REPO / "treekit.py"

# TREEKIT_BIN=path/to/treekit runs the tests that start treekit as a subprocess
# against that executable (the PyInstaller build) instead of against treekit.py
TREEKIT_BIN = os.environ.get("TREEKIT_BIN")
TREEKIT_CMD = [TREEKIT_BIN] if TREEKIT_BIN else [sys.executable, str(TREEKIT_SCRIPT)]

# make `import treekit` work when pytest is run from any directory
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


@pytest.fixture
def run(tmp_path):
    """Run the treekit script as a subprocess; returns CompletedProcess."""

    def _run(*args, cwd=None, timeout=60, env=None):
        # env: variables added to the environment of the run
        return subprocess.run(
            [*TREEKIT_CMD, *args],
            cwd=str(cwd or tmp_path),
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=timeout,
            env={**os.environ, **env} if env else None,
        )

    return _run


@pytest.fixture
def write(tmp_path):
    """Write a text file in tmp_path and return its path."""

    def _write(name, text):
        path = tmp_path / name
        path.write_text(text)
        return path

    return _write
