import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TREEKIT_SCRIPT = REPO / "treekit.py"

# make `import treekit` work when pytest is run from any directory
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


@pytest.fixture
def run(tmp_path):
    """Run the treekit script as a subprocess; returns CompletedProcess."""

    def _run(*args, cwd=None, timeout=60):
        return subprocess.run(
            [sys.executable, str(TREEKIT_SCRIPT), *args],
            cwd=str(cwd or tmp_path),
            capture_output=True,
            text=True,
            timeout=timeout,
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
