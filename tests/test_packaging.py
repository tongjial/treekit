import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "packaging" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build = load("build")


def touch(path, text="x"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class TestCleanPreviousOutput:

    def test_removes_the_earlier_output_of_the_build(self, tmp_path):
        dist = tmp_path / "dist"
        release = dist / "release"
        old = [touch(release / "treekit-v0.0.9-linux-x86_64.tar.gz"),
               touch(release / "treekit-v0.0.9-linux-x86_64.tar.gz.sha256"),
               touch(release / "treekit-v0.0.9-windows-x86_64.zip"),
               touch(dist / "treekit"), touch(dist / "treekit.exe")]
        build.clean_previous_output(dist, release)
        assert not any(path.exists() for path in old)

    def test_keeps_everything_else_in_dist(self, tmp_path):
        dist = tmp_path / "dist"
        release = dist / "release"
        keep = [touch(dist / "other_package-1.0-py3-none-any.whl"),
                touch(release / "other-project-1.0.tar.gz"),
                touch(release / "notes.txt"),
                touch(release / "nested" / "treekit-v1-linux-x86_64.tar.gz"),
                touch(dist / "sub" / "treekit")]
        build.clean_previous_output(dist, release)
        assert all(path.exists() for path in keep)

    def test_a_missing_folder_is_fine(self, tmp_path):
        build.clean_previous_output(tmp_path / "dist", tmp_path / "dist" / "release")

    @pytest.mark.skipif(sys.platform == "win32", reason="symbolic links need privileges on Windows")
    def test_a_linked_dist_does_not_delete_what_belongs_to_others(self, tmp_path):
        # dist/ points at a shared folder that has a release/ folder of another project
        shared = tmp_path / "shared"
        keep = [touch(shared / "release" / "other-project-1.0.tar.gz"),
                touch(shared / "release" / "other-project-1.0.tar.gz.sha256"),
                touch(shared / "release" / "deep" / "file.bin"),
                touch(shared / "unrelated.whl")]
        old = touch(shared / "release" / "treekit-v0.0.9-macos-arm64.tar.gz")
        dist = tmp_path / "dist"
        os.symlink(shared, dist, target_is_directory=True)
        build.clean_previous_output(dist, dist / "release")
        assert all(path.exists() for path in keep)
        assert not old.exists()                      # an old treekit archive is still removed
        assert dist.is_symlink() and (shared / "release").is_dir()

    @pytest.mark.skipif(sys.platform == "win32", reason="symbolic links need privileges on Windows")
    def test_a_linked_release_folder_is_not_followed(self, tmp_path):
        shared = tmp_path / "shared"
        keep = touch(shared / "other.txt")
        dist = tmp_path / "dist"
        dist.mkdir()
        os.symlink(shared, dist / "release", target_is_directory=True)
        build.clean_previous_output(dist, dist / "release")
        assert keep.exists() and (dist / "release").is_symlink()
