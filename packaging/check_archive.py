#!/usr/bin/env python3
"""Check the release archive written by packaging/build.py, as a user would use it.

    python packaging/check_archive.py [dist/release]

Verifies the .sha256 file, unpacks the archive into a temporary directory, checks its files,
and runs the unpacked executable: --version, a summary, a parallel summary (-c 2: this fails
if multiprocessing is not frozen correctly) and a tree file with non-ASCII taxon names.
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
NAME = re.compile(r"treekit-v(?P<version>.+)-(?P<system>linux|macos|windows)-(?P<arch>[^-]+)$")


def fail(message):
    sys.exit("error: " + message)


def run(command, **kwargs):
    return subprocess.run([str(part) for part in command], capture_output=True, timeout=300, **kwargs)


def main():
    release = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "dist" / "release")
    archives = sorted(p for p in release.iterdir() if p.name.endswith((".zip", ".tar.gz")))
    if len(archives) != 1:
        fail("expected one archive in %s, found %s" % (release, [p.name for p in archives]))
    archive = archives[0]
    stem = archive.name[:-len(".zip")] if archive.name.endswith(".zip") else archive.name[:-len(".tar.gz")]
    match = NAME.fullmatch(stem)
    if not match:
        fail("unexpected archive name: " + archive.name)
    is_windows = match["system"] == "windows"
    if is_windows != archive.name.endswith(".zip"):
        fail("Windows archives are .zip, the others .tar.gz: " + archive.name)

    # the checksum file is what users compare with
    checksum_file = Path(str(archive) + ".sha256")
    recorded = checksum_file.read_text(encoding="utf-8").split()
    actual = hashlib.sha256(archive.read_bytes()).hexdigest()
    if recorded != [actual, archive.name]:
        fail("%s does not match the archive: %s" % (checksum_file.name, recorded))

    with tempfile.TemporaryDirectory() as tmp:
        shutil.unpack_archive(str(archive), tmp)
        folder = Path(tmp) / stem
        exe = folder / ("treekit.exe" if is_windows else "treekit")
        expected = {exe.name, "LICENSE", "README.md", "THIRD_PARTY_NOTICES.txt"}
        found = {p.name for p in folder.iterdir()} if folder.is_dir() else set()
        if found != expected:
            fail("the archive holds %s, expected %s" % (sorted(found), sorted(expected)))
        if not is_windows and not os.access(exe, os.X_OK):
            fail("%s is not executable after unpacking" % exe.name)
        notices = (folder / "THIRD_PARTY_NOTICES.txt").read_text(encoding="utf-8")
        for needed in ("DendroPy", "BSD 3-Clause", "PyInstaller", "Python "):
            if needed not in notices:
                fail("THIRD_PARTY_NOTICES.txt does not mention %r" % needed)

        shown = run([exe, "--version"], text=True)
        if shown.returncode != 0 or shown.stdout.strip() != "treekit " + match["version"]:
            fail("--version printed %r, expected 'treekit %s'" % (shown.stdout + shown.stderr, match["version"]))

        single = run([exe, "summary", "-i", EXAMPLES / "test"], text=True, encoding="utf-8")
        if single.returncode != 0 or len(single.stdout.splitlines()) != 2:
            fail("summary failed: %r %r" % (single.stdout, single.stderr))

        # -c 2 starts worker processes by running the executable again
        parallel = run([exe, "summary", "-c", "2", "-i", EXAMPLES / "test",
                        EXAMPLES / "concat_6genes.trees"], text=True, encoding="utf-8")
        if parallel.returncode != 0 or len(parallel.stdout.splitlines()) != 3:
            fail("parallel summary failed: %r %r" % (parallel.stdout, parallel.stderr))

        # non-ASCII names must come out as UTF-8 whatever the console encoding is
        names = Path(tmp) / "names.tre"
        names.write_bytes("((Café_é:1,'Zoë b':1):1,C:1,D:1);\n".encode("utf-8"))
        taxa = run([exe, "taxa", "-i", names])
        lines = taxa.stdout.decode("utf-8", errors="replace").splitlines()
        if taxa.returncode != 0 or lines != ["C", "Café_é", "D", "Zoë b"]:
            fail("taxa printed %r (stderr %r)" % (lines, taxa.stderr))

    print("ok: %s (%.1f MB), checked the checksum, the files, --version, summary, -c 2 and UTF-8 names"
          % (archive.name, archive.stat().st_size / 1e6))


if __name__ == "__main__":
    main()
