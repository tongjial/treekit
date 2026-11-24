#!/usr/bin/env python3
"""Build the single-file treekit executable with PyInstaller and pack it for a release.

    python packaging/build.py [--max-glibc 2.28]

Needs PyInstaller and DendroPy in the Python that runs it (packaging/requirements-build.txt).
PyInstaller cannot cross-compile: this builds for the system it runs on. It writes, in
dist/release/, an archive named treekit-v<version>-<os>-<arch> (.zip on Windows, .tar.gz
elsewhere) and a .sha256 file next to it. The archive holds the executable, LICENSE, README.md
and THIRD_PARTY_NOTICES.txt.
"""
import argparse
import hashlib
import importlib.metadata
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
WORK = ROOT / "build" / "pyinstaller"
RELEASE = DIST / "release"

# shared libraries that belong to glibc itself: they are never bundled
GLIBC_LIBS = ("linux-vdso", "ld-linux", "libc.so", "libm.so", "libdl.so", "libpthread.so",
              "librt.so", "libutil.so", "libresolv.so", "libnsl.so", "libcrypt.so")


def read_version():
    text = (ROOT / "treekit.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__ = "([^"]+)"', text, re.MULTILINE)
    if not match:
        sys.exit("error: __version__ not found in treekit.py")
    return match.group(1)


def platform_tag():
    system = {"Linux": "linux", "Darwin": "macos", "Windows": "windows"}.get(platform.system())
    if system is None:
        sys.exit("error: unsupported system: " + platform.system())
    machine = platform.machine().lower()
    arch = {"amd64": "x86_64", "x86_64": "x86_64", "arm64": "arm64", "aarch64": "arm64"}.get(
        machine, machine)
    return system, arch


def build_executable(exe_name):
    command = [
        sys.executable, "-m", "PyInstaller",
        "--onefile",
        "--name", "treekit",
        "--noconfirm", "--clean",
        "--distpath", str(DIST),
        "--workpath", str(WORK),
        "--specpath", str(WORK),
        "--log-level", "WARN",
        str(ROOT / "treekit.py"),
    ]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, cwd=str(ROOT))
    exe = DIST / exe_name
    if not exe.is_file():
        sys.exit("error: PyInstaller did not produce " + str(exe))
    return exe


def glibc_versions(path):
    # the GLIBC_x.y symbol versions that a shared library or executable needs
    result = subprocess.run(["objdump", "-T", str(path)], capture_output=True, text=True)
    found = set(re.findall(r"GLIBC_(\d+(?:\.\d+)*)", result.stdout))
    return {tuple(int(part) for part in version.split(".")) for version in found}


def check_glibc(exe, limit):
    # The executable runs on every Linux whose glibc is at least as new as the newest
    # symbol it needs. Checked here: the bootloader, Python's own shared libraries (they
    # are bundled) and the non-glibc system libraries those depend on (also bundled).
    limit_tuple = tuple(int(part) for part in limit.split("."))
    prefix = Path(sys.base_prefix)
    files = [exe]
    files += sorted(prefix.glob("lib/libpython3*.so*"))
    files += sorted(prefix.glob("lib/python3*/lib-dynload/*.so"))
    for path in list(files[1:]):
        ldd = subprocess.run(["ldd", str(path)], capture_output=True, text=True).stdout
        for line in ldd.splitlines():
            match = re.match(r"\s*(\S+) => (/\S+)", line)
            if match and not any(match.group(1).startswith(name) for name in GLIBC_LIBS):
                library = Path(match.group(2))
                if library not in files:
                    files.append(library)

    newest = (0,)
    newest_file = None
    for path in files:
        versions = glibc_versions(path)
        if versions and max(versions) > newest:
            newest, newest_file = max(versions), path
    text = ".".join(map(str, newest))
    print("newest glibc symbol needed: %s (%s); limit %s; %d files checked"
          % (text, newest_file.name if newest_file else "-", limit, len(files)))
    if newest > limit_tuple:
        sys.exit("error: %s needs glibc %s, newer than %s: build with an older Python"
                 % (newest_file, text, limit))


def third_party_notices():
    # the DendroPy license has to travel with the executable that contains DendroPy
    dist = importlib.metadata.distribution("dendropy")
    licenses = [f for f in (dist.files or []) if re.match(r"(LICENSE|COPYING)", f.name, re.I)]
    if not licenses:
        sys.exit("error: license file of dendropy not found in its metadata")
    dendropy_license = licenses[0].read_text(encoding="utf-8").strip()
    pyinstaller = importlib.metadata.version("pyinstaller")
    return (
        "treekit is distributed under the GNU General Public License v3 (see LICENSE).\n"
        "This executable also contains the following software.\n\n"
        "Python %s\n"
        "  Python Software Foundation License: https://docs.python.org/3/license.html\n\n"
        "PyInstaller %s (bootloader)\n"
        "  GPL-2.0-or-later with an exception that allows distributing the executables it builds\n"
        "  under any license: https://github.com/pyinstaller/pyinstaller/blob/develop/COPYING.txt\n\n"
        "DendroPy %s\n"
        "  BSD 3-Clause License. License text:\n\n%s\n"
        % (platform.python_version(), pyinstaller, dist.version, dendropy_license)
    )


def make_archive(exe, name, is_windows):
    RELEASE.mkdir(parents=True, exist_ok=True)
    notices = WORK / "THIRD_PARTY_NOTICES.txt"
    notices.write_text(third_party_notices(), encoding="utf-8")
    members = [(exe, exe.name, 0o755), (ROOT / "LICENSE", "LICENSE", 0o644),
               (ROOT / "README.md", "README.md", 0o644), (notices, notices.name, 0o644)]

    if is_windows:
        archive = RELEASE / (name + ".zip")
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
            for source, arcname, _ in members:
                zf.write(source, name + "/" + arcname)
    else:
        archive = RELEASE / (name + ".tar.gz")
        with tarfile.open(archive, "w:gz") as tf:
            for source, arcname, mode in members:
                info = tf.gettarinfo(str(source), name + "/" + arcname)
                info.mode = mode            # keep the executable bit, whatever the umask
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                with open(source, "rb") as handle:
                    tf.addfile(info, handle)

    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    Path(str(archive) + ".sha256").write_text("%s  %s\n" % (digest, archive.name))
    return archive, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--max-glibc", default="2.28",
                        help="Linux only: fail if the executable needs a newer glibc. Default: 2.28")
    args = parser.parse_args()

    version = read_version()
    system, arch = platform_tag()
    is_windows = system == "windows"
    name = "treekit-v%s-%s-%s" % (version, system, arch)
    print("building %s with Python %s" % (name, platform.python_version()), flush=True)

    shutil.rmtree(DIST, ignore_errors=True)
    exe = build_executable("treekit.exe" if is_windows else "treekit")

    # the new executable must run at all, and report the version of the source
    shown = subprocess.run([str(exe), "--version"], capture_output=True, text=True, timeout=120)
    if shown.returncode != 0 or shown.stdout.strip() != "treekit " + version:
        sys.exit("error: '%s --version' printed %r (exit %d), expected 'treekit %s'"
                 % (exe.name, shown.stdout + shown.stderr, shown.returncode, version))

    if system == "linux":
        check_glibc(exe, args.max_glibc)

    archive, digest = make_archive(exe, name, is_windows)
    print("built %s (%.1f MB)\nsha256 %s" % (archive, archive.stat().st_size / 1e6, digest))


if __name__ == "__main__":
    main()
