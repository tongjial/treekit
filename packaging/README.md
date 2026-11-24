# Packaging and releases

treekit is released as single-file executables built with [PyInstaller](https://pyinstaller.org),
one per platform (PyInstaller cannot cross-compile): Linux x86_64, Windows x86_64, macOS arm64 and
macOS x86_64.

## Making a release
1. Set `__version__` in `treekit.py` (`pyproject.toml` reads it from there), commit and merge.
2. Tag the commit and push the tag: `git tag v0.2.0 && git push origin v0.2.0`.
   A tag with a hyphen (`v0.2.0-rc1`) is published as a pre-release.

`.github/workflows/release.yml` then builds on all four platforms, runs the whole test suite
against each executable (`TREEKIT_BIN`), checks each archive (`check_archive.py`), and only if all
four pass creates the GitHub release with the archives, their `.sha256` files and `SHA256SUMS`.
The workflow stops at once if the tag differs from `__version__`. Pull requests that change the
packaging run the same builds and tests without publishing.

## Building locally
```
uv run --no-project --python 3.12 --with-requirements packaging/requirements-build.txt packaging/build.py
uv run --no-project --python 3.12 python packaging/check_archive.py
TREEKIT_BIN=$PWD/dist/treekit pytest          # the tests against the executable
```
(`export UV_PYTHON_PREFERENCE=only-managed` first, as the workflow does.) The result is in `dist/release/`.

## Things to know
* **Python from uv.** On Linux, the glibc that the executable needs is that of the Python it was built
  with: Ubuntu's own Python 3.12 needs glibc 2.38, the uv one 2.17. `build.py` fails if the executable
  needs a glibc newer than `--max-glibc` (2.28).
* **`multiprocessing.freeze_support()`** at the start of `main()` is required: `-c` starts the workers
  by running the executable again, and without it they start as a new treekit and the run hangs.
* **Pinned versions.** `requirements-build.txt` pins dendropy and PyInstaller so that a tag always
  gives the same executable; change them on purpose, and keep dendropy in line with `pyproject.toml`.
* **Licenses.** treekit is GPL-3.0; the archives carry `LICENSE` and `THIRD_PARTY_NOTICES.txt`
  (Python, PyInstaller's bootloader, and the DendroPy license text, which must travel with DendroPy).
* **Not signed or notarized.** Signing needs an Apple developer account (macOS) or a code-signing
  certificate (Windows); until then users meet Gatekeeper and SmartScreen warnings (see the README).
