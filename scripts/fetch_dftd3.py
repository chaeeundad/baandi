#!/usr/bin/env python3 -I
"""Fetch the prebuilt simple-dftd3 1.7.0 shared library and its runtime deps.

Ticket D3-DISPERSION, ADR 0048. BAANDI loads the official simple-dftd3 library
(LGPL-3.0-or-later) at run time through its own MIT adapter
(`crates/dft-dispersion`, Cargo feature `dftd3`). This script downloads the
conda-forge binary packages pinned below (exact file names and SHA-256), verifies
every download, and extracts the shared libraries and the `s-dftd3` command-line
program (the independent oracle of the acceptance tests) into a directory OUTSIDE
the repository:

    <dest>/lib/libs-dftd3.dylib  (and the dependency dylibs)
    <dest>/bin/s-dftd3
    <dest>/manifest.json         (URLs, SHA-256 of each package and each extracted file)
    <dest>/licenses/<package>/   (licence files shipped in the packages)

Default destination: /Volumes/990-2TB/Develop/rust-dft-wt/vendor/dftd3. Binaries
are never committed. The binaries are not modified: their install names use
`@rpath` with `@loader_path/../lib`, so the libraries resolve each other inside
`<dest>/lib`. If a loader still misses one, set
`DYLD_FALLBACK_LIBRARY_PATH=<dest>/lib`.

conda `.conda` packages are zip files holding zstd-compressed tarballs. Python 3.9
and the system bsdtar have no zstd decoder, so the script first fetches the
pinned conda-forge `zstd` 1.5.2 package in the old `.tar.bz2` format (readable by
the standard library) and calls its `libzstd` through ctypes for decompression.

Only osx-arm64 is pinned (the Mac development host). Linux hosts (dev08) add a
pinned table for linux-64 in the same form before use.

Usage: python3 -I scripts/fetch_dftd3.py [--dest DIR] [--force]
Standard library only.
"""
import argparse
import ctypes
import hashlib
import io
import json
import os
import platform
import shutil
import sys
import tarfile
import urllib.request
import zipfile

BASE = "https://conda.anaconda.org/conda-forge"
DEFAULT_DEST = "/Volumes/990-2TB/Develop/rust-dft-wt/vendor/dftd3"

# Exact conda-forge files, resolved on 2026-11-26 from
# https://api.anaconda.org/package/conda-forge/<name> (depends of simple-dftd3
# 1.7.0 gfortran_h643d65c_0: libgcc, libgfortran, libgfortran5 >=15.3, llvm-openmp,
# mctc-lib 0.7.x, toml-f 0.5.x; mctc-lib adds jonquil 0.3.x). `_openmp_mutex` is a
# metapackage without libraries on osx-arm64 and is not needed.
PINS = {
    "osx-arm64": [
        ("simple-dftd3-1.7.0-gfortran_h643d65c_0.conda",
         "e967c2695ac796d3f119c4b4d95b1fe5ab198272448ee7efd95efd8b2c1cdf0a"),
        ("mctc-lib-0.7.0-gfortran_h643d65c_0.conda",
         "eebc93c99cfc29a1ab226bef1bf7e911e8094459ecdf4ecd0ceed90cf5870ba9"),
        ("toml-f-0.5.2-gfortran_h04795c8_0.conda",
         "6bf87f425e8aa8e2e429483576ebe4b9355555af4595835f0dd2e8ea7350ea6d"),
        ("jonquil-0.3.2-gfortran_h04795c8_0.conda",
         "0a7da05cf55fcfbec807f555a09f1c85d3bbf4a6f8bfc7b822734f853d7af846"),
        ("libgfortran-15.3.0-h07b0088_7.conda",
         "f1a85f76172612d28e9fe429be1f3fcb030a23a4e2337016b6da9a5000daa0fa"),
        ("libgfortran5-15.3.0-hfb40ef1_7.conda",
         "9353ac36b76089ec410a14365041e608b1cc1ae0c1d4a66c24e57ff701966e33"),
        ("libgcc-15.3.0-h99d2e81_7.conda",
         "ed98b42a9231d07381f0b27da90c979b2274932a9935deee76423528a15b6b76"),
        ("llvm-openmp-23.1.6-hdb3d66b_0.conda",
         "b777843f1ced814ea99adf79ece75947ef70ed1570e7c3bf5b5633133a496b9b"),
    ],
}
ZSTD_BOOTSTRAP = {
    "osx-arm64": ("zstd-1.5.2-h8128057_4.tar.bz2",
                  "664f7e6b46cfff6022504eafdef26b4ac0f388fa39d682e64da9a7dbfe33fb56",
                  "lib/libzstd.1.5.2.dylib"),
}


def host_subdir():
    sysname, machine = platform.system(), platform.machine()
    if sysname == "Darwin" and machine == "arm64":
        return "osx-arm64"
    if sysname == "Linux" and machine == "x86_64":
        return "linux-64"
    return f"{sysname}-{machine}"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def download(subdir, name, expected, cache):
    path = os.path.join(cache, name)
    if os.path.exists(path):
        with open(path, "rb") as f:
            data = f.read()
        if sha256(data) == expected:
            return data
    url = f"{BASE}/{subdir}/{name}"
    print(f"download {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=120) as r:
        data = r.read()
    got = sha256(data)
    if got != expected:
        raise SystemExit(f"SHA-256 mismatch for {name}: expected {expected}, got {got}")
    with open(path, "wb") as f:
        f.write(data)
    return data


class Zstd:
    """Streaming zstd decompression through libzstd (ctypes, public zstd.h API)."""

    class Buf(ctypes.Structure):
        _fields_ = [("src", ctypes.c_void_p), ("size", ctypes.c_size_t), ("pos", ctypes.c_size_t)]

    def __init__(self, lib_path):
        lib = ctypes.CDLL(lib_path)
        lib.ZSTD_createDStream.restype = ctypes.c_void_p
        lib.ZSTD_freeDStream.argtypes = [ctypes.c_void_p]
        lib.ZSTD_initDStream.argtypes = [ctypes.c_void_p]
        lib.ZSTD_initDStream.restype = ctypes.c_size_t
        lib.ZSTD_decompressStream.argtypes = [ctypes.c_void_p, ctypes.POINTER(self.Buf), ctypes.POINTER(self.Buf)]
        lib.ZSTD_decompressStream.restype = ctypes.c_size_t
        lib.ZSTD_isError.argtypes = [ctypes.c_size_t]
        lib.ZSTD_isError.restype = ctypes.c_uint
        lib.ZSTD_getErrorName.argtypes = [ctypes.c_size_t]
        lib.ZSTD_getErrorName.restype = ctypes.c_char_p
        self.lib = lib

    def decompress(self, data):
        lib = self.lib
        ds = lib.ZSTD_createDStream()
        if not ds:
            raise SystemExit("ZSTD_createDStream failed")
        try:
            lib.ZSTD_initDStream(ds)
            src = ctypes.create_string_buffer(data, len(data))
            inb = self.Buf(ctypes.cast(src, ctypes.c_void_p), len(data), 0)
            chunk = 1 << 20
            dst = ctypes.create_string_buffer(chunk)
            out = io.BytesIO()
            while True:
                outb = self.Buf(ctypes.cast(dst, ctypes.c_void_p), chunk, 0)
                ret = lib.ZSTD_decompressStream(ds, ctypes.byref(outb), ctypes.byref(inb))
                if lib.ZSTD_isError(ret):
                    raise SystemExit("zstd: " + lib.ZSTD_getErrorName(ret).decode())
                out.write(dst.raw[: outb.pos])
                if inb.pos >= inb.size and outb.pos < chunk:
                    break
            return out.getvalue()
        finally:
            lib.ZSTD_freeDStream(ds)


def safe_members(tar):
    for m in tar.getmembers():
        name = m.name
        if name.startswith("/") or ".." in name.split("/"):
            raise SystemExit(f"unsafe path in package: {name}")
        if m.issym() or m.islnk():
            target = m.linkname
            if target.startswith("/") or ".." in target.split("/"):
                raise SystemExit(f"unsafe link in package: {name} -> {target}")
        yield m


def keep(name):
    return (name.startswith("lib/") and (".dylib" in name or ".so" in name)) \
        or name == "bin/s-dftd3" or name.startswith("include/") \
        or name.startswith("info/licenses/")


def extract_conda(data, zstd, dest, pkg):
    files = []
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        inner = [n for n in z.namelist() if n.startswith(("pkg-", "info-")) and n.endswith(".tar.zst")]
        for n in inner:
            raw = zstd.decompress(z.read(n))
            with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as t:
                files += extract_tar(t, dest, pkg)
    return files


def extract_tar(t, dest, pkg):
    files = []
    for m in safe_members(t):
        if not keep(m.name):
            continue
        out_name = m.name
        if m.name.startswith("info/licenses/"):
            out_name = os.path.join("licenses", pkg, m.name[len("info/licenses/"):])
        out = os.path.join(dest, out_name)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        if m.issym():
            if os.path.lexists(out):
                os.remove(out)
            os.symlink(m.linkname, out)
            files.append({"path": out_name, "symlink": m.linkname})
        elif m.isfile():
            blob = t.extractfile(m).read()
            with open(out, "wb") as f:
                f.write(blob)
            os.chmod(out, m.mode & 0o755 | 0o400)
            files.append({"path": out_name, "sha256": sha256(blob)})
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dest", default=DEFAULT_DEST)
    ap.add_argument("--force", action="store_true", help="re-extract even if a manifest exists")
    args = ap.parse_args()
    subdir = host_subdir()
    if subdir not in PINS:
        raise SystemExit(f"no pinned packages for {subdir}; add a PINS table first")
    dest = os.path.abspath(args.dest)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if os.path.commonpath([dest, repo]) == repo:
        raise SystemExit("destination is inside the repository; binaries are never committed")
    manifest_path = os.path.join(dest, "manifest.json")
    if os.path.exists(manifest_path) and not args.force:
        print(f"{manifest_path} exists; use --force to re-extract", file=sys.stderr)
        print(os.path.join(dest, "lib"))
        return 0
    cache = os.path.join(dest, "pkgs")
    os.makedirs(cache, exist_ok=True)

    zname, zsha, zlib_rel = ZSTD_BOOTSTRAP[subdir]
    zdata = download(subdir, zname, zsha, cache)
    tool_dir = os.path.join(dest, "_bootstrap")
    os.makedirs(tool_dir, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(zdata), mode="r:bz2") as t:
        member = t.getmember(zlib_rel)
        blob = t.extractfile(member).read()
    zlib_path = os.path.join(tool_dir, os.path.basename(zlib_rel))
    with open(zlib_path, "wb") as f:
        f.write(blob)
    zstd = Zstd(zlib_path)

    manifest = {"platform": subdir, "base_url": BASE, "packages": [],
                "bootstrap": {"name": zname, "url": f"{BASE}/{subdir}/{zname}", "sha256": zsha,
                              "use": "libzstd for .conda decompression only"}}
    for name, expected in PINS[subdir]:
        data = download(subdir, name, expected, cache)
        pkg = name.rsplit("-", 2)[0]
        files = extract_conda(data, zstd, dest, pkg)
        manifest["packages"].append({"name": name, "url": f"{BASE}/{subdir}/{name}",
                                     "sha256": expected, "files": files})
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=1)
    shutil.rmtree(tool_dir, ignore_errors=True)
    lib = os.path.join(dest, "lib")
    print(f"BAANDI_DFTD3_LIB={os.path.join(lib, 'libs-dftd3.dylib' if subdir.startswith('osx') else 'libs-dftd3.so')}",
          file=sys.stderr)
    print(lib)
    return 0


if __name__ == "__main__":
    sys.exit(main())
