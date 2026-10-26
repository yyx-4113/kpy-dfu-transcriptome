#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_deposit.py -- build the Zenodo deposit archive for the Kpy x DFU
reproducibility package and split it into ~1.5 MiB sequential volumes so each
PUT stays under the flaky egress-proxy ceiling.

Output (written next to this script, i.e. inside 07_复现包):
  kpy_dfu_deposit.zip            -- full archive (also kept for single-file upload)
  kpy_dfu_deposit.zip.001 ...    -- sequential split volumes
  00_REASSEMBLE.txt              -- reassembly instructions

Excludes (regenerable / secrets / transient):
  .git, zenodo_token.txt, the split volumes themselves, 00_REASSEMBLE.txt,
  raw GEO counts (01_数据/原始数据, 01_数据/处理后数据, 01_数据/MR/*.gz, 01_数据/工具),
  *.npz, *.h5ad, __pycache__, .venv, venv.
"""
from __future__ import annotations

import zipfile
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHUNK = 1572864  # 1.5 MiB
ZIP_NAME = "kpy_dfu_deposit.zip"
REASM = HERE / "00_REASSEMBLE.txt"

EXCLUDE_DIRS = {
    ".git",
    "01_数据/原始数据",
    "01_数据/处理后数据",
    "01_数据/工具",
    "__pycache__",
    ".venv",
    "venv",
}
EXCLUDE_FILES = {
    "zenodo_token.txt",
    "00_REASSEMBLE.txt",
    ZIP_NAME,
}
EXCLUDE_SUFFIXES = {".npz", ".h5ad", ".pyc"}


def _excluded(rel: str) -> bool:
    parts = rel.split("/")
    if any(p in EXCLUDE_DIRS for p in parts):
        return True
    name = parts[-1]
    if name in EXCLUDE_FILES:
        return True
    if name.endswith(".gz") and rel.startswith("01_数据/MR"):
        return True
    if any(name.endswith(s) for s in EXCLUDE_SUFFIXES):
        return True
    return False


def build_zip() -> Path:
    zip_path = HERE / ZIP_NAME
    n = 0
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(HERE.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(HERE).as_posix()
            if _excluded(rel):
                continue
            zf.write(p, rel)
            n += 1
    print(f"[zip] {zip_path.name}: {n} files, {zip_path.stat().st_size:,} bytes")
    return zip_path


def split_zip(zip_path: Path):
    vols = []
    idx = 1
    with zip_path.open("rb") as src:
        while True:
            chunk = src.read(CHUNK)
            if not chunk:
                break
            vol = HERE / f"{ZIP_NAME}.{idx:03d}"
            vol.write_bytes(chunk)
            vols.append(vol)
            idx += 1
    print(f"[split] {len(vols)} volumes of {CHUNK:,} bytes each")
    reasm = (
        "Reassemble kpy_dfu_deposit.zip from its sequential volumes before unzipping.\n\n"
        "Linux / macOS / WSL:\n"
        "    cat kpy_dfu_deposit.zip.* > kpy_dfu_deposit.zip\n\n"
        "PowerShell:\n"
        "    Get-Content -Raw kpy_dfu_deposit.zip.[0-9][0-9][0-9] | "
        "Set-Content -NoNewline kpy_dfu_deposit.zip\n"
        "  (or:  cmd /c copy /b kpy_dfu_deposit.zip.001+kpy_dfu_deposit.zip.002+... "
        "kpy_dfu_deposit.zip)\n\n"
        "Then:  unzip kpy_dfu_deposit.zip\n\n"
        "The recombined file is byte-identical to the original archive.\n"
    )
    REASM.write_text(reasm, encoding="utf-8")
    print(f"[reasm] wrote {REASM.name}")


def main():
    print(f"building Zenodo deposit archive from {HERE}")
    zp = build_zip()
    split_zip(zp)
    print("done.")


if __name__ == "__main__":
    main()
