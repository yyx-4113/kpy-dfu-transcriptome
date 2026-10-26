# -*- coding: utf-8 -*-
"""
download_data.py — 下载 DFU 核心队列原始数据（幂等，带 MANIFEST 台账）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import datetime as dt
import hashlib
import os
import time
import urllib.request

BASE = _ROOT + "/01_数据"
MANIFEST = f"{BASE}/MANIFEST.csv"

JOBS = {
    "GSE134431": ["GSE134431_181016_diabetic_rna-seq_ANOVA_by_group.xlsx",
                  "GSE134431_181016_diabetic_rna-seq_results.gene.rpkm.xlsx"],
    "GSE199939": ["GSE199939_gene.tpm.matrix.annot_-_21.txt.gz"],
    "GSE80178": ["GSE80178_RAW.tar", "GSE80178_further_processed_data.txt.gz"],
    "GSE143735": ["GSE143735_RAW.tar"],
    "GSE68183": ["GSE68183_RAW.tar"],
    "GSE68186": ["GSE68186_RAW.tar"],
    "GSE165816": ["GSE165816_RAW.tar"],
    "GSE166120": ["GSE166120_RAW.tar"],
    "GSE223964": ["GSE223964_allcell.h5ad.gz"],
    "GSE231643": ["GSE231643_RAW.tar"],
    "GSE245703": ["GSE245703_RAW.tar"],
    "GSE230426": ["GSE230426_RAW.tar"],
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(gse: str, fname: str) -> tuple:
    sd = f"{gse[:-3]}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{sd}/{gse}/suppl/{fname}"
    outdir = f"{BASE}/原始数据/{gse}"
    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{fname}"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "kpy-dfu-project/1.0"})
            with urllib.request.urlopen(req, timeout=1800) as r, open(out, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return out, url, os.path.getsize(out), sha256(out)
        except Exception as e:
            print(f"  [retry {attempt+1}/3] {gse}/{fname}: {e}")
            time.sleep(5)
    raise RuntimeError(f"FAILED {gse}/{fname}")


def main():
    rows = []
    if os.path.exists(MANIFEST):
        with open(MANIFEST, encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
    done = {(r["accession"], r["file"]) for r in rows}
    for gse, files in JOBS.items():
        for fname in files:
            if (gse, fname) in done:
                print(f"[skip] {gse}/{fname} already in MANIFEST")
                continue
            t0 = time.time()
            out, url, size, digest = download(gse, fname)
            rows.append({
                "accession": gse, "file": fname, "url": url, "bytes": size,
                "sha256": digest, "downloaded_at": dt.datetime.now().isoformat(timespec="seconds"),
            })
            with open(MANIFEST, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print(f"[OK] {gse}/{fname}  {size/1e6:.1f} MB  {time.time()-t0:.0f}s  sha256={digest[:12]}")
    print("ALL DONE")


if __name__ == "__main__":
    main()
