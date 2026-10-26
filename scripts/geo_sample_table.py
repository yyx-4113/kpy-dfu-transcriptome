# -*- coding: utf-8 -*-
"""
geo_sample_table.py — 拉取候选队列的样本级分组信息 + suppl 数据文件清单
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
输出：01_数据/元数据/样本表/<GSE>_samples.csv
      01_数据/元数据/样本表/<GSE>_suppl_files.csv
"""
import csv

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import gzip
import io
import os
import re
import tarfile
import urllib.request

BASE_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
OUT_DIR = _ROOT + "/01_数据/元数据/样本表"

GSES = [
    "GSE134431", "GSE80178", "GSE68183", "GSE68186", "GSE199939",
    "GSE143735", "GSE166120", "GSE165816", "GSE230426", "GSE223964",
    "GSE245703", "GSE231643", "GSE146028", "GSE114248",
]


def series_dir(gse: str) -> str:
    return f"{gse[:-3]}nnn"


def get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "kpy-dfu-project/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def list_dir(url: str) -> list:
    html = get(url, 60).decode("utf-8", "replace")
    return re.findall(r'href="([^"][^"]*)"', html)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for gse in GSES:
        sd = series_dir(gse)
        base = f"{BASE_FTP}/{sd}/{gse}"
        # --- suppl 文件清单 ---
        suppl_rows = []
        try:
            for href in list_dir(f"{base}/suppl/"):
                if href.startswith("/") or href.startswith("http") or "hhs.gov" in href:
                    continue
                suppl_rows.append(href)
            with open(f"{OUT_DIR}/{gse}_suppl_files.csv", "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["file"])
                for r in suppl_rows:
                    w.writerow([r])
        except Exception as e:
            suppl_rows = [f"ERR {e}"]

        # --- MINiML 样本表 ---
        rows = []
        try:
            raw = get(f"{base}/miniml/{gse}_family.xml.tgz", 180)
            tf = tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz")
            for m in tf.getmembers():
                if not m.name.endswith(".xml"):
                    continue
                import xml.etree.ElementTree as ET
                data = tf.extractfile(m).read()
                NS = "{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}"
                try:
                    root = ET.fromstring(data)
                    ns = NS
                except ET.ParseError:
                    root = ET.fromstring(re.sub(rb'\s+xmlns="[^"]+"', b"", data, count=1))
                    ns = ""
                for s in root.iter(ns + "Sample"):
                    acc = s.get("accession") or (s.findtext(ns + "Accession") or "").strip()
                    title = (s.findtext(ns + "Title") or "").strip()
                    chars = []
                    for ch in s.findall(".//" + ns + "Characteristics"):
                        tag = ch.get("tag", "")
                        txt = (ch.text or "").strip()
                        chars.append(f"{tag}={txt}" if tag else txt)
                    ptype = (s.findtext(ns + "Type") or "").strip()
                    plat = ""
                    pr = s.find(ns + "Platform-Ref")
                    if pr is not None:
                        plat = pr.get("ref") or ""
                    rows.append({
                        "gsm": acc, "title": title, "type": ptype,
                        "platform": plat, "characteristics": " | ".join(chars)[:400],
                    })
        except Exception as e:
            print(f"[FAIL] {gse}: {e}")

        if rows:
            out = f"{OUT_DIR}/{gse}_samples.csv"
            with open(out, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.DictWriter(f, fieldnames=["gsm", "title", "type", "platform", "characteristics"])
                w.writeheader()
                w.writerows(rows)
            print(f"[OK] {gse}: {len(rows)} samples; suppl files: {len(suppl_rows)}")
            for r in rows[:40]:
                print(f"    {r['gsm']:<11} {r['type'][:14]:<14} {r['characteristics'][:160]}")
        print(f"    suppl: {suppl_rows[:8]}")
        print("-" * 100)


if __name__ == "__main__":
    main()
