#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deposit_zenodo.py -- one-command Zenodo deposit for the Kpy x DFU
(lysine propionylation x diabetic foot ulcer) reproducibility package (v1.0.0).

Adapted from the proven CPSP DRG-repurposing deposit script (Zenodo InvenioRDM
REST API, /api/records). Retains the robust upload/resume/publish flow and the
post-publish repository backfill (--fill-doi).

Verified gotchas handled:
  * top-level "pids" block when creating the draft -> HTTP 500 (we omit it).
  * missing metadata.publisher -> publish fails HTTP 400 (we set "Zenodo").
  * related_identifiers.relation_type must be lowercase id ("isderivedfrom").
  * large single PUT stalls through the egress proxy -> split into 1.5 MiB volumes.

Token resolution order: --token  >  env ZENODO_TOKEN  >  ./zenodo_token.txt
Never echo the token.

Usage
-----
  # rehearsal on the sandbox (separate sandbox token)
  python deposit_zenodo.py --sandbox --token <SANDBOX_TOKEN> --publish

  # production: create draft, upload split volumes, publish, then backfill the DOI
  python deposit_zenodo.py --token <PROD_TOKEN> --publish --fill-doi

  # show the exact metadata payload without touching the network
  python deposit_zenodo.py --dry-run
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import re
import subprocess
import sys
import time
import http.client
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent          # 07_复现包 (== the reproducibility repo)
PROJECT = HERE.parent                           # 丙酮酰化与糖尿病足
PROD_BASE = "https://zenodo.org"
SANDBOX_BASE = "https://sandbox.zenodo.org"

DEFAULT_FILES = ["kpy_dfu_deposit.zip"]
TOKEN_FILE = HERE / "zenodo_token.txt"

# Repository files touched by --fill-doi
MANUSCRIPT = PROJECT / "05_稿件" / "manuscript_en.md"
MANUSCRIPT_CN = PROJECT / "05_稿件" / "manuscript_cn.md"
BUILDER = PROJECT / "02_脚本" / "md_to_docx_manuscript.py"
DOCX_OUT = PROJECT / "05_稿件" / "manuscript_en.docx"
DOCX_COPY = PROJECT / "08_投稿包" / "主稿件" / "manuscript_en.docx"
SNAPP_TXT = PROJECT / "08_投稿包" / "数据可用性声明_系统提交用.txt"
CITATION = HERE / "CITATION.cff"
README = HERE / "README.md"
DOI_PLACEHOLDER = "10.5281/zenodo.XXXXXXX"

# ---------------------------------------------------------------------------
# Metadata (Kpy x DFU v1.0.0)
# ---------------------------------------------------------------------------
TITLE = (
    "Lysine propionylation (Kpy) and diabetic foot ulcer: a multi-omics "
    "re-analysis reproducibility package (v1.0.0)"
)
AFFILIATION = (
    "The Second Affiliated Hospital of Fujian University of Traditional Chinese "
    "Medicine, Fuzhou, Fujian 350003, China"
)
ORCID = "0009-0004-9698-6552"
VERSION = "v1.0.0"

DESCRIPTION = (
    "Reproducibility package for the lysine propionylation (Kpy) x diabetic foot "
    "ulcer (DFU) multi-omics re-analysis (single author; re-analysis of public "
    "data only, no wet-lab Kpy measurement executed).<br><br>"
    "<b>What it contains.</b> The machine-readable per-number audit ledger (every "
    "numeric value in the manuscript traced to its source file and extraction "
    "expression), all consolidated result tables (two-cohort Stouffer meta-analysis "
    "of Kpy-residue regulation, gene-set statistics with set-level correction, "
    "dual machine-learning hub-gene prediction with independent external "
    "validation, immune-deconvolution and pan-immune profiling, single-cell and "
    "spatial localisation of hub genes, a two-sample Mendelian-randomisation "
    "sensitivity on the Kpy-enzyme / glycolysis axis, and structure-based docking "
    "of tractable hub targets against the ChEMBL approved-drug library), the "
    "analysis scripts that generated them, the Figure 1-10 source PNGs (>=300 "
    "DPI), the supplementary tables S1-S11, and the SHA-256 MANIFEST.<br><br>"
    "<b>Honest negative findings.</b> The study reports, as its intended result, "
    "that the transcriptomic signal of Kpy-residue regulation in DFU is modest and "
    "that several hypotheses (strong immune-reprogramming enrichment, robust "
    "single-cell/MR positive signals) were not supported; these nulls are "
    "transparently documented rather than minimised. Every number traces to the "
    "files here; nothing was fabricated.<br><br>"
    "<b>Excluded by design.</b> The ~GB-scale raw GEO counts (01_数据/原始数据) and "
    "preprocessed intermediates are not bundled: they are regenerable from the 14 "
    "public GEO accessions (inventoried in Supplementary Table S2) and the "
    "documented PDB IDs / ChEMBL SMILES. The canonical source is the versioned "
    "GitHub release v1.0.0 at https://github.com/yyx-4113/kpy-dfu-transcriptome ."
    "<br><br>"
    "<b>File layout.</b> This deposit is published as sequential split volumes "
    "(kpy_dfu_deposit.zip.001, .002, ...) because the upload egress stalls on "
    "single transfers larger than ~2 MiB. Download all volumes into one folder and "
    "recombine with <code>cat kpy_dfu_deposit.zip.* &gt; kpy_dfu_deposit.zip</code> "
    "(Linux/macOS/WSL) or the equivalent PowerShell/split command (see "
    "00_REASSEMBLE.txt); the recombined file is byte-identical to the original zip. "
    "Unzip it to obtain the audit ledger, derived tables, figures, code and "
    "MANIFEST."
)

KEYWORDS = [
    "lysine propionylation",
    "Kpy",
    "diabetic foot ulcer",
    "diabetes complications",
    "multi-omics",
    "meta-analysis",
    "machine learning",
    "single-cell RNA-seq",
    "spatial transcriptomics",
    "Mendelian randomisation",
    "molecular docking",
    "drug repurposing",
    "reproducibility",
]

RELATED = [{
    "relation_type": {"id": "isderivedfrom"},
    "scheme": "url",
    "identifier": "https://github.com/yyx-4113/kpy-dfu-transcriptome/releases/tag/v1.0.0",
}]


def build_metadata() -> dict:
    today = _dt.date.today().isoformat()
    return {
        "metadata": {
            "title": TITLE,
            "description": DESCRIPTION,
            "publication_date": today,
            "resource_type": {"id": "software"},
            "creators": [
                {
                    "person_or_org": {
                        "type": "personal",
                        "name": "Yang, Yongxin",
                        "given_name": "Yongxin",
                        "family_name": "Yang",
                        "identifiers": [{"scheme": "orcid", "identifier": ORCID}],
                    },
                    "affiliations": [{"name": AFFILIATION}],
                    "role": {"id": "researcher"},
                }
            ],
            "rights": [{"id": "mit"}],
            "languages": [{"id": "eng"}],
            "keywords": KEYWORDS,
            "related_identifiers": RELATED,
            "publisher": "Zenodo",
            "version": VERSION,
        }
    }


# ---------------------------------------------------------------------------
# Thin HTTP helper (unchanged robust implementation)
# ---------------------------------------------------------------------------
def _req(base, path, token, method="GET", json_body=None, raw_body=None,
         content_type=None, timeout=300, extra_headers=None):
    url = path if str(path).startswith("http") else base + path
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json"}
    if extra_headers:
        headers.update(extra_headers)
    data = None
    if raw_body is not None:
        data = raw_body
        headers["Content-Type"] = content_type or "application/octet-stream"
    elif json_body is not None:
        data = json.dumps(json_body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = resp.read()
            if not payload:
                return resp.status, {}
            try:
                return resp.status, json.loads(payload)
            except Exception:
                return resp.status, {"_raw": payload[:800].decode("utf-8", "replace")}
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, {"error": body[:1500]}
    except (urllib.error.URLError, OSError, http.client.HTTPException) as e:
        return 0, {"error": f"network unreachable: {e!r}"}


def _fail(msg, payload=None):
    print(f"  \u2717 {msg}", file=sys.stderr)
    if payload is not None:
        print(f"    response: {json.dumps(payload)[:1500]}", file=sys.stderr)
    sys.exit(1)


def extract_version_doi(rec):
    if not isinstance(rec, dict):
        return None
    pids = rec.get("pids") or {}
    d = pids.get("doi") or {}
    if d.get("identifier"):
        return d["identifier"]
    if rec.get("doi"):
        return rec["doi"]
    l = (rec.get("links") or {}).get("doi")
    if l:
        return l.rstrip("/").split("doi.org/")[-1]
    return None


def extract_concept_doi(rec):
    if not isinstance(rec, dict):
        return None
    links = rec.get("links") or {}
    if links.get("parent_doi"):
        return links["parent_doi"].rstrip("/").split("doi.org/")[-1]
    parent = rec.get("parent") or {}
    pd = (parent.get("pids") or {}).get("doi") or {}
    if pd.get("identifier"):
        return pd["identifier"]
    if rec.get("conceptdoi"):
        return rec["conceptdoi"]
    if links.get("conceptdoi"):
        return links["conceptdoi"].rstrip("/").split("doi.org/")[-1]
    return None


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------
def create_draft(base, token):
    print("[1/6] creating empty draft record ...")
    body = {
        "access": {"record": "public", "files": "public"},
        "files": {"enabled": True},
        **build_metadata(),
    }
    status, resp = _req(base, "/api/records", token, method="POST", json_body=body)
    if status not in (200, 201):
        _fail(f"draft creation failed (HTTP {status})", resp)
    rid = resp.get("id")
    print(f"  \u2713 draft id={rid}  edit at {(resp.get('links') or {}).get('self_html', '')}")
    return resp


def draft_committed_keys(base, token, rid):
    for attempt in range(1, 5):
        status, resp = _req(base, f"/api/records/{rid}/draft/files", token, method="GET")
        if status == 200:
            ents = resp.get("entries", []) if isinstance(resp, dict) else (resp if isinstance(resp, list) else [])
            done = set()
            for f in ents:
                if f.get("status") == "completed":
                    done.add(f.get("key"))
            return done
        print(f"  ! committed-keys GET attempt {attempt} HTTP {status}; "
              f"retrying in {min(2 ** attempt, 30)}s")
        time.sleep(min(2 ** attempt, 30))
    return set()


def upload_files(base, token, draft, paths, skip_existing=True):
    rid = draft["id"]
    done = draft_committed_keys(base, token, rid) if skip_existing else set()
    failed = []
    for i, p in enumerate(paths, 1):
        p = Path(p)
        if not p.exists():
            _fail(f"file not found: {p}")
        if skip_existing and p.name in done:
            print(f"[3/6] file {i}/{len(paths)}: {p.name} already committed -- skip")
            continue
        size = p.stat().st_size
        print(f"[3/6] registering file {i}/{len(paths)}: {p.name} ({size:,} bytes) ...")
        registered = False
        for attempt in range(1, 6):
            status, resp = _req(base, f"/api/records/{rid}/draft/files", token,
                                method="POST", json_body=[{"key": p.name}])
            if status in (200, 201, 409):
                registered = True
                break
            print(f"    ! registration attempt {attempt} HTTP {status}: "
                  f"{str(resp)[:80]}; retrying in {min(2 ** attempt, 30)}s")
            time.sleep(min(2 ** attempt, 30))
        if not registered:
            print(f"  \u2717 {p.name} registration FAILED this run; will retry on next run")
            failed.append(p.name)
            continue
        content_url = f"/api/records/{rid}/draft/files/{p.name}/content"
        print(f"        uploading bytes (single PUT, retry on network reset) ...")
        data = p.read_bytes()
        ok = False
        for attempt in range(1, 6):
            status, resp = _req(base, content_url, token, method="PUT",
                                raw_body=data, content_type="application/octet-stream",
                                timeout=900)
            if status in (200, 201):
                ok = True
                break
            print(f"    ! single PUT attempt {attempt} HTTP {status}: "
                  f"{str(resp)[:100]}; retrying in {min(2 ** attempt, 30)}s")
            time.sleep(min(2 ** attempt, 30))
        if not ok:
            print(f"  \u2717 {p.name} FAILED this run (HTTP {status}); will retry on next run")
            failed.append(p.name)
            continue
        print(f"        committing ...")
        status, resp = _req(base, f"/api/records/{rid}/draft/files/{p.name}/commit",
                            token, method="POST")
        if status not in (200, 201, 202):
            print(f"  \u2717 {p.name} commit FAILED this run (HTTP {status}); will retry on next run")
            failed.append(p.name)
            continue
        print(f"  \u2713 {p.name} uploaded and committed")
        time.sleep(1)
    return failed


def update_draft_metadata(base, token, draft):
    rid = draft["id"]
    print("[5/6] updating draft metadata (re-PUT corrected schema) ...")
    status, resp = _req(base, f"/api/records/{rid}/draft", token, method="PUT",
                        json_body=build_metadata())
    if status not in (200, 201):
        _fail(f"draft metadata update failed (HTTP {status})", resp)
    print(f"  \u2713 draft metadata updated")
    return resp


def publish(base, token, draft):
    rid = draft["id"]
    print("[6/6] publishing ...")
    status, resp = _req(base, f"/api/records/{rid}/draft/actions/publish",
                        token, method="POST")
    if status not in (200, 202):
        _fail(f"publish failed (HTTP {status})", resp)
    st2, rec = _req(base, f"/api/records/{rid}", token, method="GET")
    if st2 not in (200,):
        rec = resp
    vdoi = extract_version_doi(rec)
    cdoi = extract_concept_doi(rec)
    html = (rec.get("links") or {}).get("self_html") or (rec.get("links") or {}).get("record_html")
    print(f"  \u2713 published: {html}")
    if cdoi:
        print(f"  \u2713 concept DOI : {cdoi}")
    if vdoi:
        print(f"  \u2713 version DOI : {vdoi}")
    return rec, vdoi, cdoi


def new_version(base, token, old_id, paths):
    print(f"[new-version] creating new version of record {old_id} ...")
    st, r = _req(base, f"/api/records/{old_id}/versions", token, method="POST")
    if st not in (200, 201):
        _fail(f"new-version creation failed (HTTP {st})", r)
    new_id = r.get("id")
    print(f"  \u2713 new draft id={new_id}")
    st, r = _req(base, f"/api/records/{new_id}/draft", token, method="PUT",
                 json_body=build_metadata())
    if st not in (200, 201):
        _fail(f"new-version metadata update failed (HTTP {st})", r)
    upload_files(base, token, {"id": new_id}, paths)
    rec, vdoi, cdoi = publish(base, token, {"id": new_id})
    return rec, vdoi, cdoi


# ---------------------------------------------------------------------------
# Post-publish repository backfill (--fill-doi)
# ---------------------------------------------------------------------------
EN_OLD = (
    "will be released in the public GitHub repository "
    "https://github.com/yyx-4113/kpy-dfu-transcriptome (tag v1.0.0); the "
    "repository MANIFEST records each file's path, size and SHA-256 checksum. "
    "The repository will be archived on Zenodo with a versioned DOI upon acceptance."
)
EN_NEW = (
    "are released in the public GitHub repository "
    "https://github.com/yyx-4113/kpy-dfu-transcriptome (tag v1.0.0); the "
    "repository MANIFEST records each file's path, size and SHA-256 checksum. "
    "The repository is archived on Zenodo under DOI {doi} (https://doi.org/{doi})."
)

CN_OLD = (
    "将于公共 GitHub 仓库 https://github.com/yyx-4113/kpy-dfu-transcriptome"
    "（tag v1.0.0）发布；仓库 MANIFEST 记录每个文件的路径、大小与 SHA-256 "
    "校验和。该仓库将在接受后归档于 Zenodo 并取得版本化 DOI。"
)
CN_NEW = (
    "已于公共 GitHub 仓库 https://github.com/yyx-4113/kpy-dfu-transcriptome"
    "（tag v1.0.0）发布；仓库 MANIFEST 记录每个文件的路径、大小与 SHA-256 "
    "校验和。该仓库已归档于 Zenodo，版本化 DOI：{doi}（https://doi.org/{doi}）。"
)


def _bake_docx(doi):
    if BUILDER.exists():
        print(f"  rebuilding Manuscript.docx via {BUILDER.name} ...")
        r = subprocess.run([sys.executable, str(BUILDER)], cwd=str(PROJECT))
        if r.returncode == 0 and DOCX_OUT.exists():
            DOCX_COPY.parent.mkdir(parents=True, exist_ok=True)
            DOCX_COPY.write_bytes(DOCX_OUT.read_bytes())
            print(f"  \u2713 {DOCX_OUT.name} rebuilt; copied to 08_投稿包")
        else:
            print("  ! docx rebuild reported failure; left as-is")
    else:
        print(f"  ! builder not found at {BUILDER}; skipped docx rebuild")


def fill_doi_in_repo(doi):
    print(f"[fill-doi] backfilling DOI {doi} into the repository ...")
    changed = []

    if MANUSCRIPT.exists():
        txt = MANUSCRIPT.read_text(encoding="utf-8")
        if EN_OLD in txt:
            MANUSCRIPT.write_text(txt.replace(EN_OLD, EN_NEW.format(doi=doi)), encoding="utf-8")
            changed.append("manuscript_en.md Data Availability -> DOI")
        else:
            print("  ! EN DA placeholder not found verbatim; left unchanged "
                  "(inspect manually)")

    if MANUSCRIPT_CN.exists():
        txt = MANUSCRIPT_CN.read_text(encoding="utf-8")
        if CN_OLD in txt:
            MANUSCRIPT_CN.write_text(txt.replace(CN_OLD, CN_NEW.format(doi=doi)), encoding="utf-8")
            changed.append("manuscript_cn.md 数据与材料可用性 -> DOI")
        else:
            print("  ! CN DA placeholder not found verbatim; left unchanged")

    _bake_docx(doi)

    if CITATION.exists():
        s = CITATION.read_text(encoding="utf-8")
        if "identifiers:" not in s:
            block = (
                "identifiers:\n"
                "  - type: doi\n"
                f"    value: {doi}\n"
                "    description: Concept DOI resolving to the versioned "
                "reproducibility deposit on Zenodo.\n"
            )
            s = s.replace("version: v1.0.0\n", "version: v1.0.0\n" + block, 1)
            CITATION.write_text(s, encoding="utf-8")
            changed.append("CITATION.cff identifiers: doi added")
        else:
            print("  ! CITATION.cff already carries the DOI; left unchanged")

    if README.exists():
        s = README.read_text(encoding="utf-8")
        if doi not in s:
            badge = (f"[![DOI](https://zenodo.org/badge/DOI/{doi}.svg)]"
                     f"(https://doi.org/{doi})\n")
            s = re.sub(r"^(# .*)$", r"\1\n\n" + badge, s, count=1, flags=re.M)
            README.write_text(s, encoding="utf-8")
            changed.append("README.md Zenodo badge added")
        else:
            print("  ! README.md already carries the DOI; left unchanged")

    if SNAPP_TXT.exists():
        txt = SNAPP_TXT.read_text(encoding="utf-8")
        new_txt = txt.replace(
            "will be archived on Zenodo with a versioned DOI upon acceptance.",
            f"is archived on Zenodo under DOI {doi} (https://doi.org/{doi}).")
        if new_txt != txt:
            SNAPP_TXT.write_text(new_txt, encoding="utf-8")
            changed.append("08_投稿包/数据可用性声明_系统提交用.txt -> DOI")
        else:
            print("  ! SNAPP txt 'upon acceptance' phrase not found; left unchanged")

    for c in changed:
        print(f"  \u2713 {c}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def resolve_token(args):
    if args.token:
        return args.token.strip()
    if os.environ.get("ZENODO_TOKEN"):
        return os.environ["ZENODO_TOKEN"].strip()
    if TOKEN_FILE.exists():
        return TOKEN_FILE.read_text(encoding="utf-8").strip()
    return None


def main():
    ap = argparse.ArgumentParser(
        description="One-command Zenodo deposit for the Kpy x DFU reproducibility package.")
    ap.add_argument("--token", help="Zenodo personal access token (deposit:write + deposit:actions). "
                                   "Falls back to $ZENODO_TOKEN or ./zenodo_token.txt")
    ap.add_argument("--sandbox", action="store_true",
                    help="use https://sandbox.zenodo.org (needs a separate sandbox token)")
    ap.add_argument("--files", nargs="+", default=DEFAULT_FILES,
                    help="deposit file(s) (default: kpy_dfu_deposit.zip)")
    ap.add_argument("--publish", action="store_true",
                    help="publish after upload (without it, the record stays a draft)")
    ap.add_argument("--new-version", metavar="RECORD_ID", type=int,
                    help="create a NEW VERSION of an existing record, upload --files, and publish")
    ap.add_argument("--resume", metavar="DRAFT_ID", type=int,
                    help="resume an EXISTING draft (skip create_draft)")
    ap.add_argument("--fill-doi", action="store_true",
                    help="after publishing, backfill the DOI into the manuscripts, "
                         "CITATION.cff, README, SNAPP txt and rebuild the .docx")
    ap.add_argument("--patch-doi", metavar="DOI",
                    help="ONLY backfill the DOI into the repo and rebuild the docx "
                         "(no network, no token needed)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the exact payload and plan without touching the network")
    ap.add_argument("--no-proxy", action="store_true",
                    help="disable the egress HTTP proxy and connect directly to Zenodo")
    args = ap.parse_args()

    if args.no_proxy:
        for _k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY",
                   "all_proxy", "ALL_PROXY"):
            os.environ.pop(_k, None)
        print("[no-proxy] egress HTTP proxy disabled; connecting directly to Zenodo")

    if args.patch_doi:
        fill_doi_in_repo(args.patch_doi.strip())
        return

    base = SANDBOX_BASE if args.sandbox else PROD_BASE
    if args.files == DEFAULT_FILES:
        single = HERE / "kpy_dfu_deposit.zip"
        vols = sorted(HERE.glob("kpy_dfu_deposit.zip.[0-9][0-9][0-9]"))
        reasm = HERE / "00_REASSEMBLE.txt"
        if vols:
            files = [str(p) for p in vols] + ([str(reasm)] if reasm.exists() else [])
            print(f"[info] using {len(vols)} split volumes + reassembly helper")
        elif single.exists():
            files = [str(single)]
            print("[info] using single kpy_dfu_deposit.zip")
        else:
            _fail("no deposit archive found; run make_deposit.py first")
    else:
        files = [str((HERE / f) if not Path(f).is_absolute() else Path(f)) for f in args.files]

    print("=" * 72)
    print("Zenodo deposit -- Kpy x DFU reproducibility package")
    print(f"  host    : {base}{'   (SANDBOX)' if args.sandbox else '   (PRODUCTION)'}")
    print(f"  files   : {', '.join(Path(f).name for f in files)}")
    print(f"  publish : {args.publish}   fill-doi: {args.fill_doi}   "
          f"new-version: {args.new_version or '-'}")
    print("=" * 72)

    if args.dry_run:
        print("[dry-run] metadata payload that would be POSTed to /api/records:")
        print(json.dumps(build_metadata(), indent=2, ensure_ascii=False))
        print("\n[dry-run] planned calls:")
        print(f"  POST {base}/api/records")
        for f in files:
            print(f"  POST {base}/api/records/<id>/draft/files   key={Path(f).name}")
            print(f"  PUT  {base}/api/records/<id>/draft/files/{Path(f).name}/content   "
                  f"({Path(f).stat().st_size if Path(f).exists() else '?'} bytes)")
            print(f"  POST {base}/api/records/<id>/draft/files/{Path(f).name}/commit")
        if args.publish:
            print(f"  POST {base}/api/records/<id>/draft/actions/publish")
        return

    token = resolve_token(args)
    if not token:
        _fail("no token. Pass --token, set $ZENODO_TOKEN, or create ./zenodo_token.txt (never commit it)")

    if args.new_version:
        rec, vdoi, cdoi = new_version(base, token, args.new_version, files)
        doi = cdoi or vdoi
        if args.fill_doi and doi:
            fill_doi_in_repo(doi)
    else:
        draft_id = args.resume
        if draft_id:
            draft = {"id": draft_id}
            print(f"[resume] reusing existing draft id={draft_id}")
        else:
            draft = create_draft(base, token)
            draft_id = draft["id"]
        print(f"[info] draft id={draft_id} -- upload is resumable: volumes already "
              f"committed in the draft are skipped automatically")
        failed = upload_files(base, token, draft, files)
        if failed:
            print()
            print("=" * 72)
            print(f"UPLOAD INCOMPLETE: {len(failed)}/{len(files)} volumes not yet "
                  f"committed on draft {draft_id}:")
            for k in failed:
                print(f"  - {k}")
            print("Re-run the SAME command with --resume "
                  f"{draft_id} to continue; committed volumes are skipped.")
            print("=" * 72)
            doi = None
        else:
            if args.publish:
                update_draft_metadata(base, token, draft)
                rec, vdoi, cdoi = publish(base, token, draft)
                doi = cdoi or vdoi
                if args.fill_doi and doi:
                    fill_doi_in_repo(doi)
            else:
                print(f"[done] draft left unpublished. Review and publish at: "
                      f"{(draft.get('links') or {}).get('self_html')}")
                doi = None

    print()
    print("=" * 72)
    if doi:
        print(f"RESULT_DOI={doi}")
        print(f"RESULT_URL=https://doi.org/{doi}")
    else:
        print("RESULT_DOI=(none yet - record not published)")
    print("=" * 72)


if __name__ == "__main__":
    main()
