#!/usr/bin/env python3
"""Download the inputs for the OSD-682/685/698/699 GeoMx DSP re-analysis.

Two sources are needed, because they carry different things:

  NASA OSDR  -- ISA metadata archives (factor values, mission parameters,
                sequencing parameters).  OSDR hosts *only* raw FASTQ for these
                four studies (7-11 GB per study) plus a raw-read MultiQC
                report; there are no GeneLab-processed count or differential-
                expression tables, so the processed layer has to come from GEO.

  NCBI GEO    -- GSE239336, the submitter-deposited processed layer:
                per-region Q3-normalised expression matrices, per-region
                differential-expression exports from the NanoString DSP
                Analysis Suite, the series matrix (per-ROI design + QC
                metadata), and the raw per-ROI DCC files.

Total download is roughly 12 MB.  The raw FASTQ files are deliberately not
fetched; re-processing them would require the GeoMx DND pipeline and the
instrument configuration files, which are not deposited.

Usage:  python scripts/download_data.py [--with-dcc]
"""
from __future__ import annotations

import argparse
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import common as C

OSDR_URL = (
    "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download"
    "?source=datamanager&file={fname}"
)
GEO_SUPPL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/{series}/suppl/{fname}"
)
GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/{series}/matrix/{fname}"
)

DE_CONTRASTS = list(C.PROVIDED_CONTRASTS)


def download(url: str, dest: Path, retries: int = 3) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  [skip] {dest.name} ({dest.stat().st_size/1e6:.2f} MB, present)")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "osdr-geomx-reanalysis/1.0"})
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = resp.read()
            if not data:
                raise urllib.error.URLError("empty response")
            dest.write_bytes(data)
            print(f"  [ok]   {dest.name} ({len(data)/1e6:.2f} MB)")
            return True
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            print(f"  [try {attempt}/{retries}] {dest.name}: {exc}")
            time.sleep(2 * attempt)
    print(f"  [FAIL] {dest.name}")
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--with-dcc",
        action="store_true",
        help="also unpack the 48 raw per-ROI DCC files from GSE239336_RAW.tar",
    )
    args = ap.parse_args()

    ok = True

    print("NCBI GEO -- GSE239336 processed layer")
    for study in C.STUDIES:
        name = f"{C.GEO_SERIES}_{study.code}_GeneExpression_Q3norm.txt.gz"
        ok &= download(GEO_SUPPL.format(series=C.GEO_SERIES, fname=name), C.GEO_DIR / name)
        for contrast in DE_CONTRASTS:
            name = f"{C.GEO_SERIES}_{study.code}_{contrast}_DEanalysis.txt.gz"
            ok &= download(GEO_SUPPL.format(series=C.GEO_SERIES, fname=name), C.GEO_DIR / name)

    name = f"{C.GEO_SERIES}_series_matrix.txt.gz"
    ok &= download(GEO_MATRIX.format(series=C.GEO_SERIES, fname=name), C.GEO_DIR / name)

    name = f"{C.GEO_SERIES}_RAW.tar"
    raw_tar = C.GEO_DIR / name
    ok &= download(GEO_SUPPL.format(series=C.GEO_SERIES, fname=name), raw_tar)

    if args.with_dcc and raw_tar.exists():
        dcc_dir = C.GEO_DIR / "dcc"
        dcc_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(raw_tar) as tf:
            members = [m for m in tf.getmembers() if m.isfile()]
            for m in members:
                m.name = Path(m.name).name  # flatten, guard against path traversal
            tf.extractall(dcc_dir, members=members)
        print(f"  [ok]   unpacked {len(members)} DCC files -> {dcc_dir.relative_to(C.ROOT)}")

    print("\nNASA OSDR -- ISA metadata archives")
    for study in C.STUDIES:
        name = f"{study.osd}_metadata_{study.osd}-ISA.zip"
        dest = C.OSDR_DIR / name
        if download(OSDR_URL.format(osd=study.osd, fname=name), dest):
            out = C.OSDR_DIR / study.osd
            out.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(dest) as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        continue
                    target = out / Path(info.filename).name
                    target.write_bytes(zf.read(info))
            print(f"         unpacked -> {out.relative_to(C.ROOT)}")
        else:
            ok = False

    print("\nDone." if ok else "\nCompleted with some failures.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
