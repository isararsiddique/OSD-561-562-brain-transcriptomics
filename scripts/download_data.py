#!/usr/bin/env python3
"""Download GeneLab-processed bulk RNA-seq tables for the CNS Injury EV project.

Only small processed CSV tables are downloaded (a few MB each) -- NOT the raw
FASTQ / BAM files (which are 100s of GB combined).

Datasets:
  OSD-561  ->  GLDS-556   Cerebellum (Cb), RRRM2 spaceflight, OLD/YNG
  OSD-562  ->  GLDS-557   Hippocampus (HPC), RRRM2 spaceflight, OLD/YNG

Files per study (GeneLab bulk RNA-seq processing pipeline outputs):
  - Normalized_Counts      : DESeq2 median-of-ratios normalized counts (heatmaps)
  - VST_Counts             : variance-stabilized counts (PCA, clustering, heatmaps)
  - differential_expression: log2FC / p / adj-p for every contrast (volcano plots)
  - SampleTable            : sample -> group mapping
  - contrasts              : contrast definitions

Usage:
  python scripts/download_data.py
"""
from __future__ import annotations

import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

BASE = "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={fname}"

# study id -> (OSD accession, GLDS prefix)
STUDIES = {
    "OSD-561": "GLDS-556",  # Cerebellum
    "OSD-562": "GLDS-557",  # Hippocampus
}

# suffixes of the processed tables we need (prefixed per-study with the GLDS id)
TABLE_SUFFIXES = [
    "rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
    "rna_seq_VST_Counts_GLbulkRNAseq.csv",
    "rna_seq_differential_expression_GLbulkRNAseq.csv",
    "rna_seq_SampleTable_GLbulkRNAseq.csv",
    "rna_seq_contrasts_GLbulkRNAseq.csv",
    "rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
]

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def download(url: str, dest: Path, retries: int = 3) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  [skip] {dest.name} ({dest.stat().st_size/1e6:.1f} MB, already present)")
        return True
    req = urllib.request.Request(url, headers={"User-Agent": "osdr-cns-ev/1.0"})
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = resp.read()
            dest.write_bytes(data)
            print(f"  [ok]   {dest.name} ({len(data)/1e6:.1f} MB)")
            return True
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            print(f"  [try {attempt}/{retries}] {dest.name}: {exc}")
            time.sleep(2 * attempt)
    print(f"  [FAIL] {dest.name}")
    return False


def main() -> int:
    ok = True
    for osd, glds in STUDIES.items():
        out = DATA_DIR / osd
        out.mkdir(parents=True, exist_ok=True)
        print(f"\n{osd}  ({glds})")
        for suffix in TABLE_SUFFIXES:
            fname = f"{glds}_{suffix}"
            url = BASE.format(osd=osd, fname=fname)
            ok &= download(url, out / fname)
    print("\nDone." if ok else "\nCompleted with some failures.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
