#!/usr/bin/env python3
"""Snapshot processed OSDR inputs without modifying the legacy workspace.

This is the only script that needs the original local downloads. All analysis
scripts operate on the portable inputs/ directory it creates. Local metadata
contains derived cohort fields; these are preserved and explicitly identified.
"""
import argparse
import gzip
import hashlib
import io
import json
import re
import shutil
import zipfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATASETS = {
    "525": ("RR-23", "cerebellum"),
    "561": ("RRRM-2", "cerebellum"),
    "562": ("RRRM-2", "hippocampus"),
    "613": ("RRRM-2", "left cerebral hemisphere"),
    "563": ("RR-10", "cerebellum"),
    "564": ("RR-10", "hippocampus"),
    "612": ("RR-10", "left cerebral hemisphere"),
}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=HERE.parents[1])
    parser.add_argument("--local-data", type=Path, required=True)
    parser.add_argument("--metadata-cache", type=Path)
    args = parser.parse_args()
    dest = HERE / "inputs"
    dest.mkdir(exist_ok=True)
    records = []

    def snapshot(source, target):
        # Deterministic gzip; decompressed bytes are identical to the source.
        target = dest / target
        with source.open("rb") as inp, target.open("wb") as out:
            with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0) as gz:
                shutil.copyfileobj(inp, gz)
        records.append({"source_path": str(source.resolve()),
                        "source_sha256": sha(source), "bundled_path": target.name,
                        "bundled_sha256": sha(target), "source_bytes": source.stat().st_size})

    for n in DATASETS:
        if n in ("561", "562"):
            snapshot(args.local_data / f"counts/OSD-{n}_counts.csv", f"OSD-{n}_counts.csv.gz")
            snapshot(args.local_data / f"meta/OSD-{n}_metadata.csv", f"OSD-{n}_metadata.csv.gz")
        else:
            base = args.repo / f"research/brain_integrative_2026/data/OSD-{n}"
            files = sorted(base.glob("*STAR_Unnormalized_Counts*.csv"))
            if len(files) != 1:
                raise ValueError(f"Expected exactly one STAR matrix in {base}")
            snapshot(files[0], f"OSD-{n}_counts.csv.gz")
            if n == "525":
                rsem = next(base.glob("*RSEM_Unnormalized_Counts*.csv"))
                snapshot(rsem, "OSD-525_RSEM_counts.csv.gz")
                isa = base / "OSD-525_metadata_OSD-525-ISA.zip"
                with zipfile.ZipFile(isa) as z:
                    md = pd.read_csv(io.BytesIO(z.read("s_OSD-525.txt")), sep="\t")
                md = md.rename(columns={"Sample Name": "sample"})
                md.to_csv(dest / "OSD-525_metadata.csv.gz", index=False,
                          compression={"method": "gzip", "mtime": 0})
                records.append({"source_path": str(isa.resolve()), "source_sha256": sha(isa),
                                "bundled_path": "OSD-525_metadata.csv.gz",
                                "bundled_sha256": sha(dest / "OSD-525_metadata.csv.gz"),
                                "transformation": "s_OSD-525.txt TSV to CSV; Sample Name renamed sample"})
            else:
                snapshot(args.repo / f"research/brain_integrative_2026/tables/refit/OSD-{n}_metadata.csv",
                         f"OSD-{n}_metadata.csv.gz")

    gtf = args.repo / "research/osdr_brain_paper_2026/data/Mus_musculus.GRCm39.112.gtf.gz"
    genes = []
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', parts[8]))
            genes.append({"ensembl_id": a["gene_id"], "gene_symbol": a.get("gene_name", ""),
                          "gene_biotype": a.get("gene_biotype", "")})
    pd.DataFrame(genes).to_csv(dest / "gene_annotation.csv.gz", index=False,
                              compression={"method": "gzip", "mtime": 0})
    records.append({"source_path": str(gtf.resolve()), "source_sha256": sha(gtf),
                    "bundled_path": "gene_annotation.csv.gz",
                    "bundled_sha256": sha(dest / "gene_annotation.csv.gz"),
                    "transformation": "GRCm39 Ensembl release 112 gene features; ID, symbol and biotype"})
    if args.metadata_cache:
        for p in sorted(args.metadata_cache.glob("osd_*_public_metadata.json")):
            shutil.copyfile(p, dest / p.name)
            records.append({"source_path": str(p), "bundled_path": p.name,
                            "source_sha256": sha(p), "bundled_sha256": sha(dest / p.name),
                            "retrieval_date": "2026-10-14", "source": "NASA OSDR metadata API"})
    (dest / "provenance.json").write_text(json.dumps(records, indent=2) + "\n")
    print(f"Snapshotted {len(records)} inputs into {dest}")


if __name__ == "__main__":
    main()
