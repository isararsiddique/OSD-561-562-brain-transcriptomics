#!/usr/bin/env python3
"""Compare core bundled counts with the named NASA STAR count files.

Download only public processed matrices. Snapshot the API file listing and
downloaded bytes, verify gene/sample IDs and every count, and fail on mismatch.
No files or scientific results at the source repository are modified.
"""
import gzip
import datetime
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]


def main():
    records = []
    for n in ("561", "562"):
        listing_url = f"https://osdr.nasa.gov/osdr/data/osd/files/{n}"
        with urllib.request.urlopen(listing_url, timeout=45) as response:
            listing_bytes = response.read()
        listing = json.loads(listing_bytes)
        matches = [f for f in listing["studies"][f"OSD-{n}"]["study_files"]
                   if "STAR_Unnormalized_Counts_GLbulk" in f["file_name"] and not f.get("restricted")]
        if len(matches) != 1:
            raise ValueError(f"Expected one public primary STAR matrix in OSD-{n}")
        metadata = matches[0]
        url = "https://osdr.nasa.gov" + metadata["remote_url"]
        with urllib.request.urlopen(url, timeout=45) as response:
            raw = response.read()
        if len(raw) != metadata["file_size"]:
            raise ValueError(f"Downloaded size differs from API listing: OSD-{n}")
        target = HERE / f"inputs/upstream_OSD-{n}_STAR.csv.gz"
        with target.open("wb") as stream:
            with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
                compressed.write(raw)
        (HERE / f"inputs/upstream_OSD-{n}_file_listing.json").write_bytes(listing_bytes)
        upstream = pd.read_csv(target, index_col=0)
        bundled = pd.read_csv(HERE / f"inputs/OSD-{n}_counts.csv.gz", index_col=0)
        same = set(upstream.index) == set(bundled.index) and set(upstream.columns) == set(bundled.columns)
        difference = float(np.max(np.abs(upstream.loc[bundled.index, bundled.columns].to_numpy() - bundled.to_numpy()))) if same else None
        record = {"accession": f"OSD-{n}", "official_filename": metadata["file_name"],
                  "download_url": url, "file_listing_url": listing_url,
                  "retrieval_date": datetime.datetime.now(datetime.timezone.utc).date().isoformat(), "download_sha256": hashlib.sha256(raw).hexdigest(),
                  "download_bytes": len(raw), "bundled_upstream_path": target.name,
                  "bundled_upstream_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                  "equal_gene_and_sample_ids": same, "max_count_difference": difference,
                  "n_genes": len(bundled), "n_samples": len(bundled.columns)}
        records.append(record)
        if not same or difference != 0:
            raise ValueError(f"OSD-{n} local counts do not match NASA STAR matrix: {record}")
    (HERE / "inputs/upstream_count_provenance.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
