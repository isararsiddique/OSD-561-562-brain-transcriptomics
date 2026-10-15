#!/usr/bin/env python3
"""Functional enrichment of Paper 2 junction gene sets via the g:Profiler REST API.

Background is the MEASURED gene universe (domain_scope='custom'), never the
genome: every gene in the background had at least one analysed junction cluster,
so the test asks whether flight- or housing-associated genes are unusual *among
genes this assay could see*.

Only public gene symbols leave the machine. Raw responses are written verbatim
to research/I_enrichment/raw/ so the result can be re-read without re-querying.
"""
from __future__ import annotations
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "tables"
RAW = HERE / "raw"
RAW.mkdir(parents=True, exist_ok=True)

API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC"]
UA = "Paper2-enrichment/1.0 (research use; contact via submitting author)"

QUERIES = {
    "FLT_GC_BY": "set_FLT_GC_BY.txt",
    "FLT_GC_top100": "set_FLT_GC_top100.txt",
    "GC_VIV_top100": "set_GC_VIV_top100.txt",
}


def read_list(name: str) -> list[str]:
    return [g for g in (TAB / name).read_text().split("\n") if g.strip()]


def post(payload: dict) -> dict:
    req = urllib.request.Request(
        API,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": UA},
    )
    with urllib.request.urlopen(req, timeout=300) as resp:
        return json.loads(resp.read().decode())


def main() -> None:
    background = read_list("background_measured_genes.txt")
    meta = {
        "api": API,
        "queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "organism": "mmusculus",
        "domain_scope": "custom",
        "significance_threshold_method": "g_SCS",
        "user_threshold": 0.05,
        "sources": SOURCES,
        "background_n": len(background),
        "runs": {},
    }
    frames = []
    for label, fname in QUERIES.items():
        genes = read_list(fname)
        payload = {
            "organism": "mmusculus",
            "query": {label: genes},
            "sources": SOURCES,
            "user_threshold": 0.05,
            "all_results": True,          # keep non-significant terms for honesty
            "ordered": False,
            "no_iea": False,
            "combined": False,
            "measure_underrepresentation": False,
            "domain_scope": "custom",
            "background": background,
            "significance_threshold_method": "g_SCS",
            "no_evidences": True,
        }
        res = post(payload)
        (RAW / f"gprofiler_{label}.json").write_text(json.dumps(res, indent=1))
        rows = pd.DataFrame(res.get("result", []))
        meta["runs"][label] = {
            "query_n_submitted": len(genes),
            "n_terms_returned": int(len(rows)),
            "n_significant": int(rows.significant.sum()) if len(rows) else 0,
            "effective_domain_size": (
                int(rows.effective_domain_size.iloc[0]) if len(rows) else None
            ),
            "query_size_used_by_api": int(rows.query_size.iloc[0]) if len(rows) else None,
        }
        if len(rows):
            rows = rows.assign(gene_set=label)
            frames.append(rows)
        print(label, meta["runs"][label])
        time.sleep(2)

    if frames:
        allr = pd.concat(frames, ignore_index=True)
        keep = [
            "gene_set", "source", "native", "name", "p_value", "significant",
            "term_size", "query_size", "intersection_size", "effective_domain_size",
            "precision", "recall",
        ]
        allr = allr[[c for c in keep if c in allr.columns]]
        # Observed/expected ratio under the measured background.
        allr["expected_overlap"] = (
            allr.term_size * allr.query_size / allr.effective_domain_size
        )
        allr["fold_enrichment"] = allr.intersection_size / allr.expected_overlap
        allr = allr.sort_values(["gene_set", "p_value"])
        allr.to_csv(TAB / "gprofiler_all_terms.csv.gz", index=False)
        allr[allr.significant].to_csv(TAB / "gprofiler_significant_terms.csv", index=False)
        print("\nsignificant terms per set:")
        print(allr[allr.significant].gene_set.value_counts().to_string())

    (RAW / "gprofiler_run_metadata.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
