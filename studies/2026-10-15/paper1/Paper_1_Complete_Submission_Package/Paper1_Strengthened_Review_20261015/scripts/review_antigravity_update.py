#!/usr/bin/env python3
"""Read-only checks of the Antigravity Paper 1 update; fail on known blockers.

This is a focused review gate, not a replacement for statistical validation.
It writes only to the separate correction workspace and never edits its source.
"""
import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    revision = args.root / "research/paper1_revision_2026"
    source = args.root / "research/q1_flagship_2026"
    output = args.output or revision / "antigravity_update_review_20261014"
    output.mkdir(parents=True, exist_ok=True)
    manuscript = source / "PAPER_1_MANUSCRIPT.md"
    text = manuscript.read_text()
    lines = text.splitlines()
    authors = json.loads((revision / "authors.json").read_text())
    with (source / "tables/H25_doseresponse_regressions.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    primary = next(row for row in rows if row["subset"] == "Primary_SingleTissue_25Contrasts"
                   and row["outcome"] == "IEGHSP26_abs_mean_t")
    abstract = text.split("## Abstract", 1)[1].split("**Keywords:", 1)[0]
    section = text.split("#### 2.3.2.", 1)[1].split("### 2.4.", 1)[0]
    checks = []

    def check(name, passed, observed, expected, needle):
        line = next((i for i, value in enumerate(lines, 1) if needle in value), None)
        checks.append(dict(check=name, passed=bool(passed), observed=observed,
                           expected=expected, source=str(manuscript), line=line))

    expected_names = [author["name"] for author in authors["authors"]]
    check("Supplied author list", all(name in text[:2000] for name in expected_names),
          lines[2], expected_names, "**Zainulabedin")
    check("Supplied affiliation and correspondence",
          authors["affiliation"] in text[:2000] and authors["corresponding_email"] in text[:2000],
          lines[4:6], [authors["affiliation"], authors["corresponding_email"]], "Independent Computational")
    unsupported = [claim for claim in ["received no external grant funding",
                   "All authors read and approved the final manuscript",
                   "declare no competing financial or non-financial interests"] if claim in text]
    check("Unresolved declarations stay unresolved", not unsupported, unsupported,
          "Funding, interests and author approval are unresolved in supplied author record", "## Acknowledgments")
    check("Abstract contrast count", "22 cross-deposit contrasts" not in abstract,
          "22 contrasts" if "22 cross-deposit contrasts" in abstract else "No stale 22 claim",
          {"primary": int(primary["n_contrasts"]), "all": 29}, "Across 22 cross-deposit")
    check("Abstract model specification", "eight-predictor" not in abstract,
          "eight-predictor" if "eight-predictor" in abstract else "No stale eight-predictor claim",
          "Two predictors plus intercept: cohort indicator and duration_days", "eight-predictor")
    check("Primary regression count in narrative and Table 4", not re.search(r"n\s*=\s*27|27 contrasts", section),
          sorted(set(re.findall(r"n\s*=\s*27|27 contrasts", section))),
          int(primary["n_contrasts"]), "Across the 27 contrasts")
    expected_r = f"{float(primary['r_cohort_confounded']):.4f}"
    primary_paragraph = next((value for value in section.splitlines()
                              if value.startswith("Across ") and "primary regression subset" in value), "")
    check("Primary regression coefficient in its own paragraph", expected_r in primary_paragraph,
          re.findall(r"r = [0-9.]+", primary_paragraph), expected_r, "Across the 27 contrasts")
    check("Comparator is not asserted to be unconfounded", "unconfounded `FLT vs VGC`" not in text
          and "unconfounded `FLT vs. VGC`" not in text,
          "Unconfounded FLT versus VGC language persists", "Housing/date comparator; biological truth unknown",
          "unconfounded `FLT vs. VGC`")
    check("Figure 1 tissue labels", "cerebellum, cortex, and brainstem" not in text,
          "RR-10 caption uses cortex/brainstem" if "cerebellum, cortex, and brainstem" in text else "Absent",
          "OSD-563 cerebellum; OSD-564 hippocampus; OSD-612 left cerebral hemisphere", "total plotted profile bars")
    check("Sensor induction timings are not assigned as measured harvest time",
          "diagnostic of an acute 15–60 min perimortem window" not in text,
          "Graphical abstract assigns a 15–60 min window", "No measured transcriptomic stopwatch",
          "diagnostic of an acute")
    check("Finite random-set p-value resolution acknowledged", "1 / 2001" in text,
          "1/2001 stated", "Empirical resolution floor for 2000 random sets", "resolution floor")
    check("Secondary responses not universally flat", "SRG` was not universally flat" in text,
          "Explicit qualification", "Retain the hippocampal secondary-response result", "not universally flat")
    check("Split partitions not independent replication", "statistically dependent partitions" in text,
          "Explicit dependence statement", "126 partitions reuse one cohort", "statistically dependent partitions")
    check("Anchor uncertainty limitation acknowledged", "without propagating estimation covariance" in text,
          "Plug-in uncertainty caveat", "Retain caveat and perform joint inference for stronger claims",
          "without propagating estimation covariance")
    inputs = [manuscript, source / "PAPER_1_MANUSCRIPT.tex",
              source / "tables/H25_doseresponse_regressions.csv",
              source / "scripts/run_H25_cross_deposit_doseresponse.py",
              source / "scripts/run_deep_analyses_p1_p2.py",
              source / "scripts/verify_manuscript_table_sync.py",
              source / "tables/manuscript_table_sync_audit.json"]
    report = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(), scope="Focused Paper 1 update review",
                  source_checksums={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  checks=checks, passed=sum(item["passed"] for item in checks),
                  failed=sum(not item["passed"] for item in checks), publication_ready=False,
                  primary_regression_row=primary,
                  limitation="These checks identify selected known blockers; they do not certify all tables or claims.")
    (output / "focused_review_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Focused checks: {report['passed']} passed; {report['failed']} failed. Publication ready: false.")
    print(output / "focused_review_checks.json")
    return 1 if report["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
