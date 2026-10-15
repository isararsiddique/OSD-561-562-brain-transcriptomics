"""Worker A, step 1: extract every number from Paper 1 and match it to the package result files.

Scope: MANUSCRIPT_COMPLETE.md Abstract, Introduction, Results, Discussion, Tables 1-6 and
figure legends; SUPPLEMENTARY_INFORMATION.md results text and figure legends (not the file
index Table S1).  Identifiers (accessions, mission names, gene symbols, strain names,
figure/table/equation/reference numbers) are removed before extraction.

Matching: every numeric cell of every CSV/JSON results file in results/ and extensions/
(excluding files that are themselves manuscript text or numeric look-up tables, and the
gene-level/repetition tables above 3 MB) is indexed.  A written number matches a source value
if the source rounds to it at the written precision (MATCH), or differs by at most one unit in
the last written digit (ROUNDING).  Derived forms are also tried: percent of a fraction,
100*(1 - ratio), sign flip.  Numbers with no indexed source are NO_SOURCE; the recomputation
(A_recompute_results.json) and hand checks then resolve them.  Hand-check verdicts are
applied from HAND (bottom) and override the automatic status.

Output: A_claims_check.csv with id, section, sentence, number_as_written, source_file,
source_value, status, note.
"""
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = Path(__file__).resolve().parent

SUP = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
EXCLUDE_FILES = {"manuscript_content.json", "manuscript_numeric_lookups.csv", "manuscript_sources.json",
                 "supplement_content.json", "manuscript_pre_strengthening.json", "supplement_pre_strengthening.json",
                 "validation_report.json", "strengthening_validation.json", "publication_figure_manifest.json",
                 "figure_manifest.json", "visual_review.json", "journal_artwork_validation.json",
                 "journal_format_validation.json", "supplementary_data_validation.json", "PACKAGE_FILE_MANIFEST.json",
                 "reference_number_map.json", "strengthening_reference_number_map.json",
                 "strengthening_equation_number_map.json", "manuscript_format_snapshot.json"}


# ----------------------------------------------------------------------------- text

def sections_md(path, keep_top, sup=False):
    text = path.read_text()
    out, top, sub = [], None, None
    for line in text.splitlines():
        if line.startswith("## "):
            top, sub = line[3:].strip(), None
            continue
        if line.startswith("### "):
            sub = line[4:].strip()
            continue
        if line.startswith("# ") or not line.strip() or line.startswith("!["):
            continue
        if top in keep_top:
            if sup and line.startswith("|") and ("extensions/" in line or "results/" in line or "scripts/" in line or "Files |" in line):
                continue                        # Table S1 is a file index, not results
            if line.startswith("| ---"):
                continue
            sec = top + (f" / {sub}" if sub else "")
            if top == "Figures":
                m = re.match(r"(Table \d+|Figure \d+)\.", line)
                if m:
                    sec = f"Figures / {m.group(1)}"
                    cur_tab = m.group(1)
                elif line.startswith("|"):
                    sec = f"Figures / {cur_tab}"
            if sup and re.match(r"(Supplementary Figure \d+)\.", line):
                sec = "Supplementary / " + re.match(r"(Supplementary Figure \d+)\.", line).group(1)
            out.append((("SI: " if sup else "") + sec, line))
    return out


def split_sentences(par):
    if par.startswith("|"):
        cells = [c.strip() for c in par.strip("|").split("|")]
        return [" | ".join(cells)]
    # protect decimals and abbreviations
    p = re.sub(r"(Fig|Figs|Eq|Eqs|Supplementary Fig|et al|e\.g|i\.e|vs)\.", lambda m: m.group(0).replace(".", "§"), par)
    parts = re.split(r"(?<=[.;])\s+(?=[A-Z(\u2212−0-9])", p)
    return [s.replace("§", ".").strip() for s in parts if s.strip()]


STRIP = [
    r"\[[\d,\s–\-]+\]",                                  # citations
    r"(Supplementary )?(Fig|Figs|Figure|Figures|Table|Tables|Data|Eq|Eqs|equation)\.?\s*\d+[a-z]?(\s*(,|–|-|and)\s*\d+[a-z]?)*",
    r"\b(OSD|GLDS|RR|RRRM|SpX|SpaceX)-\d+[A-Za-z]*",
    r"\bRR-\d+\b", r"C57BL/6\w*", r"B6129SF2/J", r"\bBALB/c", r"\b[Cc]obalt-57\b", r"\bSi-?28\b", r"\bCo-57\b",
    r"\blog[2₂]\b", r"\blog₂\(", r"\bL=2\b(?= before)",            # keep other L=2 mentions
    r"\(\d+\)\.?\s*$",                                    # equation numbers
    r"\b(release|version)\s+\d+", r"\bPyDESeq2 \d+\.\d+\.\d+", r"\bneuroCombat \d+\.\d+\.\d+",
    r"[A-Za-z]+\d+[A-Za-z0-9]*",                          # gene-like or code tokens: Egr1, Sensor31, HC3, PLS2, age29
    r"\b[A-Z][a-z]{0,3}-[a-z]{1,2}\d*\b",                 # Hbb-bt, Hba-a1
    r"Rep\d+", r"\b(one|two)-",
]
NUM = re.compile(r"(?<![A-Za-z0-9_/.\-])([−\-]?\d{1,3}(?:,\d{3})+(?:\.\d+)?|[−\-]?\d+(?:\.\d+)?)"
                 r"(\s*×\s*10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))?(\s*%)?")


def numbers(sentence):
    s = sentence
    for pat in STRIP:
        s = re.sub(pat, lambda m: " " * len(m.group(0)), s)
    out = []
    for m in NUM.finditer(s):
        raw = m.group(0).strip()
        base = m.group(1).replace("−", "-").replace(",", "")
        # A hyphen that directly follows a letter/space-less token is a range dash, not a minus
        if base.startswith("-") and m.start() > 0 and sentence[m.start() - 1].isdigit():
            base = base[1:]
        exp = m.group(3).translate(SUP) if m.group(3) else None
        try:
            val = float(base) * (10 ** int(exp) if exp else 1)
        except ValueError:
            continue
        dec = len(base.split(".")[1]) if "." in base else 0
        out.append(dict(raw=raw, value=val, dec=dec, exp=int(exp) if exp else None, pct=bool(m.group(4)),
                        integer="." not in base and not exp))
    return out


# ----------------------------------------------------------------------------- source index

def iter_tables():
    for root in ("results", "extensions"):
        for f in sorted((PKG / root).rglob("*")):
            if not f.is_file() or f.name in EXCLUDE_FILES:
                continue
            if "raw_sources" in f.parts or "literature" in f.parts or "journal_audit" in str(f) or "cache" in f.parts:
                continue
            if f.suffix not in (".csv", ".json") and not f.name.endswith(".csv.gz"):
                continue
            if f.stat().st_size > 3_000_000:
                continue
            yield f


def flatten_json(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten_json(v, f"{prefix}.{k}" if prefix else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten_json(v, f"{prefix}[{i}]")
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        yield prefix, float(obj)


def build_index():
    """Per-file sorted value arrays; cached in $KIROCREW_SCRATCH to avoid re-reading 6M cells."""
    import os, pickle
    cache = Path(os.environ.get("KIROCREW_SCRATCH", str(OUT))) / "A_claims_value_index.pkl"
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    index = {}
    for f in iter_tables():
        rel = str(f.relative_to(PKG))
        vals, locs = [], []
        try:
            if f.suffix == ".json":
                for k, v in flatten_json(json.loads(f.read_text())):
                    vals.append(v)
                    locs.append(k)
            else:
                df = pd.read_csv(f, low_memory=False)
                keycols = [c for c in df.columns if df[c].dtype == object][:4]
                keys = df[keycols].astype(str).agg("|".join, axis=1).to_numpy() if keycols else np.array([""] * len(df))
                num = df.select_dtypes(include=[np.number])
                arr = num.to_numpy(float)
                for j, col in enumerate(num.columns):
                    ok = np.flatnonzero(np.isfinite(arr[:, j]))
                    vals.extend(arr[ok, j].tolist())
                    locs.extend([f"row{i}[{keys[i]}]::{col}" for i in ok])
        except Exception as e:
            print("skip", rel, e)
            continue
        if vals:
            v = np.array(vals)
            o = np.argsort(v)
            index[rel] = (v[o], np.array(locs, dtype=object)[o])
    cache.write_bytes(pickle.dumps(index))
    return index


# ----------------------------------------------------------------------------- matching

PREF = {  # section keywords -> preferred files (substrings), checked in order
    "Sample identity": ["dataset_accounting", "animal_region_manifest", "metadata_conflicts", "condition_structure_R2", "real_data_structure"],
    "External datasets": ["external_data", "biological_validation", "cohort_inventory"],
    "Ground collection cohorts": ["panel_results", "gradient_results"],
    "Count models": ["count_model_panel_comparison", "panel_results", "count_model"],
    "direct flight by collection": ["flight_collection_interaction"],
    "Pairing and reproducibility": ["paired_biology", "regional_gene_replication", "saturated_age"],
    "Score removal": ["operator_geometry", "mathematics_review"],
    "Collection-score discrimination": ["ground_holdout", "exact_ground", "independent_audit"],
    "Constrained multi-sensor": ["method_benchmark/strengthened", "method_benchmark/matched_grid", "method_benchmark/method_summary", "convex"],
    "generalizes": ["method_benchmark/strengthened", "method_benchmark/matched_grid"],
    "Real-data attenuation": ["ground_anchor_sensitivity", "joint_bootstrap_summary", "OSD525_disjoint"],
    "Known-truth": ["known_truth_benchmark", "high_gain_benchmark", "fixed_sensor_noise"],
    "Count simulations": ["correction_benchmark"],
    "Expanded robustness": ["biological_validation/all_biological_effects", "biological_validation"],
    "Table 1": ["dataset_accounting"], "Table 3": ["strengthened_summary"], "Table 4": ["known_truth_benchmark_config"],
    "Table 5": ["high_gain_benchmark_config"],
    "Electronic results": ["flight_collection_interaction"], "Ground holdout": ["ground_holdout", "real_data_structure"],
    "Supplementary results": ["independent_audit", "training_expression_tertile"],
}


def search(index, files, target, tol):
    out = []
    for f in files:
        v, l = index[f]
        lo, hi = np.searchsorted(v, target - tol), np.searchsorted(v, target + tol, side="right")
        out += [(v[i], f"{f}::{l[i]}") for i in range(lo, hi)]
    return out


def preferred_files(index, section):
    prefs = []
    for k, files in PREF.items():
        if k.lower() in section.lower():
            prefs += files
    chosen = []
    for p in prefs:
        chosen += [f for f in index if p in f and f not in chosen]
    return chosen


def row_of(loc):
    # file::rowN[...]::col  -> file::rowN ;  json file::a.b[3].c -> file::a.b[3]
    parts = loc.split("::")
    if len(parts) >= 3 and parts[1].startswith("row"):
        return parts[0] + "::" + parts[1].split("[")[0]
    key = parts[1] if len(parts) > 1 else ""
    return parts[0] + "::" + (key.rsplit(".", 1)[0] if "." in key else "")


def all_candidates(index, n, section):
    """All candidates (status, how, value, loc, preferred) for one written number, capped."""
    unit = 10.0 ** (-n["dec"]) if not n["exp"] else abs(n["value"]) * 10 ** (-n["dec"])
    tries = [("direct", n["value"], unit)]
    if n["pct"]:
        tries += [("fraction", n["value"] / 100, unit / 100), ("1-ratio", 1 - n["value"] / 100, unit / 100),
                  ("ratio-1", 1 + n["value"] / 100, unit / 100)]
    tries += [("sign-flipped", -n["value"], unit)]
    pref = preferred_files(index, section)
    out = []
    for tolmult, status in ((0.5000001, "MATCH"), (1.5, "ROUNDING")):
        for how, t, u in tries:
            for v, l in search(index, pref, t, tolmult * u)[:5000]:
                out.append((status, how, v, l, True))
        if out:
            return out
    sig = len(re.sub(r"[^0-9]", "", n["raw"].split("×")[0]).lstrip("0"))
    if not n["integer"] and sig >= 3:
        for v, l in search(index, list(index), n["value"], 0.5000001 * unit)[:5000]:
            out.append(("MATCH", "direct (outside the section's expected files)", v, l, False))
    return out


TOKEN = re.compile(r"[a-z][a-z0-9_]+|osd-\d+")
SYN = {"cerebellum": ["cerebellum", "cb", "561"], "hippocampus": ["hippocampus", "hip", "562"], "circadian": ["circadian"],
       "rapid": ["rprg"], "secondary": ["srg"], "heat-shock": ["hsp"], "ridge": ["ridge"], "constrained": ["constrained"],
       "bounded": ["bounded"], "projection": ["projection"], "raw": ["raw"], "score": ["score"], "flight": ["flt"],
       "ground": ["gc", "ground"], "isS-t": ["isst"], "lar": ["lar"]}


def sentence_tokens(sent):
    t = set(TOKEN.findall(sent.lower()))
    for k, v in SYN.items():
        if k.lower() in sent.lower():
            t |= set(v)
    return t


def assign_sentence(cands_per_number, sent=""):
    """Pick, for each number, the candidate whose source row covers the most numbers of the sentence,
    then whose row key shares most tokens with the sentence; prefer summary tables over
    per-repetition, development-archive and gene-level tables."""
    cover = {}
    for k, cands in enumerate(cands_per_number):
        for c in {row_of(c[3]) for c in cands}:
            cover.setdefault(c, set()).add(k)
    toks = sentence_tokens(sent)

    def score(c):
        loc = c[3].lower()
        key_t = set(TOKEN.findall(loc.replace("_", " ")))
        overlap = len(toks & key_t)
        pen = sum(w in loc for w in ("repetition", "development_v1", "gene_results", "nb_interaction", "nb_coefficients",
                                     "heldout_scores", "inner_cv", "unit_errors", "counts.csv", "endpoint_scores", "split_results"))
        bonus = sum(w in loc for w in ("summary", "panel_results", "interactions.csv", "operator_geometry", "accounting"))
        return (c[0] == "MATCH", -pen, len(cover[row_of(c[3])]), overlap + 2 * bonus, c[1] == "direct")
    picks = []
    for k, cands in enumerate(cands_per_number):
        if not cands:
            picks.append((None, 0))
            continue
        best = max(cands, key=score)
        picks.append((best, len(cover[row_of(best[3])])))
    return picks


def main():
    pars = sections_md(PKG / "MANUSCRIPT_COMPLETE.md", {"Abstract", "Introduction", "Results", "Discussion", "Figures"})
    pars += sections_md(PKG / "SUPPLEMENTARY_INFORMATION.md", {"Supplementary results"}, sup=True)
    index = build_index()
    print("indexed files:", len(index), "values:", sum(len(v[0]) for v in index.values()))
    recompute = json.loads((OUT / "A_recompute_results.json").read_text()) if (OUT / "A_recompute_results.json").exists() else []
    rows = []
    k = 0
    for section, par in pars:
        for sent in split_sentences(par):
            nums = numbers(sent)
            cands = [all_candidates(index, n, section) for n in nums]
            picks = assign_sentence(cands, sent)
            for n, (b, cov) in zip(nums, picks):
                k += 1
                sig = len(re.sub(r"[^0-9]", "", n["raw"].split("×")[0]).lstrip("0"))
                src_file, src_val, note = "", "", ""
                if b is None:
                    status = "NO_SOURCE"
                else:
                    status, how, v, l, preferred = b
                    src_file, src_val = l, v
                    note = f"auto: {how}; source row also holds {cov - 1} other number(s) of this sentence"
                    weak = (n["integer"] or sig <= 2) and cov < 2
                    if weak or (n["integer"] and how != "direct"):
                        status, note = "NO_SOURCE", "auto: only an isolated low-precision/integer coincidence found; hand check needed"
                        src_file, src_val = "", ""
                rows.append(dict(id=f"N{k:04d}", section=section, sentence=sent, number_as_written=n["raw"],
                                 parsed_value=n["value"], integer=n["integer"], source_file=src_file, source_value=src_val,
                                 status=status, note=note))
    df = pd.DataFrame(rows)
    # hand-check verdicts (applied after the automatic pass; see A_claims_hand.py)
    df["source_value"] = df.source_value.astype(object)
    df["auto_status"] = df.status
    import A_claims_hand
    n = A_claims_hand.apply(df)
    print("hand-check rule applications:", n)
    df.drop(columns=["parsed_value", "integer"]).to_csv(OUT / "A_claims_check.csv", index=False)
    print(df.groupby(["auto_status", "status"]).size().to_string())
    print(df.status.value_counts().to_string())
    return df


import sys
sys.path.insert(0, str(OUT))

if __name__ == "__main__":
    main()
