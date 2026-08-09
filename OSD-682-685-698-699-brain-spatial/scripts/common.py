#!/usr/bin/env python3
"""Shared IO, design and statistics helpers for the OSD-682/685/698/699 re-analysis.

The four OSD accessions are four brain regions of a *single* GeoMx Digital Spatial
Profiling experiment (GEO series GSE239336):

    OSD-682  GLDS-613  CA   Cornu Ammonis 1 (CA1, hippocampus)
    OSD-685  GLDS-616  DG   Dentate Gyrus (hippocampus)
    OSD-698  GLDS-627  FCT  Frontal Cortex
    OSD-699  GLDS-628  CT   Cerebral Cortex

Design (per region): 2 x 2 factorial, n = 3
    Spaceflight : Space Flight (FLT) vs Ground Control (GC)
    Treatment   : BuOE antioxidant (TX) vs saline control (SAL)

Nothing about the design is hard-coded from assumption: group membership, ROI
identity and every QC metric are parsed out of the deposited GEO series matrix
and the OSDR ISA archives.  The only hard-coded facts are the accession
cross-reference table above and the sign convention resolved empirically by
``scripts/check_conventions.py``.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
GEO_DIR = DATA / "geo"
OSDR_DIR = DATA / "osdr"
RESULTS = ROOT / "results"
FIG_DIR = RESULTS / "figures"
PANEL_DIR = FIG_DIR / "panels"
TAB_DIR = RESULTS / "tables"
CSV_DIR = TAB_DIR / "csv"
CONFIG = ROOT / "config"

GEO_SERIES = "GSE239336"

# --------------------------------------------------------------------------- #
# Study cross-reference
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Study:
    osd: str
    glds: str
    code: str  # region code used in the GEO supplementary filenames
    region: str  # human-readable region label
    short: str  # compact label for figures
    geo_tissue: str  # value of the GEO `tissue:` characteristic


STUDIES: tuple[Study, ...] = (
    Study("OSD-682", "GLDS-613", "CA", "Cornu Ammonis 1", "CA1", "Cornu Ammonis 1"),
    Study("OSD-685", "GLDS-616", "DG", "Dentate Gyrus", "DG", "Dentate Gyrus"),
    Study("OSD-698", "GLDS-627", "FCT", "Frontal Cortex", "FCtx", "Frontal Cortex"),
    Study("OSD-699", "GLDS-628", "CT", "Cerebral Cortex", "Ctx", "Cortex"),
)

BY_CODE = {s.code: s for s in STUDIES}
BY_OSD = {s.osd: s for s in STUDIES}
REGION_CODES = [s.code for s in STUDIES]
REGION_ORDER = [s.short for s in STUDIES]

# --------------------------------------------------------------------------- #
# Experimental groups
# --------------------------------------------------------------------------- #
# Slide id -> (spaceflight level, treatment level).  On this GeoMx run each
# experimental group occupies its own slide, so slide id and group are
# equivalent (and perfectly confounded -- see README "Limitations").
SLIDE_TO_GROUP = {
    "FLT_SAL_2": ("FLT", "SAL"),
    "FLT_TX_2": ("FLT", "BuOE"),
    "GC_SAL_2": ("GC", "SAL"),
    "GC_TX_2": ("GC", "BuOE"),
}

GROUPS = ["GC_SAL", "GC_BuOE", "FLT_SAL", "FLT_BuOE"]
GROUP_LABELS = {
    "GC_SAL": "Ground / Saline",
    "GC_BuOE": "Ground / BuOE",
    "FLT_SAL": "Flight / Saline",
    "FLT_BuOE": "Flight / BuOE",
}
GROUP_COLORS = {
    "GC_SAL": "#4C72B0",
    "GC_BuOE": "#55A868",
    "FLT_SAL": "#C44E52",
    "FLT_BuOE": "#DD8452",
}

SF_LEVELS = ["GC", "FLT"]
TRT_LEVELS = ["SAL", "BuOE"]

# --------------------------------------------------------------------------- #
# Contrasts
# --------------------------------------------------------------------------- #
# Contrasts distributed with the GEO series (NanoString DSP Analysis Suite
# exports).  `numerator`/`denominator` record the orientation *after* the sign
# correction below, i.e. positive log2 = higher in `numerator`.
#
# The deposited "Log2" column orientation was resolved empirically against the
# Q3-normalised matrix by scripts/check_conventions.py; SIGN_FLIP records
# whether the deposited column has to be negated to match the stated
# numerator - denominator orientation.
PROVIDED_CONTRASTS = {
    "GCvsFLT-SAL": dict(
        label="Flight vs Ground (Saline)",
        numerator="FLT_SAL",
        denominator="GC_SAL",
        factor="Spaceflight",
        stratum="Saline",
    ),
    "SALvsBuOE-GC": dict(
        label="BuOE vs Saline (Ground)",
        numerator="GC_BuOE",
        denominator="GC_SAL",
        factor="Treatment",
        stratum="Ground",
    ),
    "SALvsBuOE-FLT": dict(
        label="BuOE vs Saline (Flight)",
        numerator="FLT_BuOE",
        denominator="FLT_SAL",
        factor="Treatment",
        stratum="Flight",
    ),
}

# Filled in from config/sign_convention.json (written by check_conventions.py).
_SIGN_CACHE: dict[str, int] | None = None

# Simple contrasts re-derived from the Q3 matrix (numerator, denominator).
DERIVED_CONTRASTS = {
    "FLT_vs_GC_in_SAL": ("FLT_SAL", "GC_SAL", "Spaceflight", "Saline"),
    "FLT_vs_GC_in_BuOE": ("FLT_BuOE", "GC_BuOE", "Spaceflight", "BuOE"),
    "BuOE_vs_SAL_in_GC": ("GC_BuOE", "GC_SAL", "Treatment", "Ground"),
    "BuOE_vs_SAL_in_FLT": ("FLT_BuOE", "FLT_SAL", "Treatment", "Flight"),
}

# DEG thresholds.  Matched to the sibling OSD-561/562 project for comparability.
FDR_CUTOFF = 0.05
LFC_CUTOFF = 1.0
# GeoMx DSP is a targeted 15.8k-probe assay on 3 ROIs/group; effect sizes are
# compressed relative to bulk RNA-seq, so a relaxed secondary threshold is also
# reported.
LFC_CUTOFF_RELAXED = 0.5


def sign_for(contrast: str) -> int:
    """Return +1 or -1: the multiplier applied to the deposited Log2 column."""
    global _SIGN_CACHE
    if _SIGN_CACHE is None:
        path = CONFIG / "sign_convention.json"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing -- run scripts/check_conventions.py first "
                "(it resolves the deposited log2 orientation empirically)."
            )
        _SIGN_CACHE = json.loads(path.read_text())["signs"]
    return int(_SIGN_CACHE[contrast])


# --------------------------------------------------------------------------- #
# Low-level readers
# --------------------------------------------------------------------------- #
def _read_text(path: Path) -> list[str]:
    """Read a possibly gzipped text file and normalise CR / CRLF line endings."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    return raw.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n").split("\n")


def q3_path(code: str) -> Path:
    return GEO_DIR / f"{GEO_SERIES}_{code}_GeneExpression_Q3norm.txt.gz"


def de_path(code: str, contrast: str) -> Path:
    return GEO_DIR / f"{GEO_SERIES}_{code}_{contrast}_DEanalysis.txt.gz"


def read_q3(code: str) -> pd.DataFrame:
    """Q3-normalised GeoMx expression matrix: genes x ROIs.

    Column names in the deposited file look like ``FLT_TX_2 | 002 | Full ROI``;
    they are rewritten to ``<slide>|<roi>`` keys that join onto the sample table.
    """
    lines = _read_text(q3_path(code))
    header = lines[0].split("\t")
    cols = []
    for c in header[1:]:
        slide, roi, _segment = [p.strip() for p in c.split("|")]
        cols.append(f"{slide}|{roi}")
    genes: list[str] = []
    rows: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        genes.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    df = pd.DataFrame(rows, index=pd.Index(genes, name="TargetName"), columns=cols)
    return df


def read_provided_de(code: str, contrast: str) -> pd.DataFrame:
    """Parse a NanoString DSP Analysis Suite differential-expression export.

    The file carries a short key/value preamble before the real header row,
    which starts with ``Target tag``.
    """
    lines = _read_text(de_path(code, contrast))
    hidx = next(i for i, l in enumerate(lines) if l.startswith("Target tag"))
    meta = {}
    for l in lines[:hidx]:
        parts = l.split("\t")
        if len(parts) >= 2 and parts[0]:
            meta[parts[0].lstrip("#").strip()] = parts[1].strip()
    header = [h.strip() for h in lines[hidx].split("\t")]
    body = [l for l in lines[hidx + 1 :] if l.strip()]
    reader = csv.reader(io.StringIO("\n".join(body)), delimiter="\t", quotechar='"')
    df = pd.DataFrame(list(reader), columns=header)

    out = pd.DataFrame(
        {
            "gene": df["Target name"].astype(str),
            "target_group": df["Target group membership/s"].astype(str),
            "log2FC_deposited": pd.to_numeric(df["Log2"], errors="coerce"),
            "pvalue": pd.to_numeric(df["Pvalue"], errors="coerce"),
            "padj_deposited": pd.to_numeric(df["Adjusted pvalue"], errors="coerce"),
        }
    )
    out.attrs["dataset_name"] = meta.get("Dataset Name", "")
    out.attrs["analysis_name"] = meta.get("Analysis Name", "")
    return out


def read_series_matrix() -> pd.DataFrame:
    """Per-ROI design + QC metadata for all 48 samples, from the GEO series matrix."""
    path = GEO_DIR / f"{GEO_SERIES}_series_matrix.txt.gz"
    lines = _read_text(path)
    fields: dict[str, list[list[str]]] = {}
    for line in lines:
        if not line.startswith("!Sample_"):
            continue
        parts = next(csv.reader(io.StringIO(line), delimiter="\t", quotechar='"'))
        fields.setdefault(parts[0], []).append(parts[1:])

    titles = fields["!Sample_title"][0]
    accessions = fields["!Sample_geo_accession"][0]
    n = len(titles)

    # `treatment:` appears twice per sample (spaceflight level, then drug level);
    # keep both by suffixing repeats.
    per: list[dict[str, str]] = [{} for _ in range(n)]
    for row in fields.get("!Sample_characteristics_ch1", []):
        for i, value in enumerate(row):
            if not value or ":" not in value:
                continue
            key, val = value.split(":", 1)
            key, val = key.strip(), val.strip()
            base = key
            k = 2
            while key in per[i]:
                key = f"{base}_{k}"
                k += 1
            per[i][key] = val

    recs = []
    for i in range(n):
        c = per[i]
        slide = c["slide id"]
        roi = c["roi number"]
        sf, trt = SLIDE_TO_GROUP[slide]
        title = titles[i]
        m = re.search(r"_rep(\d+)$", title)
        recs.append(
            {
                "sample": accessions[i],
                "title": title,
                "roi_key": f"{slide}|{roi}",
                "slide": slide,
                "roi": roi,
                "segment": c.get("segment type", ""),
                "geo_tissue": c.get("tissue", ""),
                "spaceflight": sf,
                "treatment": trt,
                "group": f"{sf}_{trt}",
                "rep": int(m.group(1)) if m else np.nan,
                "strain": c.get("strain", ""),
                "sex": c.get("Sex", c.get("sex", "")),
                "area_um2": pd.to_numeric(c.get("area"), errors="coerce"),
                "nuclei": pd.to_numeric(c.get("nuclei_counts"), errors="coerce"),
                "roi_x": pd.to_numeric(c.get("roi x coordinate"), errors="coerce"),
                "roi_y": pd.to_numeric(c.get("roi y coordinate"), errors="coerce"),
                "raw_reads": pd.to_numeric(c.get("rawreads"), errors="coerce"),
                "trimmed_reads": pd.to_numeric(c.get("trimmedreads"), errors="coerce"),
                "stitched_reads": pd.to_numeric(c.get("stitchedreads"), errors="coerce"),
                "aligned_reads": pd.to_numeric(c.get("alignedreads"), errors="coerce"),
                "dedup_reads": pd.to_numeric(c.get("deduplicatedreads"), errors="coerce"),
                "seq_saturation": pd.to_numeric(c.get("sequencingsaturation"), errors="coerce"),
            }
        )
    meta = pd.DataFrame(recs)

    tissue_to_code = {s.geo_tissue: s.code for s in STUDIES}
    meta["region_code"] = meta["geo_tissue"].map(tissue_to_code)
    unknown = meta.loc[meta["region_code"].isna(), "geo_tissue"].unique()
    if len(unknown):
        raise ValueError(f"unmapped GEO tissue label(s): {list(unknown)}")
    meta["region"] = meta["region_code"].map(lambda c: BY_CODE[c].region)
    meta["region_short"] = meta["region_code"].map(lambda c: BY_CODE[c].short)
    meta["osd"] = meta["region_code"].map(lambda c: BY_CODE[c].osd)
    meta["glds"] = meta["region_code"].map(lambda c: BY_CODE[c].glds)
    meta["align_rate"] = meta["aligned_reads"] / meta["raw_reads"]
    meta["dedup_rate"] = meta["dedup_reads"] / meta["aligned_reads"]

    order = {c: i for i, c in enumerate(REGION_CODES)}
    meta = meta.sort_values(
        by=["region_code", "group", "rep"],
        key=lambda s: s.map(order) if s.name == "region_code" else s,
    ).reset_index(drop=True)
    return meta


def read_isa_samples() -> pd.DataFrame:
    """Factor values and sequencing parameters from the four OSDR ISA archives."""
    recs = []
    for study in STUDIES:
        sdir = OSDR_DIR / study.osd
        sfile = sdir / f"s_{study.osd}.txt"
        if not sfile.exists():
            continue
        with open(sfile, newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        hdr = rows[0]

        def col(name: str) -> int | None:
            return hdr.index(name) if name in hdr else None

        i_sample = col("Sample Name")
        i_trt = col("Factor Value[Treatment]")
        i_sf = col("Factor Value[Spaceflight]")
        i_age = col("Characteristics[Age at Launch]")
        i_dur = col("Parameter Value[duration]")
        i_sex = col("Characteristics[Sex]")
        i_strain = col("Characteristics[Strain]")

        assay = next(sdir.glob("a_*.txt"), None)
        depth: dict[str, str] = {}
        rrna: dict[str, str] = {}
        raw_files: dict[str, str] = {}
        if assay is not None:
            with open(assay, newline="", encoding="utf-8-sig") as fh:
                arows = list(csv.reader(fh, delimiter="\t"))
            ah = arows[0]
            a_sample = ah.index("Sample Name")
            a_depth = ah.index("Parameter Value[Read Depth]") if "Parameter Value[Read Depth]" in ah else None
            a_rrna = ah.index("Parameter Value[rRNA Contamination]") if "Parameter Value[rRNA Contamination]" in ah else None
            a_raw = ah.index("Raw Data File") if "Raw Data File" in ah else None
            for r in arows[1:]:
                if not r or len(r) <= a_sample:
                    continue
                s = r[a_sample]
                if a_depth is not None and len(r) > a_depth:
                    depth[s] = r[a_depth]
                if a_rrna is not None and len(r) > a_rrna:
                    rrna[s] = r[a_rrna]
                if a_raw is not None and len(r) > a_raw:
                    raw_files[s] = r[a_raw]

        for r in rows[1:]:
            if not r or len(r) <= i_sample:
                continue
            s = r[i_sample]
            recs.append(
                {
                    "sample": s,
                    "osd": study.osd,
                    "glds": study.glds,
                    "region": study.region,
                    "isa_spaceflight": r[i_sf] if i_sf is not None else "",
                    "isa_treatment": r[i_trt] if i_trt is not None else "",
                    "isa_sex": r[i_sex] if i_sex is not None else "",
                    "isa_strain": r[i_strain] if i_strain is not None else "",
                    "age_at_launch_wk": pd.to_numeric(r[i_age], errors="coerce") if i_age is not None else np.nan,
                    "mission_duration_d": pd.to_numeric(r[i_dur], errors="coerce") if i_dur is not None else np.nan,
                    "isa_read_depth": pd.to_numeric(depth.get(s), errors="coerce"),
                    "isa_rrna_pct": pd.to_numeric(rrna.get(s), errors="coerce"),
                    "raw_data_files": raw_files.get(s, ""),
                }
            )
    return pd.DataFrame(recs)


def load_sample_table() -> pd.DataFrame:
    """GEO ROI metadata joined to the OSDR ISA factor values."""
    meta = read_series_matrix()
    isa = read_isa_samples()
    if not isa.empty:
        meta = meta.merge(
            isa.drop(columns=["osd", "glds", "region"]), on="sample", how="left"
        )
        # Cross-check: ISA factor values must agree with the GEO characteristics.
        sf_map = {"Space Flight": "FLT", "Ground Control": "GC"}
        trt_map = {"BuOE treated": "BuOE", "saline treated control": "SAL"}
        bad_sf = meta.loc[
            meta["isa_spaceflight"].notna()
            & (meta["isa_spaceflight"].map(sf_map) != meta["spaceflight"])
        ]
        bad_trt = meta.loc[
            meta["isa_treatment"].notna()
            & (meta["isa_treatment"].map(trt_map) != meta["treatment"])
        ]
        if len(bad_sf) or len(bad_trt):
            raise ValueError(
                "ISA factor values disagree with GEO characteristics for "
                f"{sorted(set(bad_sf['sample']) | set(bad_trt['sample']))}"
            )
    return meta


def load_expression(code: str, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (log2 Q3 matrix, per-ROI metadata) for one region, columns aligned."""
    mat = read_q3(code)
    sub = samples.loc[samples["region_code"] == code].copy()
    missing = set(sub["roi_key"]) - set(mat.columns)
    extra = set(mat.columns) - set(sub["roi_key"])
    if missing or extra:
        raise ValueError(
            f"{code}: ROI keys do not reconcile between the Q3 matrix and the "
            f"series matrix (missing={sorted(missing)}, extra={sorted(extra)})"
        )
    mat = mat[sub["roi_key"].tolist()]
    mat.columns = sub["sample"].tolist()
    floor = 1.0
    log2 = np.log2(mat.clip(lower=floor))
    return log2, sub.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs propagate."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = p[ok]
    n = vals.size
    order = np.argsort(vals, kind="mergesort")
    ranked = vals[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def two_way_anova(
    y: np.ndarray, sf: np.ndarray, trt: np.ndarray
) -> dict[str, np.ndarray]:
    """Per-gene 2x2 fixed-effects model on log2 expression.

    ``y`` is genes x samples.  ``sf`` / ``trt`` are the per-sample factor levels.
    Effect coding (-1 / +1) so that the estimated coefficients are half the
    difference between marginal means; they are doubled below to report
    interpretable log2 fold changes.

        y = b0 + b_sf * SF + b_trt * TRT + b_int * SF*TRT

    Returns log2 fold changes and t-test p-values for the two main effects and
    the interaction, plus the four simple (cell-mean) contrasts evaluated
    against the pooled residual variance of the full model.
    """
    y = np.asarray(y, dtype=float)
    x_sf = np.where(np.asarray(sf) == "FLT", 1.0, -1.0)
    x_trt = np.where(np.asarray(trt) == "BuOE", 1.0, -1.0)
    X = np.column_stack([np.ones_like(x_sf), x_sf, x_trt, x_sf * x_trt])
    n, p = X.shape
    df = n - p
    if df <= 0:
        raise ValueError(f"not enough samples ({n}) for the 2x2 model")

    XtX_inv = np.linalg.inv(X.T @ X)
    beta = y @ X @ XtX_inv.T  # genes x 4
    resid = y - beta @ X.T
    sigma2 = (resid ** 2).sum(axis=1) / df
    from scipy import stats

    res: dict[str, np.ndarray] = {}
    # Main effects and interaction. Coefficients are half-differences under
    # effect coding, so log2FC = 2 * beta.
    for name, j in (("spaceflight", 1), ("treatment", 2), ("interaction", 3)):
        se = np.sqrt(sigma2 * XtX_inv[j, j])
        with np.errstate(divide="ignore", invalid="ignore"):
            t = beta[:, j] / se
        res[f"log2FC_{name}"] = 2.0 * beta[:, j]
        res[f"se_{name}"] = 2.0 * se
        res[f"t_{name}"] = t
        res[f"p_{name}"] = 2.0 * stats.t.sf(np.abs(t), df)

    # Simple contrasts as linear combinations of the coefficients.
    simple = {
        # FLT - GC within saline:      (b0 + b_sf - b_trt - b_int) - (b0 - b_sf - b_trt + b_int)
        "FLT_vs_GC_in_SAL": np.array([0.0, 2.0, 0.0, -2.0]),
        "FLT_vs_GC_in_BuOE": np.array([0.0, 2.0, 0.0, 2.0]),
        "BuOE_vs_SAL_in_GC": np.array([0.0, 0.0, 2.0, -2.0]),
        "BuOE_vs_SAL_in_FLT": np.array([0.0, 0.0, 2.0, 2.0]),
    }
    for name, L in simple.items():
        est = beta @ L
        se = np.sqrt(sigma2 * (L @ XtX_inv @ L))
        with np.errstate(divide="ignore", invalid="ignore"):
            t = est / se
        res[f"log2FC_{name}"] = est
        res[f"se_{name}"] = se
        res[f"t_{name}"] = t
        res[f"p_{name}"] = 2.0 * stats.t.sf(np.abs(t), df)

    res["df_resid"] = np.full(y.shape[0], df, dtype=float)
    res["sigma2"] = sigma2
    res["mean_log2"] = y.mean(axis=1)
    return res


def call_degs(
    padj: np.ndarray, lfc: np.ndarray, fdr: float = FDR_CUTOFF, min_lfc: float = LFC_CUTOFF
) -> np.ndarray:
    padj = np.asarray(padj, dtype=float)
    lfc = np.asarray(lfc, dtype=float)
    return (padj < fdr) & (np.abs(lfc) >= min_lfc) & np.isfinite(padj) & np.isfinite(lfc)


def deg_counts(padj: np.ndarray, lfc: np.ndarray, **kw) -> dict[str, int]:
    sig = call_degs(padj, lfc, **kw)
    lfc = np.asarray(lfc, dtype=float)
    return {
        "n_deg": int(sig.sum()),
        "n_up": int((sig & (lfc > 0)).sum()),
        "n_down": int((sig & (lfc < 0)).sum()),
    }


# --------------------------------------------------------------------------- #
# Output helpers
# --------------------------------------------------------------------------- #
def ensure_dirs() -> None:
    for d in (RESULTS, FIG_DIR, PANEL_DIR, TAB_DIR, CSV_DIR):
        d.mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, name: str, subdir: Path | None = None) -> Path:
    target = (subdir or CSV_DIR)
    target.mkdir(parents=True, exist_ok=True)
    path = target / name
    df.to_csv(path, index=False)
    print(f"  [csv]  {path.relative_to(ROOT)}  ({len(df):,} rows)")
    return path
