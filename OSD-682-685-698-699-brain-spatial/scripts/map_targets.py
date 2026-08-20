#!/usr/bin/env python3
"""Map the CNS target panel (human) onto the mouse GeoMx target space.

Panel source: ``CNS mRNA targets.csv`` at the repository root, maintained by the
project lead. The list is human analytes with human Ensembl gene IDs, several
given as protein or assay names rather than gene symbols ("Amyloid-beta",
"GLT-1", "Tau", "Neurofilament Light (NfL)"). The data are mouse, so two things
must happen before any analysis and neither may be guessed:

  1. Resolve each Ensembl gene ID to its current human symbol. The ENSG is the
     stable identifier, so it -- not the spreadsheet label -- defines the gene.
     Any disagreement between the two is reported rather than silently resolved.
  2. Map each human gene to its mouse ortholog through Ensembl Compara, keeping
     the orthology type. Title-casing a human symbol is not a safe shortcut:
     human CXCL8 has no mouse ortholog at all, and several panel members were
     renamed between the two genomes.

Naming note: this panel is a CNS injury / neuroinflammation / neurovascular
target set. It is **not** an extracellular-vesicle panel -- the canonical EV and
exosome markers CD9, CD63 and CD81 are not part of it. They are tracked
separately as a reference set (EV_MARKERS below) purely so their measurability
can be reported; GeoMx measures tissue mRNA, which is not EV cargo.

Everything is cached under data/mapping/ so the network is hit once.

Outputs
  config/cns_targets_mapped.csv          full audit trail, one row per analyte
  config/cns_targets_mouse_symbols.txt   mouse symbols, assay-independent
  config/focus_genes.txt                 mouse symbols present on the GeoMx panel
  results/tables/csv/target_panel_mapping.csv
  results/tables/csv/target_panel_coverage_by_category.csv
  results/tables/csv/ev_marker_status.csv

Usage:  python scripts/map_targets.py [--offline]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

import common as C

PANEL_CSV = C.ROOT.parent / "CNS mRNA targets.csv"
CACHE = C.DATA / "mapping"
ENSEMBL = "https://rest.ensembl.org"
NEG_PROBE = "NegProbe-WTX"

# Canonical EV / exosome markers. Deliberately NOT part of the CNS target panel;
# carried only so their presence and behaviour can be reported separately.
EV_MARKERS = {
    "CD9": "ENSG00000010278",
    "CD63": "ENSG00000135404",
    "CD81": "ENSG00000110651",
}

# Functional grouping of the panel. Edit freely -- this only affects how results
# are aggregated and plotted, never which targets are tested.
CATEGORIES: dict[str, tuple[str, ...]] = {
    "Neuroinflammation / cytokine": (
        "CCL2", "CCL7", "CCR2", "CXCL10", "CXCL2", "CXCL8", "CXCR2", "CXCR3",
        "IL-17", "IL-18", "IL-1alpha", "IL11", "IL1B", "IL1RN", "IL22", "IL6",
        "IL6RA", "IL6ST", "TNF", "TNFAIP6", "TGFB1", "IRF1", "STAT3", "MYD88",
        "TLR4", "NFKB1", "NFKB2", "RELA", "PTGS2", "HMGB1", "AIM2", "C3",
    ),
    "Inflammasome / pyroptosis": (
        "NLRC4", "NLRP1", "NLRP2", "NLRP3", "PYCARD", "GSDMD", "ICE", "NINJ1",
    ),
    "Blood-brain barrier / endothelial": (
        "CLDN5", "TJP1", "PECAM1", "CDH5", "ICAM", "VCAM", "ATP1A2", "EDN1",
        "EDNRA", "NOS3", "KCNK3", "FOXF2", "MMP2", "MMP9", "MMP12", "ITGAV",
        "ITGB3", "ITGA2B", "AGER",
    ),
    "Coagulation / vascular": (
        "F2", "F3", "F5", "FGA", "PLAT", "PLAU", "PROC", "SERPINE1", "VKORC1",
        "MTHFR", "LPA", "APOA1", "ALB", "HBA1", "ABO", "ACE", "ACVRL1", "PDE3A",
        "SH2B3", "PLA2G7",
    ),
    "Neuronal / synaptic": (
        "CAMK2A", "CAMK2B", "GRIA1", "GRIN1", "SLC17A7", "MAP2", "ANK2", "L1CAM",
        "CACNA1C", "GABRA6", "HTR4", "ADORA1", "NPY", "RGS7", "FOS", "JUNB",
        "ATF3", "Neurogranin", "IGF1", "NGF", "FGF2",
    ),
    "Neurogenesis / cytoskeletal": (
        "DCX", "SPTAN1", "UCHL1",
    ),
    "Glial / myelin": (
        "GFAP", "AIF1", "MBP", "SLC1A3", "GLT-1", "S100B", "S100A5", "GPR17",
        "CYBB", "ALOX5", "ALOX5AP",
    ),
    "Injury biomarker": (
        "Neurofilament Light (NfL)", "Tau", "TDP43", "Amyloid-beta", "SAA1",
        "SERPINH1", "HSPD1",
    ),
    "Oxidative stress / metabolism": (
        "SOD1", "SOD2", "NOS2", "XDH", "UCP2", "SIRT1", "PRKAA2", "SMOX", "PAOX",
        "ALDH1", "ADH1B", "VDR", "FGF21", "SPHK2", "PLD1",
    ),
    "Transcription / chromatin / other": (
        "HDAC9", "KDM5D", "CASZ1", "NKX2-5", "PITX2", "TBX3", "ZFHX3", "ZCCHC14",
        "PRPF8", "CDK6", "APC", "LRCH1", "SH3PXD2A", "TSPAN2", "WNT2B",
    ),
}

# Ensembl Compara returns one-to-many orthology for a few panel members and its
# first hit is not always the canonical mouse gene. Where the intended ortholog
# is unambiguous it is pinned here; it must still exist on the assay to be used.
ORTHOLOG_PREFERENCE: dict[str, str] = {
    "ALDH1A1": "Aldh1a1",
    "HBA2": "Hba-a1",
    "SAA1": "Saa1",
}


def category_of(analyte: str) -> str:
    for cat, members in CATEGORIES.items():
        if analyte in members:
            return cat
    return "Unassigned"


def get_json(url: str, *, data: dict | None = None, retries: int = 3):
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "osdr-geomx-reanalysis/1.0",
    }
    body = json.dumps(data).encode() if data is not None else None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=body, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                return None
            time.sleep(1.5 * attempt)
        except (urllib.error.URLError, TimeoutError, OSError):
            time.sleep(1.5 * attempt)
    return None


def _cached_lookup(ensids: list[str], cache_name: str, offline: bool, label: str):
    path = CACHE / cache_name
    cached = json.loads(path.read_text()) if path.exists() else {}
    todo = [e for e in ensids if e and e not in cached]
    if todo and not offline:
        print(f"  Ensembl {label} for {len(todo)} id(s) ...")
        for i in range(0, len(todo), 100):
            chunk = todo[i : i + 100]
            res = get_json(f"{ENSEMBL}/lookup/id", data={"ids": chunk})
            if res:
                for k, v in res.items():
                    cached[k] = v or {}
            for k in chunk:
                cached.setdefault(k, {})
            time.sleep(0.3)
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cached, indent=1))
    return cached


def lookup_orthologs(ensids: list[str], offline: bool) -> dict[str, list[dict]]:
    path = CACHE / "mouse_orthologs.json"
    cached = json.loads(path.read_text()) if path.exists() else {}
    todo = [e for e in ensids if e and e not in cached]
    if todo and not offline:
        print(f"  Ensembl Compara orthology for {len(todo)} gene(s) ...")
        for n, ensid in enumerate(todo, 1):
            res = get_json(
                f"{ENSEMBL}/homology/id/human/{ensid}"
                "?target_species=mus_musculus&type=orthologues"
                "&sequence=none&format=condensed"
            )
            hits = []
            if res and res.get("data"):
                for h in res["data"][0].get("homologies", []):
                    hits.append({"id": h.get("id"), "type": h.get("type")})
            cached[ensid] = hits
            if n % 25 == 0:
                CACHE.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(cached, indent=1))
                print(f"    {n}/{len(todo)}")
            time.sleep(0.12)
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cached, indent=1))
    return cached


def read_panel() -> pd.DataFrame:
    """Read the lead-maintained panel CSV, tolerating trailing empty columns."""
    raw = pd.read_csv(PANEL_CSV, dtype=str)
    raw.columns = [str(c).strip() for c in raw.columns]
    cols = list(raw.columns)
    analyte_col = next(c for c in cols if c.lower().startswith("analyte"))
    ensid_col = next(c for c in cols if c.lower().startswith("ensid"))
    other = [c for c in cols if c not in (analyte_col, ensid_col)]
    df = pd.DataFrame(
        {
            "analyte": raw[analyte_col].astype(str).str.strip(),
            "ensid": raw[ensid_col].astype(str).str.strip(),
        }
    )
    # Any remaining column may carry a free-text note (e.g. "Iba1" beside AIF1).
    if other:
        notes = raw[other].bfill(axis=1).iloc[:, 0]
        df["panel_note"] = notes.fillna("").astype(str).str.strip()
    else:
        df["panel_note"] = ""
    df = df[(df["analyte"] != "") & (df["ensid"].str.startswith("ENSG"))]
    return df.reset_index(drop=True)


def resolve(panel: pd.DataFrame, offline: bool, targets: pd.Index) -> pd.DataFrame:
    human = _cached_lookup(panel["ensid"].tolist(), "human_lookup.json", offline, "lookup")
    orthos = lookup_orthologs(panel["ensid"].tolist(), offline)
    mouse_ids = sorted({h["id"] for v in orthos.values() for h in v if h.get("id")})
    mouse = _cached_lookup(mouse_ids, "mouse_lookup.json", offline, "mouse lookup")

    panel_lower = {t.lower(): t for t in targets}
    order = {"ortholog_one2one": 0, "ortholog_one2many": 1, "ortholog_many2many": 2}
    rows = []
    for r in panel.itertuples():
        h = human.get(r.ensid, {}) or {}
        h_sym = h.get("display_name", "")
        hits = sorted(orthos.get(r.ensid, []) or [],
                      key=lambda x: order.get(x.get("type", ""), 9))
        m_syms = []
        for hit in hits:
            info = mouse.get(hit.get("id") or "", {}) or {}
            if info.get("display_name"):
                m_syms.append((info["display_name"], hit.get("type", ""), hit.get("id")))
        pref = ORTHOLOG_PREFERENCE.get(h_sym)
        if pref:
            m_syms.sort(key=lambda x: 0 if x[0].lower() == pref.lower() else 1)

        matched, route, mtype, m_sym, m_id = "", "", "", "", ""
        for sym, typ, mid in m_syms:
            if sym.lower() in panel_lower:
                matched, route, mtype, m_sym, m_id = panel_lower[sym.lower()], "ensembl_ortholog", typ, sym, mid
                break
        if not matched and m_syms:
            m_sym, mtype, m_id = m_syms[0]
        if not matched:
            guess = h_sym.capitalize() if h_sym else ""
            if guess and guess.lower() in panel_lower:
                matched, route = panel_lower[guess.lower()], "symbol_titlecase_fallback"
                m_sym = m_sym or guess

        label = r.analyte.upper().replace("-", "").replace(" ", "")
        sym_norm = (h_sym or "").upper().replace("-", "")
        rows.append(
            {
                "analyte": r.analyte,
                "panel_note": r.panel_note,
                "ensid": r.ensid,
                "human_symbol_ensembl": h_sym,
                "human_biotype": h.get("biotype", ""),
                "label_matches_ensembl": bool(h_sym and label == sym_norm),
                "n_mouse_orthologs": len(m_syms),
                "mouse_ortholog_symbol": m_sym,
                "mouse_ortholog_id": m_id,
                "orthology_type": mtype,
                "geomx_target": matched,
                "on_geomx_panel": bool(matched),
                "mapping_route": route
                or ("no_mouse_ortholog" if not m_syms else "not_on_geomx_panel"),
                "category": category_of(r.analyte),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="use only the cached Ensembl responses")
    args = ap.parse_args()

    if not PANEL_CSV.exists():
        print(f"ERROR: panel not found at {PANEL_CSV}", file=sys.stderr)
        return 1
    panel = read_panel()
    print(f"Panel: {len(panel)} analytes, {panel['ensid'].nunique()} unique Ensembl IDs")
    print(f"  source: {PANEL_CSV.name}")

    targets = C.read_q3(C.STUDIES[0].code).drop(index=[NEG_PROBE], errors="ignore").index
    mapped = resolve(panel, args.offline, targets)

    C.CONFIG.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(C.CONFIG / "cns_targets_mapped.csv", index=False)
    print(f"  wrote {(C.CONFIG / 'cns_targets_mapped.csv').relative_to(C.ROOT)}")

    C.ensure_dirs()
    C.write_csv(mapped, "target_panel_mapping.csv")

    present = mapped.loc[mapped["on_geomx_panel"], "geomx_target"].drop_duplicates().tolist()
    (C.CONFIG / "focus_genes.txt").write_text(
        "\n".join(
            [
                "# CNS target panel, mapped to the mouse GeoMx target space.",
                "# Generated by scripts/map_targets.py from the root CNS mRNA targets.csv.",
                f"# {len(present)} of {len(mapped)} analytes are measured by this assay.",
                "",
            ]
            + sorted(present)
        )
        + "\n"
    )
    all_orth = sorted({s for s in mapped["mouse_ortholog_symbol"] if s})
    (C.CONFIG / "cns_targets_mouse_symbols.txt").write_text(
        "\n".join(
            [
                "# CNS target panel as mouse gene symbols, via Ensembl Compara.",
                "# Assay-independent: use this for any mouse dataset.",
                "# NOTE: this is a CNS injury / neuroinflammation / neurovascular panel.",
                "# It is not an EV panel; CD9/CD63/CD81 are not part of it.",
                f"# {len(all_orth)} symbols from {len(mapped)} analytes.",
                "",
            ]
            + all_orth
        )
        + "\n"
    )
    print(f"  wrote config/focus_genes.txt ({len(present)} genes)")
    print(f"  wrote config/cns_targets_mouse_symbols.txt ({len(all_orth)} symbols)")

    # EV / exosome reference markers, reported separately from the panel.
    ev_rows = []
    ev_human = _cached_lookup(list(EV_MARKERS.values()), "human_lookup.json",
                              args.offline, "lookup")
    ev_orth = lookup_orthologs(list(EV_MARKERS.values()), args.offline)
    ev_mouse_ids = [h["id"] for v in ev_orth.values() for h in v if h.get("id")]
    ev_mouse = _cached_lookup(ev_mouse_ids, "mouse_lookup.json", args.offline, "mouse lookup")
    tl = {t.lower(): t for t in targets}
    for sym, ensid in EV_MARKERS.items():
        hits = ev_orth.get(ensid, []) or []
        m = ""
        for h in hits:
            info = ev_mouse.get(h.get("id") or "", {}) or {}
            if info.get("display_name"):
                m = info["display_name"]
                break
        ev_rows.append(
            {
                "ev_marker": sym,
                "ensid": ensid,
                "human_symbol_ensembl": (ev_human.get(ensid, {}) or {}).get("display_name", ""),
                "mouse_ortholog_symbol": m,
                "in_cns_target_panel": bool((mapped["analyte"].str.upper() == sym).any()),
                "geomx_target": tl.get(m.lower(), "") if m else "",
                "measurable_on_geomx": bool(m and m.lower() in tl),
            }
        )
    ev = pd.DataFrame(ev_rows)
    C.write_csv(ev, "ev_marker_status.csv")

    cov = (
        mapped.groupby("category")
        .agg(n_analytes=("analyte", "size"), n_on_panel=("on_geomx_panel", "sum"))
        .reset_index()
    )
    cov["pct_on_panel"] = (cov["n_on_panel"] / cov["n_analytes"] * 100).round(1)
    C.write_csv(cov.sort_values("n_analytes", ascending=False),
                "target_panel_coverage_by_category.csv")

    pd.set_option("display.width", 220)
    print(f"\nMapped {len(present)}/{len(mapped)} analytes onto the GeoMx mouse panel")
    print("\nCoverage by functional category:")
    print(cov.sort_values("n_analytes", ascending=False).to_string(index=False))
    print("\nMapping routes:")
    print(mapped["mapping_route"].value_counts().to_string())

    print("\nEV / exosome reference markers (not part of the CNS panel):")
    print(ev.to_string(index=False))

    bad = mapped[~mapped["label_matches_ensembl"]]
    if len(bad):
        print(f"\nLabel differs from the Ensembl symbol ({len(bad)}) -- "
              "expected for protein/assay names, review for genuine mismatches:")
        print(bad[["analyte", "panel_note", "ensid", "human_symbol_ensembl",
                   "mouse_ortholog_symbol", "on_geomx_panel"]].to_string(index=False))

    lost = mapped[~mapped["on_geomx_panel"]]
    if len(lost):
        print(f"\nNot measured by this assay ({len(lost)}):")
        print(lost[["analyte", "human_symbol_ensembl", "mouse_ortholog_symbol",
                    "mapping_route"]].to_string(index=False))
    unassigned = mapped[mapped["category"] == "Unassigned"]
    if len(unassigned):
        print(f"\nWARNING: {len(unassigned)} analyte(s) have no functional category:")
        print(unassigned[["analyte", "human_symbol_ensembl"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
