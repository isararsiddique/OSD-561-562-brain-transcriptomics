#!/usr/bin/env python3
"""Map the CNS / EV target panel (human) onto the mouse GeoMx target space.

The supplied panel lists human analytes with human Ensembl gene IDs, several of
them protein or assay names rather than gene symbols ("Amyloid-beta", "GLT-1",
"Iba1", "SBDP (SNTF)", "Neurofilament Light (NfL)").  The data are mouse.  Two
things therefore have to happen before any analysis, and neither may be guessed:

  1. Resolve each Ensembl gene ID to its current human symbol.  The ENSG is the
     stable identifier, so it -- not the spreadsheet label -- defines the gene.
     Where the two disagree the discrepancy is recorded rather than silently
     resolved.
  2. Map each human gene to its mouse ortholog through Ensembl Compara, keeping
     the orthology type.  Title-casing a human symbol is *not* a safe shortcut:
     human CXCL8 has no mouse ortholog at all, and several panel members were
     renamed between the two genomes.

Everything is cached under data/mapping/ so the network is hit once.

Outputs
  config/cns_ev_targets_mapped.csv    full audit trail, one row per analyte
  config/focus_genes.txt              mouse symbols present on the GeoMx panel
  results/tables/csv/target_panel_mapping.csv
  results/tables/csv/target_panel_coverage_by_category.csv

Usage:  python scripts/map_targets.py [--offline]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

import common as C

PANEL_CSV = (
    C.ROOT.parent
    / "CNS mRNA targets & Others in EV.xlsx - organized targets - "
      "CNS mRNA targets & Others in EV.xlsx - organized targets.csv"
)
CACHE = C.DATA / "mapping"
ENSEMBL = "https://rest.ensembl.org"
NEG_PROBE = "NegProbe-WTX"

# Functional grouping of the panel. Edit freely -- it only affects how results
# are aggregated and plotted, never which genes are tested.
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
    "Glial / myelin": (
        "GFAP", "Iba1", "MBP", "SLC1A3", "GLT-1", "S100B", "S100A5", "GPR17",
        "CYBB", "ALOX5", "ALOX5AP",
    ),
    "Injury biomarker": (
        "Neurofilament Light (NfL)", "Tau", "UCHL1", "TDP43", "Amyloid-beta",
        "SBDP (SNTF)", "Serum amyloid alpha (SAA)", "SERPINH1", "HSPD1",
    ),
    "Oxidative stress / metabolism": (
        "SOD1", "SOD2", "NOS2", "XDH", "UCP2", "SIRT1", "PRKAA2", "SMOX", "PAOX",
        "ALDH1", "ADH1B", "VDR", "FGF21", "SPHK2", "PLD1",
    ),
    "Transcription / chromatin / other": (
        "HDAC9", "KDM5D", "CASZ1", "NKX2-5", "PITX2", "TBX3", "ZFHX3", "ZCCHC14",
        "PRPF8", "CDK6", "APC-CC1", "LRCH1", "SH3PXD2A", "TSPAN2", "WNT2B",
    ),
}


def category_of(analyte: str) -> str:
    for cat, members in CATEGORIES.items():
        if analyte in members:
            return cat
    return "Unassigned"


# --------------------------------------------------------------------------- #
# Curated corrections to the supplied spreadsheet
# --------------------------------------------------------------------------- #
# Resolving every ENSID against Ensembl showed that a handful point to genes
# unrelated to their label. Left uncorrected, the analysis would report results
# for the wrong gene under a recognised biomarker name -- so each is corrected
# here with the intended gene, the expected symbol it must resolve to, and the
# reason. The expected symbol is verified against Ensembl at run time; any
# disagreement is reported rather than trusted.
#
#   analyte -> (corrected ENSID, expected human symbol, reason)
CORRECTIONS: dict[str, tuple[str, str, str]] = {
    "Iba1": (
        "ENSG00000204472", "AIF1",
        "supplied ENSG00000153406 is NMRAL1; Iba1 is the protein name for AIF1",
    ),
    "Neurogranin": (
        "ENSG00000154146", "NRGN",
        "supplied ENSG00000101191 is DIDO1; neurogranin is NRGN",
    ),
    "SBDP (SNTF)": (
        "ENSG00000197694", "SPTAN1",
        "supplied ENSG00000077279 is DCX; SBDP/SNTF are alpha-II spectrin "
        "breakdown products, i.e. SPTAN1",
    ),
    "Serum amyloid alpha (SAA)": (
        "ENSG00000173432", "SAA1",
        "supplied ENSG00000154803 is FLCN; serum amyloid A is SAA1",
    ),
    "APC-CC1": (
        "ENSG00000134982", "APC",
        "supplied ENSG00000135982 returns no gene; APC is ENSG00000134982 "
        "(transposed digits)",
    ),
    "PYCARD": (
        "ENSG00000103490", "PYCARD",
        "supplied ENSG00000103483 returns no gene; PYCARD (ASC) is ENSG00000103490",
    ),
}

# Ensembl Compara returns one-to-many orthology for a few panel members and the
# first hit is not the canonical mouse gene. Where the intended mouse ortholog
# is unambiguous it is pinned here, and it still has to exist on the GeoMx panel
# to be used.
ORTHOLOG_PREFERENCE: dict[str, str] = {
    "ALDH1A1": "Aldh1a1",
    "HBA2": "Hba-a1",
    "SAA1": "Saa1",
}


def get_json(url: str, *, data: dict | None = None, retries: int = 3):
    headers = {"Content-Type": "application/json", "Accept": "application/json",
               "User-Agent": "osdr-geomx-reanalysis/1.0"}
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


def lookup_symbols(ensids: list[str], offline: bool) -> dict[str, dict]:
    path = CACHE / "human_lookup.json"
    cached = json.loads(path.read_text()) if path.exists() else {}
    todo = [e for e in ensids if e not in cached]
    if todo and not offline:
        print(f"  Ensembl lookup for {len(todo)} gene IDs ...")
        for i in range(0, len(todo), 100):
            chunk = todo[i : i + 100]
            res = get_json(f"{ENSEMBL}/lookup/id", data={"ids": chunk})
            if res is None:
                print("    lookup batch failed")
                continue
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
    todo = [e for e in ensids if e not in cached]
    if todo and not offline:
        print(f"  Ensembl Compara orthology for {len(todo)} genes ...")
        for n, ensid in enumerate(todo, 1):
            url = (
                f"{ENSEMBL}/homology/id/human/{ensid}"
                "?target_species=mus_musculus&type=orthologues"
                "&sequence=none&format=condensed"
            )
            res = get_json(url)
            hits = []
            if res and res.get("data"):
                for h in res["data"][0].get("homologies", []):
                    hits.append({"id": h.get("id"), "type": h.get("type")})
            cached[ensid] = hits
            if n % 25 == 0:
                print(f"    {n}/{len(todo)}")
                CACHE.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(cached, indent=1))
            time.sleep(0.12)
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cached, indent=1))
    return cached


def mouse_symbols(mouse_ids: list[str], offline: bool) -> dict[str, dict]:
    path = CACHE / "mouse_lookup.json"
    cached = json.loads(path.read_text()) if path.exists() else {}
    todo = [m for m in mouse_ids if m and m not in cached]
    if todo and not offline:
        print(f"  Ensembl lookup for {len(todo)} mouse gene IDs ...")
        for i in range(0, len(todo), 100):
            chunk = todo[i : i + 100]
            res = get_json(f"{ENSEMBL}/lookup/id", data={"ids": chunk})
            if res is None:
                continue
            for k, v in res.items():
                cached[k] = v or {}
            for k in chunk:
                cached.setdefault(k, {})
            time.sleep(0.3)
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cached, indent=1))
    return cached


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="use only the cached Ensembl responses")
    args = ap.parse_args()

    if not PANEL_CSV.exists():
        print(f"ERROR: target panel not found at {PANEL_CSV}", file=sys.stderr)
        return 1
    panel = pd.read_csv(PANEL_CSV)
    panel.columns = [c.strip() for c in panel.columns]
    panel = panel.rename(columns={"Analyte": "analyte", "ENSID": "ensid"})
    panel["analyte"] = panel["analyte"].astype(str).str.strip()
    panel["ensid"] = panel["ensid"].astype(str).str.strip()
    print(f"Panel: {len(panel)} analytes, {panel['ensid'].nunique()} unique Ensembl IDs")

    panel["ensid_supplied"] = panel["ensid"]
    panel["correction_reason"] = ""
    panel["expected_symbol"] = ""
    for i, r in panel.iterrows():
        fix = CORRECTIONS.get(r["analyte"])
        if fix:
            panel.at[i, "ensid"] = fix[0]
            panel.at[i, "expected_symbol"] = fix[1]
            panel.at[i, "correction_reason"] = fix[2]
    n_fixed = int((panel["ensid"] != panel["ensid_supplied"]).sum())
    if n_fixed:
        print(f"  applied {n_fixed} curated ENSID corrections (see CORRECTIONS in this script)")

    ensids = panel["ensid"].tolist()
    human = lookup_symbols(ensids, args.offline)
    orthos = lookup_orthologs(ensids, args.offline)
    all_mouse_ids = sorted({h["id"] for v in orthos.values() for h in v if h.get("id")})
    mouse = mouse_symbols(all_mouse_ids, args.offline)

    targets = C.read_q3(C.STUDIES[0].code).drop(index=[NEG_PROBE], errors="ignore").index
    panel_lower = {t.lower(): t for t in targets}

    rows = []
    for r in panel.itertuples():
        h = human.get(r.ensid, {}) or {}
        h_sym = h.get("display_name", "")
        hits = orthos.get(r.ensid, []) or []
        # Prefer one-to-one orthology, then one-to-many, then anything.
        order = {"ortholog_one2one": 0, "ortholog_one2many": 1, "ortholog_many2many": 2}
        hits = sorted(hits, key=lambda x: order.get(x.get("type", ""), 9))
        m_syms = []
        for hit in hits:
            info = mouse.get(hit.get("id") or "", {}) or {}
            sym = info.get("display_name")
            if sym:
                m_syms.append((sym, hit.get("type", ""), hit.get("id")))
        # Pin the canonical mouse gene where Compara's first hit is a paralog.
        pref = ORTHOLOG_PREFERENCE.get(h_sym)
        if pref:
            m_syms.sort(key=lambda x: 0 if x[0].lower() == pref.lower() else 1)

        matched, route, mtype, m_sym, m_id = "", "", "", "", ""
        for sym, typ, mid in m_syms:
            if sym.lower() in panel_lower:
                matched = panel_lower[sym.lower()]
                route = "ensembl_ortholog"
                mtype, m_sym, m_id = typ, sym, mid
                break
        if not matched and m_syms:
            m_sym, mtype, m_id = m_syms[0]
        if not matched:
            # Last resort: title-cased human symbol, flagged as such so it can
            # be audited or dropped.
            guess = h_sym.capitalize() if h_sym else ""
            if guess and guess.lower() in panel_lower:
                matched = panel_lower[guess.lower()]
                route = "symbol_titlecase_fallback"
                m_sym = m_sym or guess

        rows.append(
            {
                "analyte": r.analyte,
                "ensid": r.ensid,
                "ensid_supplied": r.ensid_supplied,
                "ensid_corrected": bool(r.ensid != r.ensid_supplied),
                "correction_reason": r.correction_reason,
                "correction_verified": (
                    bool(h_sym and r.expected_symbol and h_sym.upper() == r.expected_symbol.upper())
                    if r.ensid != r.ensid_supplied else None
                ),
                "human_symbol_ensembl": h_sym,
                "human_biotype": h.get("biotype", ""),
                "label_matches_ensembl": bool(
                    h_sym and r.analyte.upper().replace("-", "") == h_sym.upper().replace("-", "")
                ),
                "n_mouse_orthologs": len(m_syms),
                "mouse_ortholog_symbol": m_sym,
                "mouse_ortholog_id": m_id,
                "orthology_type": mtype,
                "geomx_target": matched,
                "on_geomx_panel": bool(matched),
                "mapping_route": route or ("no_mouse_ortholog" if not m_syms else "not_on_geomx_panel"),
                "category": category_of(r.analyte),
            }
        )

    mapped = pd.DataFrame(rows)
    C.CONFIG.mkdir(parents=True, exist_ok=True)
    out = C.CONFIG / "cns_ev_targets_mapped.csv"
    mapped.to_csv(out, index=False)
    print(f"  wrote {out.relative_to(C.ROOT)}")

    C.ensure_dirs()
    C.write_csv(mapped, "target_panel_mapping.csv")

    present = mapped.loc[mapped["on_geomx_panel"], "geomx_target"].drop_duplicates().tolist()
    focus = C.CONFIG / "focus_genes.txt"
    header = [
        "# CNS / EV mRNA target panel, mapped to mouse GeoMx targets.",
        "# Generated by scripts/map_targets.py -- edit cns_ev_targets_mapped.csv",
        "# or the source spreadsheet rather than this file.",
        f"# {len(present)} of {len(mapped)} analytes are measured by this assay.",
        "",
    ]
    focus.write_text("\n".join(header + sorted(present)) + "\n")
    print(f"  wrote {focus.relative_to(C.ROOT)} ({len(present)} genes)")

    # Full mouse ortholog list, independent of GeoMx panel membership, so the
    # sibling bulk RNA-seq project (which measures ~26k genes) can use the same
    # panel without redoing the mapping.
    all_orthologs = sorted(
        {s for s in mapped["mouse_ortholog_symbol"].dropna().tolist() if s}
    )
    ortho_path = C.CONFIG / "cns_ev_targets_mouse_symbols.txt"
    ortho_path.write_text(
        "\n".join(
            [
                "# CNS / EV target panel as mouse gene symbols.",
                "# Generated by scripts/map_targets.py from the project target",
                "# spreadsheet via Ensembl Compara orthology, including curated",
                "# ENSID corrections. Not filtered by assay content, so this is the",
                "# list to use for any mouse dataset.",
                f"# {len(all_orthologs)} symbols from {len(mapped)} analytes.",
                "",
            ]
            + all_orthologs
        )
        + "\n"
    )
    print(f"  wrote {ortho_path.relative_to(C.ROOT)} ({len(all_orthologs)} symbols)")

    cov = (
        mapped.groupby("category")
        .agg(n_analytes=("analyte", "size"), n_on_panel=("on_geomx_panel", "sum"))
        .reset_index()
    )
    cov["pct_on_panel"] = (cov["n_on_panel"] / cov["n_analytes"] * 100).round(1)
    cov = cov.sort_values("n_analytes", ascending=False)
    C.write_csv(cov, "target_panel_coverage_by_category.csv")

    pd.set_option("display.width", 220)
    print(f"\nMapped {len(present)}/{len(mapped)} analytes onto the GeoMx mouse panel")
    print("\nCoverage by functional category:")
    print(cov.to_string(index=False))
    print("\nMapping routes:")
    print(mapped["mapping_route"].value_counts().to_string())

    fixed = mapped[mapped["ensid_corrected"]]
    if len(fixed):
        print("\nCurated ENSID corrections (verified against Ensembl):")
        print(
            fixed[["analyte", "ensid_supplied", "ensid", "human_symbol_ensembl",
                   "correction_verified", "mouse_ortholog_symbol", "on_geomx_panel"]]
            .to_string(index=False)
        )
        unver = fixed[fixed["correction_verified"] == False]  # noqa: E712
        if len(unver):
            print("  WARNING: these corrections did not resolve to the expected symbol:")
            print(unver[["analyte", "ensid", "human_symbol_ensembl"]].to_string(index=False))

    bad = mapped[~mapped["label_matches_ensembl"] & ~mapped["ensid_corrected"]]
    if len(bad):
        print(f"\nLabel differs from the Ensembl symbol, left as supplied ({len(bad)}):")
        print(bad[["analyte", "ensid", "human_symbol_ensembl", "mouse_ortholog_symbol",
                   "on_geomx_panel"]].to_string(index=False))
    lost = mapped[~mapped["on_geomx_panel"]]
    if len(lost):
        print(f"\nNot measured by this assay ({len(lost)}):")
        print(lost[["analyte", "human_symbol_ensembl", "mouse_ortholog_symbol",
                    "mapping_route"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
