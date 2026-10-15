"""Worker A, step 3: protocol-factor cross-tabulation for the 7 core and 4 external deposits.

Reads ONLY deposited metadata in the read-only source package and writes
A_design_matrix.csv (one row per distinct combination of design and protocol
factors, with library and animal counts) plus A_contrast_aliasing.csv (for every
contrast the manuscript uses: which protocol factors are matched, aliased or
partially confounded between the two compared groups).

Group labels come from the deposited Source Name (as in the paper). Protocol
factors come from the deposited ISA fields; nothing is inferred except where a
column says 'protocol text' (taken verbatim from the deposited study protocols).
"""
from pathlib import Path
import json
import pandas as pd

PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = Path(__file__).resolve().parent
EXT = PKG / "extensions/strengthening_20261015/external_data"


def first(row, names):
    for n in names:
        if n in row and pd.notna(row[n]) and str(row[n]).strip():
            return str(row[n]).strip()
    return "not deposited"


def core_rows():
    spec = {"561": ("RRRM-2", "cerebellum (right)"), "562": ("RRRM-2", "hippocampus (right)"),
            "613": ("RRRM-2", "left cerebral hemisphere"), "563": ("RR-10", "cerebellum (right)"),
            "564": ("RR-10", "hippocampus (right)"), "612": ("RR-10", "left cerebral hemisphere"),
            "525": ("RR-23", "cerebellum")}
    rows = []
    for n, (mission, region) in spec.items():
        m = pd.read_csv(PKG / f"inputs/OSD-{n}_metadata.csv.gz", dtype=str)
        for _, r in m.iterrows():
            sn = r["Source Name"]
            if mission == "RRRM-2":
                flight = "FLT" if "_FL_" in sn else "GC"
                cohort = "ISS-T" if "_ISS_" in sn else "LAR"
                age = first(r, ["Factor Value[Age]"]) + " wk (launch); " + first(
                    r, ["Parameter Value[Age at euthanasia]", "Characteristics[Age at Euthanasia]"]) + " wk at euthanasia"
                design_loc = "ISS (per Source Name)" if (flight == "FLT" and cohort == "ISS-T") else "Earth (per Source Name)"
            elif mission == "RR-23":
                tok = sn.split("_")[-1]
                flight = {"F": "FLT", "G": "HGC", "V": "VGC"}[tok[0]]
                cohort = "LAR (all groups returned/handled on Earth)"
                age = first(r, ["Characteristics[Age at Launch]"]) + " wk at launch"
                design_loc = "Earth (FLT returned live)"
            else:
                flight = "FLT" if "_FL-" in sn else "GC"
                cohort = "carcass (euthanized at end of mission, frozen)"
                age = first(r, ["Characteristics[Age at Launch]"]) + " wk launch; " + first(
                    r, ["Parameter Value[Age at Euthanasia]", "Characteristics[Age at Euthanasia]"]) + " wk euthanasia"
                design_loc = ("ISS (protocol text: flight mice euthanized after 28-29 d in microgravity)"
                              if flight == "FLT" else "Earth (ISSES ground control; same schedule)")
            rows.append(dict(
                deposit=f"OSD-{n}", mission=mission, region=region,
                exposure=flight, collection_cohort=cohort,
                euthanasia_method=first(r, ["Parameter Value[Euthanasia Method]"]),
                deposited_euthanasia_location=first(r, ["Factor Value[Euthanasia Location]"]),
                design_euthanasia_location=design_loc,
                sample_preservation=first(r, ["Parameter Value[Sample Preservation Method]"]),
                carcass_preservation=first(r, ["Parameter Value[Carcass Preservation Method]"]),
                dissection=first(r, ["Comment[Source Description]"]),
                euthanasia_date=first(r, ["Comment[Euthanasia Date]", "Comment[Euthanasia Date (dd-mmm-yyyy)]"]),
                dissection_date=first(r, ["Comment[BSP Dissection Date]", "Comment[BSP Dissection Date (dd-mmm-yyyy)]"]),
                recovery_on_earth_days=first(r, ["Parameter Value[Recovery Duration on Earth]"]),
                habitat=first(r, ["Parameter Value[habitat]"]),
                sex=first(r, ["Characteristics[Sex]"]), strain=first(r, ["Characteristics[Strain]"]),
                age=age, rin="not deposited", unit_type="animal (technical libraries summed)",
                animal=sn, library=r["sample"]))
    return rows


def external_rows():
    prot = {"202": ("Earth: Co-57 gamma x hindlimb unloading", "not deposited in sample or protocol fields"),
            "352": ("RR-3 flight vs ground", "ketamine/xylazine/acepromazine (Parameter Value[Euthanasia Chemical])"),
            "374": ("Earth: 1 Gy protons vs sham", None), "611": ("Earth: Si-28 0.3/0.6 Gy vs sham", None)}
    rows = []
    for n, (ctx, eu) in prot.items():
        m = pd.read_csv(EXT / f"OSD-{n}/metadata.csv", dtype=str)
        for _, r in m.iterrows():
            if n == "202":
                exposure = f"{r['radiation']} | {r['hindlimb_unloading']}"
                cohort = f"collected {r['collection_after_treatment']} month after treatment"
                dissect, loc = "not deposited", "Earth"
            elif n == "352":
                exposure = "FLT" if r["spaceflight"] == "Space Flight" else "GC"
                cohort = f"{r['duration']} d mission"
                dissect = "previously dissected carcass, frozen carcass (BSP dissection on Earth; protocol text)"
                loc = "not deposited"
            else:
                exposure = f"{r['radiation']} {r['radiation_dose_gy']} Gy"
                cohort = f"collected {r['collection_after_treatment']} week after exposure"
                dissect, loc = "not deposited", "Earth"
            rows.append(dict(
                deposit=f"OSD-{n}", mission=r["mission"], region=r["tissue"], exposure=exposure,
                collection_cohort=cohort,
                euthanasia_method=eu if eu else (r["euthanasia"] if pd.notna(r["euthanasia"]) else "not deposited"),
                deposited_euthanasia_location="not deposited", design_euthanasia_location=loc,
                sample_preservation=r["preservation"], carcass_preservation="not deposited",
                dissection=dissect, euthanasia_date=r["sacrifice_date"] if pd.notna(r["sacrifice_date"]) else "not deposited",
                dissection_date="not deposited", recovery_on_earth_days="not deposited",
                habitat=r["habitat"] if pd.notna(r["habitat"]) else "not deposited",
                sex=r["sex"], strain=r["strain"], age=f"{r['age']} {r['age_unit']}",
                rin=r["rin"] if pd.notna(r["rin"]) else "not deposited",
                unit_type=r["replicate_unit"] + ("" if r["primary_analysis_included"] == "True" else " [NOT in SampleTable; sensitivity only]"),
                animal=r["independent_unit_id"], library=r["sample_id"]))
    return rows


FACTORS = ["deposit", "mission", "region", "exposure", "collection_cohort", "euthanasia_method",
           "deposited_euthanasia_location", "design_euthanasia_location", "sample_preservation", "preservation_class",
           "carcass_preservation", "dissection", "euthanasia_date", "dissection_date",
           "recovery_on_earth_days", "habitat", "sex", "strain", "age", "unit_type"]


def bulk_preservation(v):
    """The LAR/RR-10 deposits use one string for two assays; the bulk RNA-seq step is the LN2 part."""
    if isinstance(v, str) and "for bulkRNAseq" in v and "liquid nitrogen" in v.lower():
        return "Liquid Nitrogen (bulk RNA-seq; spatial assay used OCT on dry ice)"
    return v


def main():
    lib = pd.DataFrame(core_rows() + external_rows())
    lib["bulk_RNAseq_preservation"] = lib.sample_preservation.map(bulk_preservation)
    lib["preservation_class"] = lib.bulk_RNAseq_preservation.str.split(" ").str[0:2].str.join(" ")
    lib.to_csv(OUT / "A_design_library_level.csv", index=False)
    # Collapse euthanasia dates to the set per cell so cells stay readable.
    g = lib.copy()
    tab = (g.groupby([f for f in FACTORS if f != "euthanasia_date"], dropna=False)
             .agg(n_libraries=("library", "size"), n_units=("animal", "nunique"),
                  euthanasia_dates=("euthanasia_date", lambda s: "; ".join(sorted(set(s)))),
                  rin_values=("rin", lambda s: "; ".join(sorted(set(s)))))
             .reset_index())
    tab.to_csv(OUT / "A_design_matrix.csv", index=False)
    print(tab[["deposit", "exposure", "collection_cohort", "euthanasia_method", "deposited_euthanasia_location",
               "sample_preservation", "n_libraries", "n_units", "euthanasia_dates"]].to_string())

    # Aliasing per manuscript contrast: compare factor value sets between groups.
    def sel(dep, **kw):
        s = lib[lib.deposit.isin(dep if isinstance(dep, list) else [dep])]
        for k, v in kw.items():
            s = s[s[k].isin(v if isinstance(v, list) else [v])]
        return s
    rr2 = ["OSD-561", "OSD-562"]
    contrasts = [
        ("GC ISS-T minus GC LAR (ground 'collection' effect), Cb+Hip", sel(rr2, exposure="GC", collection_cohort="ISS-T"), sel(rr2, exposure="GC", collection_cohort="LAR")),
        ("FLT ISS-T minus FLT LAR, Cb+Hip", sel(rr2, exposure="FLT", collection_cohort="ISS-T"), sel(rr2, exposure="FLT", collection_cohort="LAR")),
        ("FLT minus GC within ISS-T, Cb+Hip", sel(rr2, exposure="FLT", collection_cohort="ISS-T"), sel(rr2, exposure="GC", collection_cohort="ISS-T")),
        ("FLT minus GC within LAR, Cb+Hip", sel(rr2, exposure="FLT", collection_cohort="LAR"), sel(rr2, exposure="GC", collection_cohort="LAR")),
        ("OSD-613 FLT minus GC (ISS-T, left hemisphere)", sel("OSD-613", exposure="FLT"), sel("OSD-613", exposure="GC")),
        ("RR-10 FLT minus GC (OSD-563 cerebellum, OSD-564 hippocampus)", sel(["OSD-563", "OSD-564"], exposure="FLT"), sel(["OSD-563", "OSD-564"], exposure="GC")),
        ("RR-10 FLT minus GC (OSD-612 left hemisphere)", sel("OSD-612", exposure="FLT"), sel("OSD-612", exposure="GC")),
        ("RR-23 HGC minus VGC (OSD-525 ground anchor)", sel("OSD-525", exposure="HGC"), sel("OSD-525", exposure="VGC")),
        ("RR-23 FLT minus HGC (OSD-525)", sel("OSD-525", exposure="FLT"), sel("OSD-525", exposure="HGC")),
        ("RR-23 FLT minus VGC (OSD-525)", sel("OSD-525", exposure="FLT"), sel("OSD-525", exposure="VGC")),
        ("RR-3 FLT minus GC (OSD-352)", sel("OSD-352", exposure="FLT"), sel("OSD-352", exposure="GC")),
    ]
    check = ["euthanasia_method", "deposited_euthanasia_location", "design_euthanasia_location",
             "preservation_class", "carcass_preservation", "dissection", "euthanasia_date", "dissection_date",
             "recovery_on_earth_days", "habitat", "sex", "strain", "age", "collection_cohort"]
    out = []
    for name, a, b in contrasts:
        for f in check:
            va, vb = set(a[f]), set(b[f])
            if va == vb and len(va) == 1:
                status = "matched" if list(va)[0] != "not deposited" else "not deposited (cannot assess)"
            elif va == vb:
                status = "same levels in both groups (balanced/blocked or adjusted)"
            elif not (va & vb):
                status = "ALIASED (disjoint values)"
            else:
                status = "partially overlapping (unbalanced)"
            out.append(dict(contrast=name, factor=f, status=status, group1_values=" | ".join(sorted(va)),
                            group2_values=" | ".join(sorted(vb)), n_units_group1=a.animal.nunique(),
                            n_units_group2=b.animal.nunique()))
    al = pd.DataFrame(out)
    al.to_csv(OUT / "A_contrast_aliasing.csv", index=False)
    print(al[al.status.str.startswith(("ALIASED", "partially"))][["contrast", "factor", "status"]].to_string())


def location_conflicts():
    """Same animal, two regional deposits: compare deposited euthanasia location, RFID and date."""
    a = pd.read_csv(PKG / "inputs/OSD-561_metadata.csv.gz", dtype=str).set_index("Source Name")
    b = pd.read_csv(PKG / "inputs/OSD-562_metadata.csv.gz", dtype=str).set_index("Source Name")
    L = "Factor Value[Euthanasia Location]"
    rows = []
    for sn in sorted(set(a.index) | set(b.index)):
        if "_ISS_" not in sn:
            continue
        exp = "On ISS" if "_FL_" in sn else "On Earth"
        ca = a.loc[sn, L] if sn in a.index else None
        cb = b.loc[sn, L] if sn in b.index else None
        if (ca is not None and ca != exp) or (cb is not None and cb != exp):
            rows.append(dict(animal=sn, expected_from_Source_Name=exp, OSD561_cerebellum=ca, OSD562_hippocampus=cb,
                             same_RFID=(a.loc[sn, "Comment[RFID]"].upper() == b.loc[sn, "Comment[RFID]"].upper()) if (ca is not None and cb is not None) else None,
                             same_euthanasia_date=(a.loc[sn, "Comment[Euthanasia Date]"] == b.loc[sn, "Comment[Euthanasia Date]"]) if (ca is not None and cb is not None) else None,
                             age=(a.loc[sn, "Factor Value[Age]"] if ca is not None else b.loc[sn, "Factor Value[Age]"]),
                             verdict=("deposits contradict each other for the same animal" if (ca is not None and cb is not None and ca != cb)
                                      else "both deposits disagree with Source Name")))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "A_location_conflicts.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
    location_conflicts()
