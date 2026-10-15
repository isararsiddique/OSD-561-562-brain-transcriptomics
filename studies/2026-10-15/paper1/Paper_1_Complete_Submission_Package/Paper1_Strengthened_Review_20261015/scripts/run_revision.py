#!/usr/bin/env python3
"""Reanalyse seven bulk brain deposits using only bundled processed inputs.

All transformations are sensitivity analyses. Removed contrast energy is not
an estimate of the fraction of false positives or true biology recovered.
"""
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from statistics_core import (Contrast, animal_permutations, bh, lambda_diagnostic,
                             learn_anchor, normalize, permutation_p, transform)

HERE = Path(__file__).resolve().parents[1]
DATASETS = {"525": ("RR-23", "cerebellum"), "561": ("RRRM-2", "cerebellum"),
            "562": ("RRRM-2", "hippocampus"), "613": ("RRRM-2", "left cerebral hemisphere"),
            "563": ("RR-10", "cerebellum"), "564": ("RR-10", "hippocampus"),
            "612": ("RR-10", "left cerebral hemisphere")}
PANELS = {
    "rPRG": "Fos Fosb Junb Egr1 Egr2 Egr4 Nr4a1 Nr4a3 Dusp1 Dusp5 Dusp6 Arc Npas4 Btg2 Ier2 Gadd45b Gadd45g Klf4 Klf2 Zfp36".split(),
    "dPRG": "Jun Jund Fosl1 Fosl2 Nr4a2 Egr3 Atf3 Csrnp1 Per1 Sik1 Rgs2 Trib1 Maff Myc".split(),
    "SRG": "Bdnf Nptx2 Homer1 Scg2 Inhba Vgf Crh Pcsk1 Serpinb2 Plat Ptgs2 Sst Npy Igf1".split(),
    "HSP": "Hspa1a Hspa1b Hspa5 Hspb1 Hsph1 Dnajb1 Dnaja1 Hsp90aa1 Bag3 Hspa8 Hspe1".split(),
    "Glia_Complement": "Gfap Aqp4 Aldh1l1 S100b Gja1 Cx3cr1 P2ry12 Tmem119 Aif1 C1qa C1qb C1qc Trem2 Tyrobp Csf1r Itgam Cd68 Ctss".split(),
    "BBB_Endothelium": "Cldn5 Ocln Tjp1 Pecam1 Cdh5 Flt1 Kdr Tek Tie1 Eng Vwf Slc2a1 Abcb1a Abcg2 Mfsd2a Prom1 Podxl Cav1".split(),
    "Circadian": "Arntl Clock Npas2 Per1 Per2 Per3 Cry1 Cry2 Nr1d1 Nr1d2 Dbp Tef Hlf".split(),
}
PANELS["Sensor31"] = PANELS["rPRG"] + PANELS["HSP"]


def check_input_integrity():
    records = json.loads((HERE / "inputs/provenance.json").read_text())
    for record in records:
        path = HERE / "inputs" / record["bundled_path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["bundled_sha256"]:
            raise ValueError(f"Input checksum mismatch: {path}")
    upstream = HERE / "inputs/upstream_count_provenance.json"
    if upstream.exists():
        for record in json.loads(upstream.read_text()):
            path = HERE / "inputs" / record["bundled_upstream_path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != record["bundled_upstream_sha256"]:
                raise ValueError(f"Upstream input checksum mismatch: {path}")


def field(row, names, default=""):
    for name in names:
        value = row.get(name, "")
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return default


def load_dataset(n, estimator="STAR"):
    inp = HERE / "inputs"
    suffix = "RSEM_" if estimator == "RSEM" else ""
    counts = pd.read_csv(inp / f"OSD-{n}_{suffix}counts.csv.gz", index_col=0)
    md = pd.read_csv(inp / f"OSD-{n}_metadata.csv.gz", dtype=str).set_index("sample")
    if not md.index.is_unique:
        raise ValueError(f"Duplicate metadata sample IDs in OSD-{n}")
    absent = set(counts.columns) - set(md.index)
    if absent:
        raise ValueError(f"Counts missing metadata: OSD-{n}: {absent}")
    md = md.loc[counts.columns]
    mission, tissue = DATASETS[n]
    rows = []
    for sample, row in md.iterrows():
        animal = field(row, ["Source Name"])
        if not animal:
            raise ValueError(f"Missing animal ID for {sample}")
        # Use deposited Source Name to identify shared animals and design cohort.
        # This does not resolve disagreements with the deposited location field.
        if mission == "RRRM-2":
            flight = "FLT" if "_FL_" in animal else "GC"
            cohort = "ISS_T" if "_ISS_" in animal else "LAR"
            age = field(row, ["Factor Value[Age]", "age"], "unknown")
        elif mission == "RR-23":
            flight = {"F": "FLT", "G": "HGC", "V": "VGC"}[animal.split("_")[-1][0]]
            cohort, age = flight, "single age"
        else:
            flight = "FLT" if "_FL-" in animal else "GC"
            cohort, age = "matched collection", "single age"
        location = field(row, ["Factor Value[Euthanasia Location]"])
        expected = "On ISS" if flight == "FLT" and cohort == "ISS_T" else "On Earth"
        conflict = bool(location and mission == "RRRM-2" and location != expected)
        rows.append({"accession": f"OSD-{n}", "mission": mission, "tissue": tissue,
                     "sample": sample, "animal": animal, "flight": flight, "cohort": cohort,
                     "age": age, "deposited_location": location, "design_location": expected,
                     "location_conflict": conflict, "group_basis": "Source Name / cohort identifier",
                     "strain": field(row, ["Characteristics[Strain]"]),
                     "sex": field(row, ["Characteristics[Sex]"]),
                     "euthanasia_date": field(row, ["Comment[Euthanasia Date]",
                                                   "Comment[Euthanasia Date (dd-mmm-yyyy)]"]),
                     "euthanasia_method": field(row, ["Parameter Value[Euthanasia Method]"]),
                     "deposited_material": field(row, ["Characteristics[Material Type]"])})
    manifest = pd.DataFrame(rows).set_index("sample")
    # Count independent libraries for each animal but sum them into one matrix column.
    # Apply the SAME collapse to STAR and RSEM, including every technical replicate.
    for col in ["flight", "cohort", "age", "tissue"]:
        if (manifest.groupby("animal")[col].nunique() > 1).any():
            raise ValueError(f"Conflicting {col} among technical replicates in OSD-{n}")
    c = counts.T.groupby(manifest.animal, sort=True).sum().T
    meta = manifest.reset_index().drop_duplicates("animal").set_index("animal", drop=False).loc[c.columns]
    meta["n_libraries"] = manifest.groupby("animal").size().reindex(meta.index)
    c.columns = [f"OSD-{n}:{a}" for a in c.columns]
    meta.index = c.columns
    y, sf = normalize(c)
    return {"genes": y.index, "y": y, "meta": meta, "libraries": manifest.reset_index(),
            "sf": pd.Series(sf, index=y.columns), "counts": c}


def pooled(a, b):
    common = a["y"].index.intersection(b["y"].index)
    return {"y": pd.concat([a["y"].loc[common], b["y"].loc[common]], axis=1),
            "meta": pd.concat([a["meta"], b["meta"]])}


def nuisance(meta):
    columns = []
    if meta.age.nunique() > 1:
        columns.append((meta.age == "29").to_numpy(float))
    if meta.tissue.nunique() > 1:
        columns.extend(pd.get_dummies(meta.tissue, drop_first=True).to_numpy(float).T)
    return np.column_stack(columns) if columns else None


def indices(genes, annotation):
    # Avoid silent alias substitution and document missing panel genes.
    mapping = annotation.set_index("ensembl_id").gene_symbol.to_dict()
    return {p: np.flatnonzero(np.isin([mapping.get(g, "") for g in genes], syms))
            for p, syms in PANELS.items()}


def analyse_contrast(name, data, mask, group_field, positive, annotation,
                     result_rows, gradient_rows, gene_rows, plans, sensitivity="design_cohort"):
    meta = data["meta"].loc[mask].copy()
    y = data["y"].loc[:, meta.index]
    z = nuisance(meta)
    group = (meta[group_field] == positive).to_numpy(int)
    blocks = ("age",) if meta.age.nunique() > 1 else ()
    labels, method, total = animal_permutations(meta, group, blocks)
    c = Contrast(y.to_numpy(), z)
    beta, se, t, pgene, df = c.cluster_fit(group, meta.animal.to_numpy())
    panel = indices(y.index, annotation)
    union = np.unique(np.concatenate(list(panel.values())))
    small = Contrast(y.to_numpy()[union], z)
    perm_beta = small.permutation_betas(labels)
    position = {i: j for j, i in enumerate(union)}
    panel_null = {}
    # Percentile cluster bootstrap intervals, conditional on normalization and
    # sampling exchangeability within design label and age. All animal regions
    # are resampled together; these intervals do not establish causal effects.
    panel_names = list(panel)
    averages = np.vstack([y.to_numpy()[panel[p]].mean(axis=0) for p in panel_names])
    rng = np.random.default_rng(561562)
    per = meta.assign(label=group).drop_duplicates("animal")
    strata = [f.animal.tolist() for _, f in per.groupby(["label", "age"], sort=True)]
    boots = []
    for _ in range(1000):
        chosen = list(itertools.chain.from_iterable(rng.choice(s, len(s), replace=True) for s in strata))
        idx = np.concatenate([np.flatnonzero(meta.animal.to_numpy() == a) for a in chosen])
        boots.append(Contrast(averages[:, idx], None if z is None else z[idx]).fit(group[idx])[0])
    boots = np.array(boots)
    for p, idx in panel.items():
        null = perm_beta[[position[i] for i in idx]].mean(axis=0)
        panel_null[p] = null
        obs = beta[idx].mean()
        result_rows.append({"contrast": name, "sensitivity": sensitivity, "panel": p,
                            "n_genes": len(idx), "n_libraries_after_collapse": len(meta),
                            "n_animals": meta.animal.nunique(), "n_positive_animals": meta.loc[group == 1].animal.nunique(),
                            "mean_log2_effect": obs, "mean_cluster_t_exploratory": t[idx].mean(),
                            "bootstrap_lower95": float(np.quantile(boots[:, panel_names.index(p)], .025)),
                            "bootstrap_upper95": float(np.quantile(boots[:, panel_names.index(p)], .975)),
                            "n_animal_bootstraps": len(boots),
                            "p_animal_permutation": float(permutation_p(obs, null, method)),
                            "permutation_method": method, "n_assignments_evaluated": labels.shape[1],
                            "total_possible_assignments": total, "n_conflicting_libraries": int(meta.location_conflict.sum())})
    obs = beta[panel["rPRG"]].mean() - beta[panel["SRG"]].mean()
    null = panel_null["rPRG"] - panel_null["SRG"]
    gradient_rows.append({"contrast": name, "sensitivity": sensitivity,
                          "rPRG_minus_SRG_mean_log2_effect": obs,
                          "bootstrap_lower95": float(np.quantile(boots[:, panel_names.index("rPRG")] - boots[:, panel_names.index("SRG")], .025)),
                          "bootstrap_upper95": float(np.quantile(boots[:, panel_names.index("rPRG")] - boots[:, panel_names.index("SRG")], .975)),
                          "p_animal_permutation": float(permutation_p(obs, null, method))})
    plans.append({"contrast": name, "sensitivity": sensitivity, "positive_level": positive,
                  "factor": group_field, "blocks": ",".join(blocks),
                  "animals": "|".join(sorted(meta.animal.unique())),
                  "n_animals": meta.animal.nunique(), "n_samples": len(meta),
                  "causal_identification": "not established; observational cohort or exposure contrast"})
    genes = pd.DataFrame({"contrast": name, "sensitivity": sensitivity, "ensembl_id": y.index,
                          "beta_log2_effect": beta, "cluster_CR1_se": se,
                          "cluster_t": t, "p_cluster_approximate": pgene,
                          "q_BH_within_contrast": bh(pgene), "cluster_df": df})
    gene_rows.append(genes)


def anchor_evaluation(name, training, train_mask, train_field, train_positive,
                      target, target_mask, target_field, target_positive, annotation,
                      output, module_output, bootstrap=200):
    tm = training["meta"].loc[train_mask]
    em = target["meta"].loc[target_mask]
    common = training["y"].index.intersection(target["y"].index)
    yt = training["y"].loc[common, tm.index].to_numpy()
    ye = target["y"].loc[common, em.index].to_numpy()
    zt, ze = nuisance(tm), nuisance(em)
    gt = (tm[train_field] == train_positive).to_numpy(int)
    ge = (em[target_field] == target_positive).to_numpy(int)
    panels = indices(common, annotation)
    sensor = panels["Sensor31"]
    anchor = learn_anchor(yt, gt, zt, sensor)
    beta_raw, traw = Contrast(ye, ze).fit(ge)
    outside = np.setdiff1d(np.arange(len(common)), sensor)
    energy0 = np.sum(beta_raw[outside]**2)
    panelenergy0 = np.sum(beta_raw[sensor]**2)
    overlap = sorted(set(tm.animal) & set(em.animal))
    blocks = ("age",) if em.age.nunique() > 1 else ()
    perm, perm_method, _ = animal_permutations(em, ge, blocks)
    union = np.unique(np.concatenate(list(panels.values())))
    pos = {i: j for j, i in enumerate(union)}
    summaries = []
    for mode in ("raw", "projection", "score", "bounded_score"):
        corrected = transform(ye, anchor, mode)
        fit = Contrast(corrected, ze)
        beta, t = fit.fit(ge)
        summaries.append({"comparison": name, "mode": mode,
                          "n_anchor_animals": tm.animal.nunique(), "n_target_animals": em.animal.nunique(),
                          "n_shared_animals": len(overlap), "shared_animals": "|".join(overlap),
                          "n_aligned_genes": len(common), "n_sensor_genes": len(sensor),
                          "sensor_effect_energy_removed_pct": 100 * (1 - np.sum(beta[sensor]**2) / panelenergy0),
                          "outside_sensor_effect_energy_removed_pct": 100 * (1 - np.sum(beta[outside]**2) / energy0),
                          "lambda_diagnostic_OLS_conditional_on_anchor": lambda_diagnostic(t[outside], fit.df),
                          "inference_scope": "fixed-anchor sensitivity; anchor estimation uncertainty not in conditional t"})
        perm_beta = Contrast(corrected[union], ze).permutation_betas(perm)
        for panel, idx in panels.items():
            obs = beta[idx].mean()
            null = perm_beta[[pos[i] for i in idx]].mean(axis=0)
            module_output.append({"comparison": name, "mode": mode, "panel": panel,
                                  "mean_log2_effect": obs, "n_genes": len(idx),
                                  "p_permutation_conditional_on_fixed_anchor": float(permutation_p(obs, null, perm_method)),
                                  "target_overlaps_anchor": bool(overlap)})
    # Cluster bootstrap ground animals within cohort and age, keeping all regions.
    # Quantifies anchor uncertainty conditional on normalized counts and fixed target.
    rng = np.random.default_rng(561562)
    per = tm.assign(label=gt).drop_duplicates("animal")
    strata = [f.animal.tolist() for _, f in per.groupby(["label", "age"], sort=True)]
    effects = {"projection": [], "score": [], "bounded_score": []}
    for _ in range(bootstrap):
        selected = list(itertools.chain.from_iterable(rng.choice(s, len(s), replace=True) for s in strata))
        idx = np.concatenate([np.flatnonzero(tm.animal.to_numpy() == a) for a in selected])
        try:
            fitted = learn_anchor(yt[:, idx], gt[idx], None if zt is None else zt[idx], sensor)
        except ValueError:
            continue
        for mode in effects:
            beta, _ = Contrast(transform(ye, fitted, mode), ze).fit(ge)
            effects[mode].append(100 * (1 - np.sum(beta[outside]**2) / energy0))
    for row in summaries:
        values = effects.get(row["mode"], [])
        row["n_anchor_bootstraps"] = len(values)
        row["anchor_bootstrap_energy_lower95"] = np.quantile(values, .025) if values else np.nan
        row["anchor_bootstrap_energy_upper95"] = np.quantile(values, .975) if values else np.nan
    output.extend(summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=200)
    args = parser.parse_args()
    check_input_integrity()
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    annotation = pd.read_csv(HERE / "inputs/gene_annotation.csv.gz")
    datasets = {n: load_dataset(n) for n in DATASETS}
    libs = pd.concat([d["libraries"] for d in datasets.values()], ignore_index=True)
    libs.to_csv(out / "sample_manifest.csv", index=False)
    conflicts = libs.loc[libs.location_conflict].copy()
    conflicts["adjudication"] = "unresolved: do not infer terminal location from design label"
    conflicts.to_csv(out / "metadata_conflicts.csv", index=False)
    animals = pd.concat([d["meta"] for d in datasets.values()])
    animals.to_csv(out / "animal_region_manifest.csv", index_label="analysis_sample")
    pd.concat([d["sf"].rename(f"OSD-{n}") for n, d in datasets.items()]).rename("size_factor").to_csv(out / "size_factors.csv")
    pd.DataFrame([{"accession": f"OSD-{n}", "mission": DATASETS[n][0], "tissue": DATASETS[n][1],
                   "n_deposited_libraries_in_count_matrix": len(d["libraries"]),
                   "n_animals_in_region": len(d["meta"]), "n_filtered_genes": len(d["y"])}
                  for n, d in datasets.items()]).to_csv(out / "dataset_accounting.csv", index=False)
    panel_rows = [{"panel": p, "symbol": symbol, "legacy_defined_exploratory_panel": True,
                   "evidence_scope": "gene-level expression; no direct kinetic, isoform or inflammatory phenotype measurement"}
                  for p, syms in PANELS.items() for symbol in syms]
    pd.DataFrame(panel_rows).to_csv(out / "panel_definitions.csv", index=False)
    pool = pooled(datasets["561"], datasets["562"])
    rows, gradients, genes, plans = [], [], [], []
    conflict_animals = set(conflicts.animal)
    for label, data in [("RRRM2_pooled", pool), ("OSD561_cerebellum", datasets["561"]),
                        ("OSD562_hippocampus", datasets["562"])]:
        for flight in ("GC", "FLT"):
            mask = data["meta"].flight == flight
            name = f"{label}_{flight}_ISST_vs_LAR"
            analyse_contrast(name, data, mask, "cohort", "ISS_T", annotation, rows, gradients, genes, plans)
            safe = mask & ~data["meta"].animal.isin(conflict_animals)
            analyse_contrast(name, data, safe, "cohort", "ISS_T", annotation, rows, gradients, genes, plans,
                             "exclude_all_regions_of_conflicting_animals")
        # Location-field sensitivity is only meaningful within a single region;
        # contradictory animal-level location labels cannot be pooled consistently.
        if label != "RRRM2_pooled":
            for cohort in ("ISS_T", "LAR"):
                analyse_contrast(f"{label}_FLT_vs_GC_within_{cohort}", data,
                                 data["meta"].cohort == cohort, "flight", "FLT",
                                 annotation, rows, gradients, genes, plans)
                analyse_contrast(f"{label}_FLT_vs_GC_within_{cohort}", data,
                                 (data["meta"].cohort == cohort) & ~data["meta"].animal.isin(conflict_animals),
                                 "flight", "FLT", annotation, rows, gradients, genes, plans,
                                 "exclude_all_regions_of_conflicting_animals")
            analyse_contrast(f"{label}_FLT_location_ISS_vs_Earth", data, data["meta"].flight == "FLT",
                             "deposited_location", "On ISS", annotation, rows, gradients, genes, plans,
                             "deposited_location_within_region")
    for n in ("563", "564", "612", "613"):
        data = datasets[n]
        analyse_contrast(f"OSD{n}_FLT_vs_GC", data, np.ones(len(data["meta"]), bool), "flight", "FLT",
                         annotation, rows, gradients, genes, plans)
    d525 = datasets["525"]
    for a, b in [("HGC", "VGC"), ("FLT", "HGC"), ("FLT", "VGC")]:
        analyse_contrast(f"OSD525_{a}_vs_{b}", d525, d525["meta"].flight.isin([a, b]), "flight", a,
                         annotation, rows, gradients, genes, plans)
    for filename, frame in [("panel_results.csv", pd.DataFrame(rows)), ("gradient_results.csv", pd.DataFrame(gradients))]:
        # Two separate families: all panel tests; all rapid-vs-secondary gradients.
        frame["q_BH_all_tests_in_family"] = bh(frame.p_animal_permutation)
        frame.to_csv(out / filename, index=False)
    pd.concat(genes, ignore_index=True).merge(annotation, on="ensembl_id", how="left").to_csv(
        out / "gene_results_exploratory.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    pd.DataFrame(plans).to_csv(out / "contrast_manifest.csv", index=False)
    print(f"Reanalysed {len(plans)} contrasts/sensitivities, {len(rows)} panel tests", flush=True)

    evals, modules = [], []
    for train_n in ("561", "562"):
        training = datasets[train_n]
        for target_n in ("561", "562", "563", "564", "612", "613"):
            target = datasets[target_n]
            recovery = target_n in ("561", "562")
            anchor_evaluation(f"anchor_OSD{train_n}_to_OSD{target_n}" + ("_FLT_recovery" if recovery else "_FLT_vs_GC"),
                              training, training["meta"].flight == "GC", "cohort", "ISS_T",
                              target, target["meta"].flight == "FLT" if recovery else np.ones(len(target["meta"]), bool),
                              "cohort" if recovery else "flight", "ISS_T" if recovery else "FLT",
                              annotation, evals, modules, bootstrap=args.bootstrap)
    # RR-23 within-cohort anchor shares HGC animals with target; mark explicitly.
    anchor_evaluation("OSD525_shared_HGC_anchor", d525, d525["meta"].flight.isin(["HGC", "VGC"]), "flight", "HGC",
                      d525, d525["meta"].flight.isin(["FLT", "HGC"]), "flight", "FLT",
                      annotation, evals, modules, bootstrap=args.bootstrap)
    for train_n in ("561", "562"):
        training = datasets[train_n]
        anchor_evaluation(f"anchor_OSD{train_n}_to_OSD525", training, training["meta"].flight == "GC", "cohort", "ISS_T",
                          d525, d525["meta"].flight.isin(["FLT", "HGC"]), "flight", "FLT",
                          annotation, evals, modules, bootstrap=args.bootstrap)
    pd.DataFrame(evals).to_csv(out / "ground_anchor_sensitivity.csv", index=False)
    module_frame = pd.DataFrame(modules)
    module_frame["q_BH_conditional_sensitivity_family"] = bh(module_frame.p_permutation_conditional_on_fixed_anchor)
    module_frame.to_csv(out / "ground_anchor_panel_sensitivity.csv", index=False)
    print("Ground anchor comparisons complete", flush=True)

    # 126 overlapping partitions of the same nine HGC animals. No independent-n claim.
    split_rows = []
    y = d525["y"].to_numpy()
    m = d525["meta"]
    sensor = indices(d525["y"].index, annotation)["Sensor31"]
    outside = np.setdiff1d(np.arange(len(y)), sensor)
    hgc = np.flatnonzero(m.flight.to_numpy() == "HGC")
    vgc = np.flatnonzero(m.flight.to_numpy() == "VGC")
    flt = np.flatnonzero(m.flight.to_numpy() == "FLT")
    for split, a in enumerate(itertools.combinations(hgc, 4), 1):
        a = np.array(a)
        b = np.setdiff1d(hgc, a)
        ti, ei = np.r_[a, vgc], np.r_[flt, b]
        anchor = learn_anchor(y[:, ti], np.r_[np.ones(len(a)), np.zeros(len(vgc))], None, sensor)
        ge = np.r_[np.ones(len(flt)), np.zeros(len(b))]
        beta0, _ = Contrast(y[:, ei]).fit(ge)
        for mode in ("projection", "score", "bounded_score"):
            fit = Contrast(transform(y[:, ei], anchor, mode))
            beta, t = fit.fit(ge)
            split_rows.append({"split": split, "mode": mode, "training_HGC_animals": "|".join(m.iloc[a].animal),
                               "evaluation_HGC_animals": "|".join(m.iloc[b].animal), "n_shared_animals": 0,
                               "sensor_energy_removed_pct": 100 * (1 - np.sum(beta[sensor]**2) / np.sum(beta0[sensor]**2)),
                               "outside_sensor_energy_removed_pct": 100 * (1 - np.sum(beta[outside]**2) / np.sum(beta0[outside]**2)),
                               "lambda_conditional_diagnostic": lambda_diagnostic(t[outside], fit.df),
                               "independent_replication_unit": "one RR-23 cohort; partitions overlap"})
    pd.DataFrame(split_rows).to_csv(out / "OSD525_disjoint_anchor_partitions.csv", index=False)

    # Same-library estimator check, using all 29 libraries for both estimators.
    rsem = load_dataset("525", "RSEM")
    common = d525["y"].index.intersection(rsem["y"].index)
    g = (m.flight == "HGC").to_numpy(int)
    use = m.flight.isin(["HGC", "VGC"]).to_numpy()
    bstar, _ = Contrast(d525["y"].loc[common].to_numpy()[:, use]).fit(g[use])
    brsem, _ = Contrast(rsem["y"].loc[common].to_numpy()[:, use]).fit(g[use])
    idx = indices(common, annotation)["Sensor31"]
    pd.DataFrame([{"comparison": "OSD525_STAR_vs_RSEM_same_29_libraries", "n_collapsed_animals": len(m),
                   "sensor_beta_pearson_r": np.corrcoef(bstar[idx], brsem[idx])[0, 1],
                   "genome_beta_pearson_r": np.corrcoef(bstar, brsem)[0, 1],
                   "independent_validation": False}]).to_csv(out / "estimator_sensitivity.csv", index=False)
    run = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
           "scipy": scipy.__version__, "seed": 561562, "bootstrap": args.bootstrap,
           "datasets": list(DATASETS), "n_libraries": len(libs), "n_unique_source_animals": libs.animal.nunique(),
           "core_561_562_libraries": len(pool["meta"]), "core_561_562_animals": pool["meta"].animal.nunique(),
           "location_conflict_libraries": len(conflicts), "location_conflict_animals": len(conflict_animals),
           "status": "exploratory revision; unresolved metadata and external validation remain"}
    (out / "run_summary.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
