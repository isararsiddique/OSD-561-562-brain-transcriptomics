"""Parts 1 and 3: dominance test, four RRRM-2 contrasts, and published signatures (pooled)."""
import json
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from g_common import BUILD, PANELS, bh, dump, exact_panel_test, filter_and_normalise, load_deposit, panel_ids, symbol_map

regions = {}
for acc, tis in (("OSD-561", "CB"), ("OSD-562", "HPC")):
    cnt, am = load_deposit(acc)
    expr, sf = filter_and_normalise(cnt)
    am = am.copy()
    am["tissue"] = tis
    am["animal_id"] = am.index
    am.index = [f"{tis}|{a}" for a in am.index]
    expr.columns = am.index
    regions[tis] = (expr, am)
    print(acc, "animals", len(am), "genes kept", expr.shape[0], "libs>1:", int((am.n_libraries > 1).sum()))

common = regions["CB"][0].index.intersection(regions["HPC"][0].index)
pooled_expr = pd.concat([regions["CB"][0].loc[common], regions["HPC"][0].loc[common]], axis=1)
pooled_am = pd.concat([regions["CB"][1], regions["HPC"][1]])
regions["POOLED"] = (pooled_expr, pooled_am)
print("pooled genes", len(common))

PAN = {k: v for k, v in panel_ids().items() if k in PANELS}

CONTRASTS = {
    # name: (group1 (flight, cohort), group0)
    "DOM": (("GC", "ISS_T"), ("FLT", "LAR")),
    "P_G": (("GC", "ISS_T"), ("GC", "LAR")),
    "F_L": (("FLT", "LAR"), ("GC", "LAR")),
    "F_I": (("FLT", "ISS_T"), ("GC", "ISS_T")),
    "P_F": (("FLT", "ISS_T"), ("FLT", "LAR")),
}


def run(region, contrast, panels):
    expr, am = regions[region]
    g1, g0 = CONTRASTS[contrast]
    sel1 = (am.flight == g1[0]) & (am.cohort == g1[1])
    sel0 = (am.flight == g0[0]) & (am.cohort == g0[1])
    sub = am[sel1 | sel0]
    g = sel1[sub.index].astype(float).values
    return exact_panel_test(expr[sub.index], sub, g, panels, pooled=(region == "POOLED"))


results = {"part1_dominance": [], "part1_rRPG_contrasts": [], "part3": {}}
for region in ("CB", "HPC", "POOLED"):
    r = run(region, "DOM", PAN)
    for p in PANELS:
        results["part1_dominance"].append({"region": region, "panel": p, **r[p]})
q = bh([x["p_exact"] for x in results["part1_dominance"]])
for x, qq in zip(results["part1_dominance"], q):
    x["q_BH24"] = float(qq)

for region in ("CB", "HPC", "POOLED"):
    for c in ("P_G", "F_L", "F_I", "P_F"):
        r = run(region, c, {"rPRG": PAN["rPRG"]})["rPRG"]
        results["part1_rRPG_contrasts"].append({"region": region, "contrast": c, **r})

# ---------------------------------------------------------------- Part 3
sig = json.load(open(f"{BUILD}/research/F_signatures/signatures.json"))["signatures"]
sm = symbol_map()
for name, s in sig.items():
    variants = {}
    if s["map_by"] == "ensembl":
        all_ids, single_ids = [], []
        for gdef in s["genes"]:
            ids = re.findall(r"ENSMUSG\d+", str(gdef["ensembl"]))
            all_ids += ids
            if len(ids) == 1:
                single_ids += ids
        variants["primary_all_listed_ids"] = list(dict.fromkeys(all_ids))
        variants["sens_single_id_only"] = list(dict.fromkeys(single_ids))
    else:
        ids = []
        for gdef in s["genes"]:
            ids += sm.get(gdef["symbol"], [])
        variants["primary_exact_symbol"] = list(dict.fromkeys(ids))
    out = {"n_listed": len(s["genes"])}
    for vname, ids in variants.items():
        passing = [i for i in ids if i in regions["POOLED"][0].index]
        pan = {name: passing}
        rp = run("POOLED", "P_G", pan)[name]
        rd = run("POOLED", "DOM", pan)[name]
        out[vname] = {"n_ids": len(ids), "n_pass_pooled": len(passing),
                      "P_G": {k: rp[k] for k in ("effect", "p_exact", "n_assign")},
                      "DOM": {k: rd[k] for k in ("effect", "p_exact", "n_assign")}}
    results["part3"][name] = out

dump(results, "G_part13.json")
pd.DataFrame(results["part1_dominance"]).to_csv(f"{BUILD}/research/G_independent/G_part1_dominance.csv", index=False)
print(pd.DataFrame(results["part1_dominance"]).to_string())
print(pd.DataFrame(results["part1_rRPG_contrasts"]).to_string())
print(json.dumps(results["part3"], indent=1))
