"""Post-comparison DIAGNOSTIC only (written after reading Supplementary Table 2/3).
Tests which preprocessing variant explains the ~0.001 log2 differences. Does not change the
frozen primary numbers in G_results.json."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import g_common as gc

TARGET = {  # (region, contrast, panel): paper value from SI Tables 2-3
    ("CB", "DOM", "rPRG"): 0.503, ("CB", "DOM", "HSP"): 0.210, ("CB", "DOM", "Glia_Complement"): -0.314,
    ("CB", "P_G", "rPRG"): 0.704, ("CB", "F_L", "rPRG"): 0.188, ("CB", "F_I", "rPRG"): -0.206,
    ("CB", "P_F", "rPRG"): 0.332, ("HPC", "F_L", "rPRG"): 0.281, ("HPC", "F_I", "rPRG"): 0.164,
    ("HPC", "P_F", "rPRG"): 1.004, ("HPC", "DOM", "rPRG"): 0.837,
}
TP = {("CB", "DOM", "dPRG"): 0.44, ("CB", "DOM", "HSP"): 0.19}
C = {"DOM": (("GC", "ISS_T"), ("FLT", "LAR")), "P_G": (("GC", "ISS_T"), ("GC", "LAR")),
     "F_L": (("FLT", "LAR"), ("GC", "LAR")), "F_I": (("FLT", "ISS_T"), ("GC", "ISS_T")),
     "P_F": (("FLT", "ISS_T"), ("FLT", "LAR"))}


def normalise(cnt, variant):
    keep = cnt.mean(axis=1) >= 10
    sfbase = cnt if variant in ("sf_prefilter", "sf_prefilter_geo1") else cnt.loc[keep]
    pos = (sfbase > 0).all(axis=1)
    lr = np.log(sfbase.loc[pos])
    sf = np.exp(lr.sub(lr.mean(axis=1), axis=0).median(axis=0))
    if variant.endswith("geo1"):
        sf = sf / np.exp(np.log(sf).mean())
    c = cnt.loc[keep]
    return np.log2(c.div(sf, axis=1) + 1)


PAN = gc.panel_ids()
for variant in ("primary", "geo1", "sf_prefilter", "sf_prefilter_geo1", "per_contrast_subset"):
    rows = []
    for (acc, tis) in (("OSD-561", "CB"), ("OSD-562", "HPC")):
        cnt, am = gc.load_deposit(acc)
        am["animal_id"] = am.index
        for (reg, con, pan), val in {**TARGET, **{k: None for k in TP}}.items():
            if reg != tis:
                continue
            g1, g0 = C[con]
            s1 = (am.flight == g1[0]) & (am.cohort == g1[1])
            s0 = (am.flight == g0[0]) & (am.cohort == g0[1])
            sub = am[s1 | s0]
            if variant == "per_contrast_subset":
                Y = normalise(cnt[sub.index], "primary")
            else:
                Y = normalise(cnt, "geo1" if variant == "geo1" else variant)[sub.index]
            r = gc.exact_panel_test(Y, sub, s1[sub.index].astype(float).values, {pan: PAN[pan]}, False)[pan]
            rows.append((reg, con, pan, round(r["effect"], 4), val, round(r["p_exact"], 4), TP.get((reg, con, pan))))
    df = pd.DataFrame(rows, columns=["reg", "con", "panel", "effect", "paper", "p", "paper_p"])
    ok = (df.dropna(subset=["paper"]).eval("abs(effect-paper)") < 0.0005 + 1e-9).sum()
    print(f"== {variant}: {ok}/{df.paper.notna().sum()} effects match at 3 dp")
    print(df.to_string(index=False))
