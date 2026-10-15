"""Post-comparison DIAGNOSTIC for Part 2 (written after reading the Results text and Fig. 4b).
Variant A: S, N_tr, N_te estimated inside each split from TRAINING animals only, then averaged.
Variant B (RR-23): observed cosine with the gene filter fixed from all deposit animals.
Does not alter the frozen primary values in G_results.json."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import g_part2_lib as L  # noqa: E402

out = {}
for name, acc in L.SETTINGS.items():
    cnt, am, sp, lab = L.setup(acc)
    ec_w, snr_w, ec_s, snr_s, cw_B = [], [], [], [], []
    keep_all = cnt.mean(axis=1) >= 10
    for sid, sub in sp.groupby("split_id"):
        tr = list(sub[sub.role == "train"].animal)
        te = list(sub[sub.role == "test"].animal)
        Ytr = L.norm(cnt[tr])
        b, s2, c_tr = L.coef(Ytr, am.loc[tr, "age"].values, lab[tr].values)
        _, _, c_te = L.coef(L.norm(cnt[te], ref=cnt[tr]).iloc[:5], am.loc[te, "age"].values, lab[te].values)
        for idx, ec, sn in ((b.index, ec_w, snr_w), ([i for i in b.index if i in L.SENSOR], ec_s, snr_s)):
            bb, ss = b[idx].values, s2[idx].values
            Ntr, Nte = ss.sum() / c_tr, ss.sum() / c_te
            S = bb @ bb - Ntr
            sn.append(S / Ntr)
            ec.append(S / np.sqrt((S + Ntr) * (S + Nte)) if S > 0 else 0.0)
        if acc == "OSD-525":
            ctr, cte = cnt.loc[keep_all, tr], cnt.loc[keep_all, te]
            Ytr2 = L.norm_nofilter(ctr)
            Yte2 = L.norm_nofilter(cte, ref=ctr)
            b1, _, _ = L.coef(Ytr2, am.loc[tr, "age"].values, lab[tr].values)
            b2, _, _ = L.coef(Yte2, am.loc[te, "age"].values, lab[te].values)
            cw_B.append(L.cos(b1.values, b2.values))
    out[name] = {"A_expected_cos_whole": float(np.mean(ec_w)), "A_SNR_whole_mean": float(np.mean(snr_w)),
                 "A_SNR_whole_median": float(np.median(snr_w)),
                 "A_expected_cos_sensor": float(np.mean(ec_s)), "A_SNR_sensor_mean": float(np.mean(snr_s)),
                 "A_SNR_sensor_median": float(np.median(snr_s))}
    if cw_B:
        out[name]["B_observed_cos_whole_depositfilter"] = float(np.mean(cw_B))
    print(name, out[name])
pd.Series(out).to_json(L.OUT + "/G_diag_part2.json", indent=1)
