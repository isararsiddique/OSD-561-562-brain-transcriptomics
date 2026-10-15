#!/usr/bin/env python3
"""Check the fixed-sensor isotropic-noise high-dimensional gain identity."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(231200)
G,N,R = 1200,24,1000
r = N-1
values=[]
for _ in range(R):
    y=rng.normal(size=(G,N))
    y-=y.mean(axis=1,keepdims=True)
    s=y[0]
    gamma=y@s/(s@s)
    values.append(float(gamma@gamma))
values=np.asarray(values)
theoretical=1+(G-1)/(r-2)
variance=2*(G-1)*(G-1+r-2)/((r-2)**2*(r-4))
report=dict(seed=231200,genes=G,samples=N,residual_rank=r,repetitions=R,
            fixed_sensor_description="Unit first-coordinate sensor, fixed before data generation",
            population_gamma_l2_norm=1.,theoretical_mean_estimated_gain_squared=theoretical,
            empirical_mean_estimated_gain_squared=float(values.mean()),
            theoretical_mc_se=float(np.sqrt(variance/R)),empirical_mc_se=float(values.std(ddof=1)/np.sqrt(R)),
            mean_agreement_standard_errors=float(abs(values.mean()-theoretical)/np.sqrt(variance/R)),
            within_five_theoretical_standard_errors=bool(abs(values.mean()-theoretical)<5*np.sqrt(variance/R)),
            limitation="The expectation assumes fixed w independent of the Gaussian data. A sensor estimated from ground contrasts has a different sampling law.")
(HERE/"fixed_sensor_noise_inflation_results.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
if not report["within_five_theoretical_standard_errors"]:raise SystemExit(1)
