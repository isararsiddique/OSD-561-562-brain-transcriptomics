#!/usr/bin/env python3
"""Independent dense-matrix and geometric tests of the claimed operators.

Expected values are calculated from matrix products, dense SVD and the
projection of a Euclidean ball, rather than from the production formula.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from statistics_core import bounded_gamma, transform  # noqa: E402
from run_validation_extension import spectral_gain  # noqa: E402


def main():
    rng = np.random.default_rng(151061)
    errors = {name: [] for name in [
        "dense_svd_vs_gamma_norm", "production_spectral_gain_vs_dense_svd",
        "bounded_spectral_gain_above_cap", "bounded_idempotence",
        "weighted_sensor_removal", "sensor_null_biology_preservation",
        "production_transform_vs_dense", "bounded_gamma_vs_ball_projection",
        "noise_eigenvalue_bound_excess", "outside_sensor_norm_identity",
        "partial_operator_svd_vs_formula",
    ]}
    examples = []
    trials = 0
    for dim in (2, 3, 8, 31, 80):
        for length in (0., .3, 1., 8., 20.):
            for cap in (1., 1.1, 2., 10.):
                w = rng.normal(size=dim)
                w /= np.linalg.norm(w)
                z = rng.normal(size=dim)
                z -= w * (w @ z)
                z *= length / np.linalg.norm(z) if length else 0
                gamma = w + z
                A = np.eye(dim) - np.outer(gamma, w)
                dense_gain = np.linalg.svd(A, compute_uv=False)[0]
                errors["dense_svd_vs_gamma_norm"].append(abs(dense_gain - np.linalg.norm(gamma)))
                errors["production_spectral_gain_vs_dense_svd"].append(abs(spectral_gain(w, gamma) - dense_gain))
                gb = bounded_gamma(w, gamma, max_gain=cap)
                B = np.eye(dim) - np.outer(gb, w)
                gainb = np.linalg.svd(B, compute_uv=False)[0]
                errors["bounded_spectral_gain_above_cap"].append(max(0., gainb-cap))
                errors["bounded_idempotence"].append(np.linalg.norm(B @ B-B, ord="fro"))
                errors["weighted_sensor_removal"].append(np.linalg.norm(w @ B))
                biology = rng.normal(size=dim)
                biology -= w * (w @ biology)
                errors["sensor_null_biology_preservation"].append(np.linalg.norm(B @ biology-biology))
                y = rng.normal(size=(dim, 11))
                v = rng.normal(size=dim)
                v /= np.linalg.norm(v)
                if cap == 2:
                    errors["production_transform_vs_dense"].append(np.linalg.norm(transform(y, (v,w,gamma), "bounded_score")-B@y))
                # Orthogonal projection onto ball in the affine hyperplane w'c=1.
                r = np.sqrt(cap*cap-1)
                expected = w + (z if length <= r else z*r/length)
                errors["bounded_gamma_vs_ball_projection"].append(np.linalg.norm(gb-expected))
                U = rng.normal(size=(dim, dim))
                covariance = U@U.T + np.eye(dim)*.1
                before = np.linalg.eigvalsh(covariance)[-1]
                after = np.linalg.eigvalsh(B @ covariance @ B.T)[-1]
                errors["noise_eigenvalue_bound_excess"].append(max(0., after-cap*cap*before))
                for alpha in (.1,.5,1):
                    partial = np.eye(dim)-alpha*np.outer(gamma,w)
                    errors["partial_operator_svd_vs_formula"].append(abs(spectral_gain(w,gamma,alpha)-np.linalg.svd(partial,compute_uv=False)[0]))
                trials += 1
    # Sensor supported on the first two coordinates: the outside-row norm has
    # an independent exact expression sqrt(1 + ||gamma_outside||^2).
    for dim in (3,8,31):
        w = np.zeros(dim)
        w[:2] = [np.sqrt(.3),np.sqrt(.7)]
        gamma = rng.normal(size=dim)
        gamma += w*(1-w@gamma)
        A = np.eye(dim)-np.outer(gamma,w)
        outside = A[2:]
        gain = np.linalg.svd(outside,compute_uv=False)[0]
        errors["outside_sensor_norm_identity"].append(abs(gain-np.sqrt(1+np.dot(gamma[2:],gamma[2:]))))
    # Counterexample: cancellation of the scalar sensor score does not mean
    # all sensor coordinates are damped or that biology is preserved.
    w = np.array([1.,0.,0.])
    gamma = np.array([1.,2.,3.])
    A = np.eye(3)-np.outer(gamma,w)
    b = np.array([1.,0.,0.])
    examples.append({"description":"A removes weighted sensor score but increases sensor-panel energy", "sensor_indices":[0,1], "b":b.tolist(), "A_b":(A@b).tolist(), "sensor_energy_before":float(np.dot(b[:2],b[:2])), "sensor_energy_after":float(np.dot((A@b)[:2],(A@b)[:2])), "weighted_sensor_after":float(w@A@b)})
    # In dimension one the only admissible idempotent projection is zero;
    # the G>=2 qualification cannot be omitted from ||A||=||gamma||.
    one_d = float(np.linalg.svd(np.zeros((1,1)),compute_uv=False)[0])
    report = {"seed":151061,"random_geometry_trials":trials,"gene_dimensions":[2,3,8,31,80],"caps":[1.,1.1,2.,10.],"maximum_absolute_errors":{k:float(max(v)) if v else None for k,v in errors.items()},"counterexamples":examples,"one_dimensional_operator_norm":one_d,"all_numerical_identities_pass":all(not v or max(v)<1e-8 for v in errors.values()),"interpretation":"Matrix identities and norm bounds only. They do not prove causal identification, signal recovery or biological FDR control."}
    (HERE/"independent_operator_test_results.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    if not report["all_numerical_identities_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
