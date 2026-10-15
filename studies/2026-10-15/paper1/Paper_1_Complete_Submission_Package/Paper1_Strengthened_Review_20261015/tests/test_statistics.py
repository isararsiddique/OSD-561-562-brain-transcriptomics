"""Independent regression references and scientific misuse regression checks."""
import itertools
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from statistics_core import (Contrast, animal_permutations, bh, bounded_gamma, learn_anchor,
                             permutation_p, transform)
from run_revision import check_input_integrity, load_dataset


class StatisticalChecks(unittest.TestCase):
    def test_fwl_matches_full_design_lstsq_with_rank_deficient_nuisance(self):
        rng = np.random.default_rng(932)
        y = rng.normal(size=(12, 15))
        g = np.r_[np.zeros(7), np.ones(8)]
        z = rng.normal(size=(15, 1))
        duplicated = np.column_stack([z, z, np.ones(15)])
        beta, _ = Contrast(y, duplicated).fit(g)
        x = np.column_stack([np.ones(15), z, g])
        reference = np.linalg.lstsq(x, y.T, rcond=None)[0][-1]
        np.testing.assert_allclose(beta, reference, atol=1e-12)
        self.assertEqual(Contrast(y, duplicated).df, 12)

    def test_cluster_sandwich_matches_full_matrix_calculation(self):
        rng = np.random.default_rng(33)
        n = 18
        animal = np.repeat(np.arange(9), 2)
        g = np.repeat(np.r_[np.zeros(4), np.ones(5)], 2)
        z = np.tile([0, 1], 9)[:, None]
        y = rng.normal(size=(3, n)) + rng.normal(size=(3, 9))[:, animal]
        beta, se, _, _, _ = Contrast(y, z).cluster_fit(g, animal)
        x = np.column_stack([np.ones(n), z, g])
        coefs = np.linalg.lstsq(x, y.T, rcond=None)[0]
        resid = y.T - x @ coefs
        bread = np.linalg.inv(x.T @ x)
        expected = []
        for gene in range(3):
            meat = np.zeros((3, 3))
            for a in np.unique(animal):
                u = x[animal == a].T @ resid[animal == a, gene]
                meat += np.outer(u, u)
            expected.append(np.sqrt((bread @ meat @ bread)[-1, -1] * 9 / 8 * 17 / 15))
        np.testing.assert_allclose(beta, coefs[-1], atol=1e-12)
        np.testing.assert_allclose(se, expected, atol=1e-12)

    def test_permutations_keep_regions_together_and_preserve_age_counts(self):
        meta = pd.DataFrame({"animal": list("aabbccddeeff"), "age": [12] * 6 + [29] * 6})
        g = np.repeat([1, 0, 0, 1, 1, 0], 2)
        labels, method, total = animal_permutations(meta, g)
        self.assertEqual(total, 9)
        self.assertTrue(method.startswith("complete"))
        np.testing.assert_array_equal(labels[::2], labels[1::2])
        np.testing.assert_array_equal(labels[:6].sum(axis=0), np.full(9, 2))
        np.testing.assert_array_equal(labels[6:].sum(axis=0), np.full(9, 4))
        self.assertTrue(np.any(np.all(labels == g[:, None], axis=0)))
        g[1] = 0
        with self.assertRaisesRegex(ValueError, "Inconsistent label"):
            animal_permutations(meta, g)

    def test_unidentifiable_contrast_errors_instead_of_zero_effect(self):
        with self.assertRaisesRegex(ValueError, "not identifiable"):
            Contrast(np.arange(24).reshape(3, 8)).fit(np.ones(8))

    def test_permutation_p_values_have_correct_finite_resolution(self):
        self.assertEqual(permutation_p(3, np.array([3., 2, 1]), "complete"), 1 / 3)
        self.assertEqual(permutation_p(4, np.array([3., 2, 1]), "Monte Carlo"), 1 / 4)
        np.testing.assert_allclose(bh([.01, .02, .9]), [.03, .03, .9])

    def test_expression_projection_projects_beta_but_not_t(self):
        rng = np.random.default_rng(12)
        y = rng.normal(size=(20, 18))
        g = np.r_[np.zeros(9), np.ones(9)]
        v = rng.normal(size=20)
        v /= np.linalg.norm(v)
        anchor = (v, v, v)
        y2 = transform(y, anchor, "projection")
        b, t = Contrast(y).fit(g)
        b2, t2 = Contrast(y2).fit(g)
        np.testing.assert_allclose(b2, b - v * (v @ b), atol=1e-12)
        np.testing.assert_allclose(v @ y2, np.zeros(18), atol=1e-12)
        self.assertGreater(np.linalg.norm(t2 - (t - v * (v @ t))), 1e-3)

    def test_score_operator_is_oblique_and_removes_sensor_score(self):
        rng = np.random.default_rng(13)
        y = rng.normal(size=(20, 18))
        g = np.r_[np.zeros(9), np.ones(9)]
        anchor = learn_anchor(y, g, None, np.arange(5))
        v, w, gamma = anchor
        self.assertAlmostEqual(float(w @ gamma), 1, places=12)
        y2 = transform(y, anchor, "score")
        np.testing.assert_allclose(w @ y2, np.zeros(18), atol=1e-12)
        self.assertGreater(np.linalg.norm(v @ y2), 1e-3)

    def test_rejects_nonfinite_expression_and_nuisance(self):
        with self.assertRaisesRegex(ValueError, "finite gene-by-sample"):
            Contrast(np.full((3, 8), np.nan))
        with self.assertRaisesRegex(ValueError, "Nuisance design must be finite"):
            Contrast(np.ones((3, 8)), np.full((8, 1), np.inf))

    def test_exact_fit_nonzero_effect_is_not_zero_statistic(self):
        g = np.r_[np.zeros(4), np.ones(4)]
        _, t = Contrast(np.vstack([g, -g, np.zeros(8)])).fit(g)
        self.assertTrue(np.isposinf(t[0]) or t[0] > 1e7)
        self.assertTrue(np.isneginf(t[1]) or t[1] < -1e7)
        self.assertEqual(t[2], 0)

    def test_bounded_operator_dense_svd_and_sensor_null_effect(self):
        rng = np.random.default_rng(762)
        for n in [3, 10, 50]:
            w = rng.normal(size=n); w /= np.linalg.norm(w)
            h = rng.normal(size=n); h -= w * (w @ h); h *= 20 / np.linalg.norm(h)
            gl = bounded_gamma(w, w+h)
            a = np.eye(n)-np.outer(gl, w)
            self.assertAlmostEqual(np.linalg.svd(a, compute_uv=False)[0], 2, places=11)
            np.testing.assert_allclose(a@a, a, atol=1e-12)
            b = rng.normal(size=n); b -= w*(w@b)
            np.testing.assert_allclose(a@b, b, atol=1e-12)

    def test_core_manifest_and_both_estimators_collapse_same_libraries(self):
        check_input_integrity()
        a, b = load_dataset("561"), load_dataset("562")
        self.assertEqual(len(a["libraries"]) + len(b["libraries"]), 53)
        self.assertEqual(len(set(a["meta"].animal) | set(b["meta"].animal)), 29)
        star, rsem = load_dataset("525"), load_dataset("525", "RSEM")
        self.assertEqual(len(star["libraries"]), 29)
        self.assertEqual(len(star["meta"]), 26)
        np.testing.assert_array_equal(star["meta"].n_libraries, rsem["meta"].n_libraries)
        self.assertEqual(star["meta"].n_libraries.sum(), 29)


if __name__ == "__main__":
    unittest.main()
