"""Copy this self-contained extension to a fresh folder and verify identical results."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import gzip
import hashlib
import numpy as np

B = Path(__file__).resolve().parents[1]


def main():
    target = B / "tmp/offline_reproduction_check"
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    for folder in ["scripts", "inputs"]:
        shutil.copytree(B / folder, target / folder)
    shutil.copy2(B / "RETROSPECTIVE_PLAN_2026-10-15.md", target / "RETROSPECTIVE_PLAN_2026-10-15.md")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for script in ["analyse.py", "verify_independently.py", "build_outputs.py"]:
        subprocess.run([sys.executable, "scripts/" + script], cwd=target, env=env, capture_output=True, text=True, check=True)
    assert json.loads((B / "tables/summary.json").read_text()) == json.loads((target / "tables/summary.json").read_text())
    assert json.loads((target / "tables/independent_numerical_checks.json").read_text())["all_checks_passed"]
    compared = []
    for candidate in sorted((target / "tables").glob("*.csv*")):
        if candidate.name in ["source_hashes_before_after.csv", "independent_numerical_checks.csv"]:
            continue
        original = B / "tables" / candidate.name
        reader = gzip.open if candidate.suffix == ".gz" else open
        with reader(original, "rb") as a, reader(candidate, "rb") as b:
            assert hashlib.sha256(a.read()).digest() == hashlib.sha256(b.read()).digest(), candidate.name
        compared.append(candidate.name)
    # These archives were produced locally by this run; cluster-id labels are
    # NumPy object strings, so the complete array comparison loads those labels.
    with np.load(B / "tables/animal_bootstrap_cluster_effects.npz", allow_pickle=True) as a, np.load(target / "tables/animal_bootstrap_cluster_effects.npz", allow_pickle=True) as b:
        assert a.files == b.files
        for key in a.files:
            assert np.array_equal(a[key], b[key]), key
    assert (B / "SCIENTIFIC_REPORT.md").read_text() == (target / "SCIENTIFIC_REPORT.md").read_text()
    assert (B / "figures/Figure_Multiway_Transport_caption.txt").read_text() == (target / "figures/Figure_Multiway_Transport_caption.txt").read_text()
    output = dict(status="PASS", self_contained_copy_reproduces_complete_summary=True,
                  independent_numerical_checks_pass_in_copied_folder=True,
                  all_three_reproduction_commands_pass=True,
                  identical_scientific_csv_tables=compared,
                  identical_bootstrap_arrays=True,
                  identical_report_and_caption=True,
                  source_repository_or_other_extension_needed=False,
                  temporary_copy_removed=True)
    shutil.rmtree(target)
    (B / "tables/offline_reproduction_check.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
