"""Verify reproduction in a standalone beta/extensions layout, then remove test copy."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys

B = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    args = parser.parse_args()
    original = args.base.resolve()
    tmp = B / "tmp/portable_layout_check"
    beta = tmp / "beta"
    extension = beta / "extensions/independent_validation"
    if tmp.exists():
        shutil.rmtree(tmp)
    extension.mkdir(parents=True)
    for folder in ["scripts", "inputs"]:
        shutil.copytree(B / folder, extension / folder)
    shutil.copy2(B / "RETROSPECTIVE_ANALYSIS_PLAN.md", extension / "RETROSPECTIVE_ANALYSIS_PLAN.md")
    for relative in ["inputs/primary_all26_event_definitions.csv.gz", "inputs/all_animal_junction_counts.csv.gz",
                     "inputs/animal_metadata.csv", "results/all26_animal_fractions.csv.gz", "results/measurement_gates.csv.gz"]:
        path = beta / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original / relative, path)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for script in ["analyze_validation.py", "verify_independently.py"]:
        result = subprocess.run([sys.executable, "scripts/" + script, "--base", "../.."], cwd=extension,
                                env=env, capture_output=True, text=True, check=True)
        assert result.returncode == 0
    expected = json.loads((B / "tables/summary.json").read_text())
    reproduced = json.loads((extension / "tables/summary.json").read_text())
    assert expected == reproduced
    independent = json.loads((extension / "tables/independent_checks.json").read_text())
    assert independent["all_checks_passed"]
    result = dict(portable_layout_check_passed=True,
                  layout="beta/extensions/independent_validation; --base ../..",
                  external_source_repository_needed=False,
                  frozen_RRRM2_inputs_inside_extension=True,
                  required_beta_files=5,
                  complete_summary_identical=True,
                  independent_checks_passed_in_copied_layout=True,
                  temporary_copy_removed=True)
    shutil.rmtree(tmp)
    (B / "tables/portable_layout_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
