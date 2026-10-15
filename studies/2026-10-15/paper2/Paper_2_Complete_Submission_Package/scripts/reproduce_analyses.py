"""Offline reproduction of every statistical analysis in the primary package (no network needed).

Runs the analysis and independent-verification steps in order. Source files are frozen; the cohort
search and download scripts are provenance records and are not rerun.
"""
from pathlib import Path
import subprocess, sys, os, time, json

B = Path(__file__).resolve().parents[1]
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
steps = []


def run(path, *args):
    start = time.time()
    print("RUN", path, flush=True)
    subprocess.run([sys.executable, str(B / path), *map(str, args)], check=True)
    steps.append({"script": str(path), "seconds": round(time.time() - start, 2)})


for s in ["analyse.py", "verify.py", "summarise.py"]:
    run("scripts/" + s)
run("extensions/gene_omnibus/scripts/reproduce.py", "--base", B)
for s in ["select_clusters.py", "analyse.py", "verify_independently.py"]:
    run("extensions/multiway_splicing/scripts/" + s, "--base", B / "inputs")
run("extensions/multiway_splicing/scripts/build_outputs.py")
for s in ["analyse_abundance.py", "verify_abundance.py"]:
    run("extensions/abundance_integration/scripts/" + s, "--base", B)
run("extensions/abundance_integration/scripts/build_figure_report.py")
for s in ["analyze_validation.py", "verify_independently.py"]:
    run("extensions/independent_validation/scripts/" + s, "--base", B)
run("extensions/independent_validation/scripts/build_report_and_figure.py")
for s in ["analyse_robustness.py", "verify_robustness.py"]:
    run("extensions/multiway_robustness/scripts/" + s, "--base", B)
run("extensions/multiway_robustness/scripts/build_report_figure.py")
for s in ["analyse.py", "verify_independently.py", "build_outputs.py"]:
    run("extensions/multiway_transport/scripts/" + s)
for s in ["analyse_context.py", "verify_context.py"]:
    run("extensions/transcript_context/scripts/" + s, "--base", B)
run("extensions/transcript_context/scripts/build_context_outputs.py")
for s in ["freeze_rr10_counts.py", "verify_and_manifest.py"]:
    run("extensions/cohort_search/scripts/" + s)
for s in ["analyse.py", "verify.py", "build_report_figure.py"]:
    run("extensions/effect_stability/scripts/" + s)
for s in ["analyse_benchmark.py", "verify_independently.py", "build_report_figure.py"]:
    run("extensions/count_model_benchmark/scripts/" + s)
for s in ["analyse.py", "verify_independently.py", "build_report_figure.py"]:
    run("extensions/confounding_sensitivity/scripts/" + s)
for s in ["analyse.py", "verify_independently.py", "build_report_figure.py"]:
    run("extensions/comparator_transport/scripts/" + s)
(B / "results" / "analysis_reproduction.json").write_text(
    json.dumps({"status": "PASS", "network_required": False, "steps": steps,
                "elapsed_seconds": sum(s["seconds"] for s in steps)}, indent=2))
print("All analyses reproduced.")
