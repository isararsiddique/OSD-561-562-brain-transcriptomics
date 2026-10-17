"""Portable CLI, dependency checks, input integrity and sequential execution."""
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from pipeline.bundle import digest, export_inputs, verify_inputs
from pipeline.steps import STEPS, FIGURES, VERIFICATION

CODE = Path(__file__).resolve().parent
PACKAGES = ["numpy", "pandas", "scipy", "statsmodels", "scikit-learn", "pydeseq2",
            "matplotlib", "seaborn", "h5py", "requests", "Pillow", "pypdf", "pdfplumber"]


def dependencies(rscript):
    if sys.version_info < (3, 12):
        raise RuntimeError("Python 3.12 or newer is required; recorded environment uses 3.12.")
    versions, missing = {}, []
    for package in PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            missing.append(package)
    if missing:
        raise RuntimeError("Missing Python packages: " + ", ".join(missing) + ". Install pipeline/requirements.txt.")
    versions = dict(sorted((dist.metadata["Name"], dist.version) for dist in importlib.metadata.distributions() if dist.metadata.get("Name")))
    binary = shutil.which(rscript)
    if not binary:
        raise RuntimeError(f"Rscript not found: {rscript}. Install R or supply --rscript.")
    binary = str(Path(binary).resolve())
    expr = 'stopifnot(requireNamespace("limma",quietly=TRUE),requireNamespace("edgeR",quietly=TRUE),requireNamespace("statmod",quietly=TRUE));cat(paste("R",getRversion(),"limma",packageVersion("limma"),"edgeR",packageVersion("edgeR"),"statmod",packageVersion("statmod")))'
    result = subprocess.run([binary, "--vanilla", "-e", expr], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("R dependency check failed: " + result.stderr.strip())
    return {"Python": platform.python_version(), "packages": versions, "R": result.stdout.strip()}, binary


def signature(inputs, environment):
    h = hashlib.sha256()
    h.update((inputs / "input_manifest.json").read_bytes())
    h.update(json.dumps(environment, sort_keys=True).encode())
    return h.hexdigest()


def step_signatures(base):
    # A changed stage invalidates that stage and every later dependency.
    h = hashlib.sha256(base.encode())
    for name in ["paths.py", "steps.py"]:
        h.update((CODE / name).read_bytes())
    values = {}
    for step in STEPS:
        h.update(step.name.encode())
        h.update((CODE / step.script).read_bytes())
        if step.name == "integrated_figures":
            h.update((CODE / "stages/deep/workflow_figure.py").read_bytes())
        values[step.name] = h.hexdigest()
    return values


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def initialize(inputs, output):
    for group in ["core", "extended", "deep"]:
        for folder in ["sources", "tables", "prepared", "figures", "logs"]:
            (output / group / folder).mkdir(parents=True, exist_ok=True)
        source = inputs / "sources" / group
        for path in source.rglob("*"):
            if path.is_file():
                target = output / group / "sources" / path.relative_to(source)
                target.parent.mkdir(parents=True, exist_ok=True)
                # This derived download ledger is rewritten with portable paths.
                if target.exists() and path.name == "terminal_region_download_provenance.json":
                    continue
                if target.exists() and digest(target) != digest(path):
                    raise ValueError(f"Cached output source differs from the input snapshot: {target}")
                if not target.exists():
                    shutil.copy2(path, target)
    for folder in ["logs", "metadata", "deep/depth_sensitivity/tables", "deep/depth_sensitivity/prepared"]:
        (output / folder).mkdir(parents=True, exist_ok=True)


def execute(step, env, output, rscript):
    args = [a.format(deep=str(output / "deep")) for a in step.args]
    command = [rscript, "--vanilla"] if step.script.endswith(".R") else [sys.executable, "-u"]
    command += [str(CODE / step.script), *args]
    log = output / "logs" / (step.name + ".log")
    started = time.monotonic()
    print(f"RUN {step.name}  (log: {log})", flush=True)
    with log.open("w") as handle:
        subprocess.run(command, cwd=output, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    seconds = round(time.monotonic() - started, 2)
    print(f"DONE {step.name}  {seconds}s", flush=True)
    return seconds


def main(argv=None):
    if (argv or sys.argv[1:])[:1] == ['download-inputs']:
        from pipeline.acquire import main as acquire_main
        acquire_main((argv or sys.argv[1:])[1:])
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="action", required=True)
    subs.add_parser('download-inputs', help='Acquire the checksum-recorded inputs from primary sources; use download-inputs --help for options.')
    exp = subs.add_parser("export-inputs", help="Export a self-contained input bundle from this project snapshot.")
    exp.add_argument("--workspace", type=Path, default=Path.cwd())
    exp.add_argument("--destination", type=Path, required=True)
    subs.add_parser("list", help="List every stage in execution order.")
    for name in ["run", "check", "verify", "figures", "report"]:
        sub = subs.add_parser(name)
        sub.add_argument("--inputs", type=Path, default=Path("pipeline_inputs"))
        sub.add_argument("--output", type=Path, default=Path("pipeline_results"))
        sub.add_argument("--rscript", default="Rscript")
        sub.add_argument("--threads", type=int, default=4)
        if name == "run":
            sub.add_argument("--resume", action="store_true", help="Continue successfully completed stages with the same inputs, code and package versions.")
            sub.add_argument("--through", choices=[s.name for s in STEPS], help="Stop after a named stage; omitted means the complete analysis.")
            sub.add_argument("--allow-network", action="store_true", help="Explicitly permit primary repository requests. Default reproduction uses cached sources.")
    args = parser.parse_args(argv)
    try:
        if args.action == "export-inputs":
            export_inputs(args.workspace, args.destination)
            return 0
        if args.action == "list":
            for i, s in enumerate(STEPS, 1):
                print(f"{i:02} {s.name:30} {s.script} {' '.join(s.args)}")
            return 0
        if args.threads < 1:
            raise ValueError("--threads must be positive")
        inputs, output = args.inputs.resolve(), args.output.resolve()
        if output == inputs or output.is_relative_to(inputs) or inputs.is_relative_to(output):
            raise ValueError("Input and output directories must be separate, nonnested folders")
        manifest = verify_inputs(inputs)
        environment, rscript = dependencies(args.rscript)
        print(environment["R"], flush=True)
        if args.action == "check":
            print("Input and runtime checks passed. No analysis was run.")
            return 0
        if args.action != "run" and not (output / "metadata/pipeline_state.json").exists():
            raise ValueError("No pipeline run found in --output. Run the full pipeline first.")
        if args.action == "run" and not args.resume and output.exists() and any(output.iterdir()):
            raise ValueError("Output folder is not empty. Use a new folder or --resume.")
        output.mkdir(parents=True, exist_ok=True)
        lock = output / ".pipeline.lock"
        # Concurrent writers would corrupt model files and inference families.
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(fd, "w") as handle:
            handle.write(str(os.getpid()))
        try:
            current_signature = signature(inputs, environment)
            state_path = output / "metadata/pipeline_state.json"
            state = json.loads(state_path.read_text()) if state_path.exists() else {"completed": {}}
            if state.get("signature") not in {None, current_signature}:
                raise ValueError("Inputs or runtime changed. Use a new output folder for a clean rerun.")
            if args.action == "run" and args.resume and (output / "metadata/output_manifest.json").exists():
                old_inventory = json.loads((output / "metadata/output_manifest.json").read_text())
                for record in old_inventory["files"]:
                    path = output / record["path"]
                    if not path.is_file() or path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
                        raise ValueError(f"Saved output changed or disappeared: {record['path']}. Use a fresh output folder.")
                (output / "metadata/output_manifest.json").unlink()
            initialize(inputs, output)
            write_json(output / "metadata/environment.json", environment)
            write_json(output / "metadata/input_manifest.json", manifest)
            env = os.environ.copy()
            env.update(OSDR_INPUT_ROOT=str(inputs), OSDR_OUTPUT_ROOT=str(output),
                       OSDR_CORE_DIR=str(output / "core"), OSDR_EXTENDED_DIR=str(output / "extended"),
                       OSDR_DEEP_DIR=str(output / "deep"), OSDR_THREADS=str(args.threads),
                       OSDR_ALLOW_NETWORK="1" if getattr(args, "allow_network", False) else "0",
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1", PYTHONHASHSEED="0",
                       MPLCONFIGDIR=str(output / "metadata/matplotlib"),
                       PYTHONPATH=str(CODE.parent) + os.pathsep + env.get("PYTHONPATH", ""))
            state["signature"] = current_signature
            fingerprints = step_signatures(current_signature)
            selected = STEPS if args.action == "run" else VERIFICATION if args.action == "verify" else FIGURES if args.action == "figures" else [next(s for s in STEPS if s.name == "report")]
            for i, step in enumerate(selected, 1):
                print(f"[{i}/{len(selected)}]", end=" ", flush=True)
                cached = state["completed"].get(step.name, {})
                if args.action == "run" and args.resume and cached.get("fingerprint") == fingerprints[step.name]:
                    print(f"SKIP {step.name} (previously completed)", flush=True)
                else:
                    if step.name in {"bulk_deseq2", "bulk_rsem_sensitivity"}:
                        mode = "STAR" if step.name == "bulk_deseq2" else "rounded_RSEM"
                        for path in (output / "core/prepared").glob(f"OSD-*_{mode}_dds.pkl"):
                            path.unlink()
                    state["active_step"] = step.name
                    state["status"] = "running"
                    if args.action == "run":
                        position = next(j for j, s in enumerate(STEPS) if s.name == step.name)
                        for dependent in STEPS[position:]:
                            state["completed"].pop(dependent.name, None)
                    write_json(state_path, state)
                    try:
                        elapsed = execute(step, env, output, rscript)
                    except subprocess.CalledProcessError:
                        state["status"] = "failed"
                        write_json(state_path, state)
                        raise
                    state["completed"][step.name] = {"seconds": elapsed, "fingerprint": fingerprints[step.name]}
                    write_json(state_path, state)
                if step.name == getattr(args, "through", None):
                    break
            state["active_step"] = None
            state["status"] = "complete" if all(state["completed"].get(s.name, {}).get("fingerprint") == fingerprints[s.name] for s in STEPS) else "partial"
            write_json(state_path, state)
            inventory(output)
            print(f"Pipeline {state['status']}. Results: {output}", flush=True)
        finally:
            lock.unlink(missing_ok=True)
        return 0
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def inventory(output):
    records = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name not in {"output_manifest.json", ".pipeline.lock"} and "matplotlib" not in path.parts:
            records.append({"path": path.relative_to(output).as_posix(), "bytes": path.stat().st_size, "sha256": digest(path)})
    write_json(output / "metadata/output_manifest.json", {"files": records})
