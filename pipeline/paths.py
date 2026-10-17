"""Shared paths for Python stages; supplied by the pipeline runner."""
import os
from pathlib import Path


def results_root():
    value = os.environ.get("OSDR_OUTPUT_ROOT")
    if not value:
        raise RuntimeError("Run this stage through run_pipeline.py.")
    return Path(value)


def stage_root(stage):
    if stage not in {"core", "extended", "deep"}:
        raise ValueError(stage)
    return results_root() / stage


def inputs_root():
    return Path(os.environ["OSDR_INPUT_ROOT"]) / "project"


def code_root():
    return Path(__file__).resolve().parent


def publication_font():
    from matplotlib.font_manager import fontManager
    available = {f.name for f in fontManager.ttflist}
    return next(f for f in ("Arial", "Liberation Sans", "DejaVu Sans") if f in available)


def threads():
    return int(os.environ.get("OSDR_THREADS", "1"))


def network_allowed():
    return os.environ.get("OSDR_ALLOW_NETWORK") == "1"


def record_path(path):
    """Store generated provenance paths relative to the output directory."""
    return Path(path).resolve().relative_to(results_root()).as_posix()


def resolve_record(value):
    path = Path(value)
    return path if path.is_absolute() else results_root() / path


def request_get(*args, **kwargs):
    if not network_allowed():
        raise RuntimeError("A cached source is missing or invalid. Re-export the input bundle, or explicitly use --allow-network.")
    import requests
    return requests.get(*args, **kwargs)


def request_post(*args, **kwargs):
    if not network_allowed():
        raise RuntimeError("A cached source is missing or invalid. Re-export the input bundle, or explicitly use --allow-network.")
    import requests
    return requests.post(*args, **kwargs)
