"""Export only the inputs used by this study; verify every file by SHA256."""
import hashlib
import json
import shutil
from pathlib import Path

STUDY_DIRS = {"core": "analysis_20261016", "extended": "analysis_rd_20261016",
              "deep": "analysis_validation_20261017"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def export_inputs(workspace, destination):
    workspace, destination = Path(workspace).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError(f"Input destination already exists: {destination}")
    selected = {}
    def add(source, target):
        source = workspace / source
        if not source.is_file():
            raise FileNotFoundError(source)
        selected[target] = source
    for sid, glds in [(561, 556), (562, 557)]:
        for suffix in ["STAR_Unnormalized_Counts_GLbulkRNAseq.csv",
                       "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
                       "differential_expression_GLbulkRNAseq.csv",
                       "processed_md5sum_GLbulkRNAseq.tsv"]:
            name = f"GLDS-{glds}_rna_seq_{suffix}"
            add(f"data/OSD-{sid}/{name}", f"project/data/OSD-{sid}/{name}")
        name = f"OSD-{sid}_metadata_OSD-{sid}-ISA.zip"
        add(f"data/OSD-{sid}/{name}", f"project/data/OSD-{sid}/{name}")
    for code in ["CA", "DG", "FCT", "CT"]:
        name = f"GSE239336_{code}_GeneExpression_Q3norm.txt.gz"
        relative = f"OSD-682-685-698-699-brain-spatial/data/geo/{name}"
        add(relative, "project/" + relative)
    for sid, glds in [(352, 352), (563, 568), (564, 569)]:
        relative = f"research/brain_integrative_2026/data/OSD-{sid}/GLDS-{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv"
        add(relative, "project/" + relative)
    for name in ["GLDS-352_snRNA-Seq_filtered_feature_bc_matrix.h5",
                 "multiomics_brain_metadata_all_cells.csv"]:
        relative = f"research/brain_integrative_2026/data/OSD-352/{name}"
        add(relative, "project/" + relative)
    relative = "research/brain_integrative_2026/tables/rr3_metadata_download.json"
    add(relative, "project/" + relative)
    add("config/cns_targets_mapped.csv", "project/config/cns_targets_mapped.csv")
    add("CNS mRNA targets.csv", "project/CNS mRNA targets.csv")
    add("analysis_20261016/tables/verified_panel_mapping.csv", "reference/verified_panel_mapping.csv")
    # Cached primary source snapshots, not generated analysis tables or models.
    source_patterns = {
        "core": ["OSD-*_meta.json", "OSD-*_files.json", "*.download", "MGI_*.rpt", "ensembl_*.json", "ensembl_*.error.txt"],
        "extended": ["OSD-352_sample_ISA.csv", "OSD-563_sample_ISA.csv", "OSD-564_sample_ISA.csv",
                     "OSD-536_sample_ISA.csv", "OSD-352/*", "OSD-563/*", "OSD-564/*", "OSD-536/*.RCC.gz",
                     "OSD-738/*", "GSE295428_current*SOFT.txt", "GSE295428_sorted_cells/*",
                     "OSD-92[7-9]_*.json", "OSD-93[0-2]_*.json", "visium_archive_inventory.json",
                     "visium_data/*/*", "mh.all.v2025.1.Mm.symbols.gmt", "Mus_musculus.GRCm38.102.gtf.gz"],
        "deep": ["OSD-*_files.json", "OSD-*_meta.json", "*ISA.zip", "s_OSD-*.txt",
                 "*STAR_Unnormalized_Counts*.csv", "*runsheet.csv", "terminal_region_download_provenance.json"],
    }
    for stage, patterns in source_patterns.items():
        base = workspace / STUDY_DIRS[stage] / "sources"
        for pattern in patterns:
            hits = list(base.glob(pattern))
            if not hits:
                raise FileNotFoundError(f"Required source pattern has no files: {base / pattern}")
            for path in hits:
                if path.is_file():
                    selected[f"sources/{stage}/{path.relative_to(base).as_posix()}"] = path
    destination.mkdir(parents=True)
    records = []
    for target, source in sorted(selected.items()):
        out = destination / target
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, out)
        records.append({"path": target, "bytes": out.stat().st_size, "sha256": digest(out)})
    manifest = {"schema_version": 1, "snapshot_date": "2026-10-17", "files": records,
                "scope": "Primary deposited counts, metadata, annotation, supplied target panel and audited orthology reference; no fitted models or result tables."}
    (destination / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Exported {len(records)} input files ({sum(r['bytes'] for r in records):,} bytes) to {destination}", flush=True)


def verify_inputs(root):
    root = Path(root).resolve()
    manifest = json.loads((root / "input_manifest.json").read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("files"):
        raise ValueError("Unsupported or empty input manifest")
    names = set()
    for row in manifest["files"]:
        path = (root / row["path"]).resolve()
        if not path.is_relative_to(root) or row["path"] in names:
            raise ValueError(f"Invalid or duplicate manifest path: {row['path']}")
        names.add(row["path"])
        if path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
            raise ValueError(f"Input checksum mismatch: {row['path']}")
    print(f"Verified {len(names)} input SHA256 checksums", flush=True)
    return manifest
