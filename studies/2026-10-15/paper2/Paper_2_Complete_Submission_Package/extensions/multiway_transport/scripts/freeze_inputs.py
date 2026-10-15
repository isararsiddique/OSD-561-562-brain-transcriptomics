"""One-time self-contained frozen-input snapshot; analysis has no repository dependency."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import pandas as pd

B = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="Verified npj package root, used only for initial freezing")
    root = parser.parse_args().base.resolve()
    files = {
        "inputs/all_animal_junction_counts.csv.gz": "RR23_animal_junction_counts.csv.gz",
        "inputs/animal_metadata.csv": "RR23_animal_metadata.csv",
        "extensions/multiway_splicing/tables/selected_cluster_definitions.csv": "RR23_frozen_cluster_definitions.csv",
        "extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz": "RR23_frozen_alternative_compositions.csv.gz",
        "extensions/multiway_splicing/tables/frozen_features.npz": "RR23_frozen_features.npz",
        "extensions/multiway_splicing/tables/selection_summary.json": "RR23_original_selection_summary.json",
        "extensions/multiway_splicing/annotation/junction_gene_lookup.csv.gz": "junction_gene_lookup.csv.gz",
        "extensions/multiway_splicing/annotation/source.json": "annotation_source.json",
        "extensions/independent_validation/inputs/OSD-561_animal_junction_counts.csv.gz": "RRRM2_animal_junction_counts.csv.gz",
        "extensions/independent_validation/inputs/OSD-561_sample_decisions.csv": "RRRM2_animal_metadata.csv",
        "extensions/independent_validation/inputs/frozen_input_provenance.json": "RRRM2_original_input_provenance.json",
        "extensions/independent_validation/inputs/RRRM2_publisher_provenance.csv": "RRRM2_publisher_provenance.csv",
        "extensions/independent_validation/inputs/OSD-561_publisher_inventory.json": "OSD-561_publisher_inventory.json",
        "extensions/independent_validation/inputs/GLDS-556_rna_seq_processed_md5sum_GLbulkRNAseq.tsv": "OSD-561_publisher_MD5.tsv",
        "extensions/independent_validation/tables/summary.json": "prior_binary_transport_summary.json",
    }
    rows = []
    for relative, name in files.items():
        source, target = root / relative, B / "inputs" / name
        target.parent.mkdir(exist_ok=True)
        original_hash = digest(source)
        shutil.copy2(source, target)
        assert digest(target) == original_hash
        rows.append(dict(file="inputs/" + name, original_package_relative_path=relative,
                         bytes=target.stat().st_size, sha256=original_hash, original_source_unchanged=True))
    provenance = json.loads((B / "inputs/RRRM2_original_input_provenance.json").read_text())
    assert provenance["collapsed_matrix_independently_matches_all_retained_STAR_files"]
    assert provenance["original_junction_files_verified"] == 26 and provenance["retained_junction_files"] == 19
    lookup = json.loads((B / "inputs/annotation_source.json").read_text())
    assert lookup["counts_sha256"] == digest(B / "inputs/RR23_animal_junction_counts.csv.gz")
    pd.DataFrame(rows).to_csv(B / "inputs/frozen_input_manifest.csv", index=False)
    out = dict(freeze_date="2026-10-15", complete_self_contained_analysis_inputs=True,
               original_sources_unchanged=True, frozen_files=len(rows),
               RRRM2_original_junction_files_publisher_MD5_verified=26,
               RRRM2_retained_STAR_columns_independently_reconstructed=19,
               frozen_cluster_family=923, sources_recorded_as_package_relative_paths=True)
    (B / "inputs/freeze_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
