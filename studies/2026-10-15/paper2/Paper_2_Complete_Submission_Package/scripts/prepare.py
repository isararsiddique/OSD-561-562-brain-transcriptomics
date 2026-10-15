"""Prepare all 26 RR23 animals from official local STAR junction inputs."""
from pathlib import Path
import hashlib
import json
import re
import shutil

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def main():
    inputs = BASE / "inputs"
    provenance = []
    def freeze(path, name, role):
        destination = inputs / name
        shutil.copy2(path, destination)
        provenance.append(dict(file=name, source=str(path.relative_to(ROOT)), sha256=hashlib.sha256(destination.read_bytes()).hexdigest(), bytes=destination.stat().st_size, role=role))
        return destination

    meta_path = freeze(ROOT / "research/brain_integrative_2026/tables/refit/OSD-525_metadata.csv", "animal_metadata.csv", "ISA-derived source identity and group metadata")
    freeze(ROOT / "research/brain_integrative_2026/data/OSD-525/GLDS-525_rna_seq_qc_metrics_GLbulkRNAseq.csv", "library_qc.csv", "official library-level QC; technical repeats do not add animals")
    checksum_path = freeze(ROOT / "research/brain_integrative_2026/data/OSD-525/GLDS-525_rna_seq_processed_md5sum_GLbulkRNAseq.tsv", "publisher_md5.tsv", "official processed-file checksum sidecar")
    freeze(ROOT / "research/priority_osdr_study_2026/study_1/tables/frozen_event_definitions.csv.gz", "event_definitions.csv.gz", "historical pooled-coverage-defined 6277 local junction choices")
    prior = ROOT / "research/priority_osdr_study_2026/study_1"
    for name in ["cohort_summary.csv", "replication_summary.csv", "selected_candidate_evidence.csv", "OSD-525_frozen_event_results.csv.gz", "OSD-563_frozen_event_results.csv.gz", "OSD-352_frozen_event_results.csv.gz", "OSD-561_frozen_event_results.csv.gz"]:
        freeze(prior / "tables" / name, "prior_" + name, "previously executed results; not new replication in this extension")
    for name in ["requirements.lock.txt"]:
        shutil.copy2(ROOT / "research/paper2_revision_20261014" / name, BASE / name)
    freeze(ROOT / "research/paper2_revision_20261014/config/journal_target.json", "journal_target.json", "previously live-verified SJR category/year")
    meta = pd.read_csv(meta_path).set_index("sample", drop=False)
    assert len(meta) == 26 and meta["Source Name"].is_unique
    md5 = dict(pd.read_csv(checksum_path, sep="\t", header=None).itertuples(index=False, name=None))
    series, qc = {}, []
    files = sorted((ROOT / "research/osdr_brain_paper_2026/data/OSD-525").glob("*SJ.out.tab"))
    assert len(files) == 29
    for path in files:
        filename = path.name.split("_rna_seq_", 1)[1]
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        assert md5[filename] == digest, path
        library = filename.replace("_GLbulkRNAseq", "").removesuffix("_SJ.out.tab")
        sample = re.sub(r"_techrep[0-9]+$", "_techrep1", library)
        assert sample in meta.index
        d = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "strand", "motif", "annotated", "unique", "multi", "overhang"], dtype={"chr": str})
        good = d.strand.isin([1, 2]) & d.motif.between(1, 6) & d.annotated.eq(1) & d.overhang.ge(8) & d["unique"].gt(0) & d.chr.isin([str(x) for x in range(1, 20)] + ["X", "Y", "MT"])
        retained = d.loc[good]
        keys = retained.chr + ":" + retained.strand.astype(str) + ":" + retained.start.astype(str) + ":" + retained.end.astype(str)
        values = pd.Series(retained["unique"].to_numpy(), index=keys).groupby(level=0).sum()
        series[sample] = series[sample].add(values, fill_value=0) if sample in series else values
        qc.append(dict(library=library, sample=sample, source_animal=meta.loc[sample, "Source Name"], group=meta.loc[sample, "group"], retained_junction_rows=int(good.sum()), retained_junction_reads=int(retained["unique"].sum()), md5=digest, official_md5_match=True, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), source=str(path.relative_to(ROOT))))
    counts = pd.DataFrame(series).fillna(0).astype("int64").loc[:, meta.index]
    counts.to_csv(inputs / "all_animal_junction_counts.csv.gz", compression="gzip")
    pd.DataFrame(qc).to_csv(inputs / "junction_library_manifest.csv", index=False)
    for name in ["all_animal_junction_counts.csv.gz", "junction_library_manifest.csv"]:
        path = inputs / name
        provenance.append(dict(file=name, source="research/osdr_brain_paper_2026/data/OSD-525/", sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size, role="prepared from 29 official-MD5-matching STAR junction files; technical libraries summed into 26 animals"))
    pd.DataFrame(provenance).to_csv(inputs / "checksums.csv", index=False)
    (BASE / "config/extension_plan_sha256.txt").write_text(hashlib.sha256((BASE / "config/EXTENSION_PLAN.md").read_bytes()).hexdigest() + "\n")
    print(json.dumps(dict(animals=len(meta), groups=meta.group.value_counts().to_dict(), junction_rows=len(counts), source_libraries=len(files), publisher_md5_matches=len(qc)), indent=2), flush=True)


if __name__ == "__main__":
    main()
