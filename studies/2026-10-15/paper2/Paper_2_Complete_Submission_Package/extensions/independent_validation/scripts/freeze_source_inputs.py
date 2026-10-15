"""One-time source freeze and raw STAR check; not needed for portable reproduction."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

B = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.source_root.resolve()
    old = root / "research/priority_osdr_study_2026/study_1"
    out = B / "inputs"
    out.mkdir(exist_ok=True)
    files = [(old / "source_data/OSD-561_animal_junction_counts.csv.gz", "OSD-561_animal_junction_counts.csv.gz"),
             (old / "source_data/OSD-561_sample_decisions.csv", "OSD-561_sample_decisions.csv"),
             (old / "metadata/OSD-561_files.json", "OSD-561_publisher_inventory.json"),
             (old / "data/OSD-561/GLDS-556_rna_seq_processed_md5sum_GLbulkRNAseq.tsv", "GLDS-556_rna_seq_processed_md5sum_GLbulkRNAseq.tsv")]
    frozen = []
    for src, name in files:
        shutil.copy2(src, out / name)
        raw = src.read_bytes()
        frozen.append(dict(file=name, original_source_path=str(src.relative_to(root)), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == frozen[-1]["sha256"]
    matrix = pd.read_csv(out / files[0][1], index_col=0)
    md5 = pd.read_csv(out / files[3][1], sep="\t", header=None, names=["file", "md5"]).set_index("file").md5.to_dict()
    manifest = pd.read_csv(old / "tables/download_manifest.csv")
    raw_files = manifest.loc[manifest.accession.eq("OSD-561") & manifest.file.str.endswith("_SJ.out.tab")].copy()
    verified, full_prov = 0, []
    for row in raw_files.itertuples():
        path = root / row.path
        raw = path.read_bytes()
        actual_md5, actual_sha = hashlib.md5(raw).hexdigest(), hashlib.sha256(raw).hexdigest()
        published = md5.get(row.file, md5.get(row.file.removeprefix("GLDS-556_rna_seq_")))
        assert published is not None and published == actual_md5 == row.md5
        assert actual_sha == row.sha256 and len(raw) == row.bytes
        sample = row.file.split("_rna_seq_", 1)[1].removesuffix("_SJ.out.tab")
        used = sample in matrix.columns
        match = False
        if used:
            d = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "strand", "motif", "annotated", "unique", "multi", "overhang"], dtype={"chr": str})
            retain = d.strand.isin([1, 2]) & d.motif.between(1, 6) & d.annotated.eq(1) & d.overhang.ge(8) & d.unique.gt(0) & d.chr.isin([str(x) for x in range(1, 20)] + ["X", "Y", "MT"])
            keep = d.loc[retain]
            key = keep.chr + ":" + keep.strand.astype(str) + ":" + keep.start.astype(str) + ":" + keep.end.astype(str)
            counts = pd.Series(keep.unique.to_numpy(), index=key).groupby(level=0).sum()
            assert set(counts.index).issubset(matrix.index)
            np.testing.assert_array_equal(counts.reindex(matrix.index, fill_value=0).to_numpy(), matrix[sample].to_numpy())
            match = True
            verified += 1
        full_prov.append(dict(accession=row.accession, file=row.file, sample=sample, publisher_url=row.url,
                              bytes=len(raw), sha256=actual_sha, publisher_md5=published, actual_md5=actual_md5,
                              publisher_md5_match=True, included_in_collapsed_matrix=used,
                              independently_matches_collapsed_column=match))
    assert verified == matrix.shape[1] == 19
    pd.DataFrame(full_prov).to_csv(out / "RRRM2_publisher_provenance.csv", index=False)
    provenance = dict(freeze_date="2026-10-15", accession="OSD-561", source_study="RRRM2 cerebellum",
                      publisher_study_url="https://osdr.nasa.gov/bio/repo/data/studies/OSD-561",
                      publisher_inventory_url="https://osdr.nasa.gov/geode-py/ws/studies/OSD-561",
                      frozen_inputs=frozen, original_junction_files_verified=len(full_prov),
                      retained_junction_files=verified,
                      collapsed_matrix_independently_matches_all_retained_STAR_files=True,
                      collapsed_matrix_shape=list(matrix.shape),
                      filters="unique reads; annotated=1; strand1/2; canonical motif1..6; overhang>=8; unique>0; chromosomes1..19,X,Y,MT",
                      raw_junction_files_needed_for_portable_reproduction=False)
    (out / "frozen_input_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
