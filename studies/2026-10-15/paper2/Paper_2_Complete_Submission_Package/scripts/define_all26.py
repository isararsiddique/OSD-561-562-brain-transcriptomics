"""Create group-label-invariant features and exact GTF intron annotation."""
from pathlib import Path
from collections import defaultdict
import gzip
import hashlib
import re

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def main():
    counts = pd.read_csv(BASE / "inputs/all_animal_junction_counts.csv.gz", index_col=0)
    assert counts.shape[1] == 26
    selected = counts.loc[(counts >= 5).sum(axis=1) >= 3]
    parts = selected.index.str.split(":", expand=True)
    info = pd.DataFrame(dict(chr=parts.get_level_values(0), strand=parts.get_level_values(1), start=parts.get_level_values(2), end=parts.get_level_values(3)), index=selected.index)
    rows = []
    for boundary in ["start", "end"]:
        for key, ids in info.groupby(["chr", "strand", boundary], sort=True).groups.items():
            if len(ids) != 2:
                continue
            first, second = sorted(ids)
            x = selected.loc[[first, second]].to_numpy()
            total = x.sum(axis=0)
            if (total < 20).any():
                continue
            frac = x[0] / total
            if not .05 < frac.mean() < .95:
                continue
            chromosome, strand, position = key
            event_id = "|".join([chromosome, strand, position, boundary, first, second])
            rows.append(dict(event_id=event_id, chr=chromosome, strand=int(strand), boundary=int(position), boundary_type=boundary, junction_a=first, junction_b=second, pooled_mean_fraction=float(frac.mean()), minimum_pair_reads=float(total.min())))
    definitions = pd.DataFrame(rows).sort_values("event_id").reset_index(drop=True)
    assert definitions.event_id.is_unique
    needed = set(definitions.junction_a) | set(definitions.junction_b)
    exons = defaultdict(list)
    transcript_info, symbols = {}, {}
    gtf = ROOT / "research/osdr_brain_paper_2026/data/Mus_musculus.GRCm39.112.gtf.gz"
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            gene = re.search(r'gene_id "([^"]+)"', parts[8])[1]
            transcript = re.search(r'transcript_id "([^"]+)"', parts[8])[1]
            symbol = re.search(r'gene_name "([^"]+)"', parts[8])
            if symbol:
                symbols[gene] = symbol[1]
            transcript_info[transcript] = (parts[0], 1 if parts[6] == "+" else 2, gene)
            exons[transcript].append((int(parts[3]), int(parts[4])))
    junction_genes = defaultdict(set)
    for transcript, coordinates in exons.items():
        chromosome, strand, gene = transcript_info[transcript]
        ordered = sorted(set(coordinates))
        for first, second in zip(ordered, ordered[1:]):
            key = f"{chromosome}:{strand}:{first[1]+1}:{second[0]-1}"
            if key in needed:
                junction_genes[key].add(gene)
    gids, names, mapping_counts = [], [], []
    for row in definitions.itertuples():
        genes = junction_genes[row.junction_a] & junction_genes[row.junction_b]
        mapping_counts.append(len(genes))
        gene = next(iter(genes)) if len(genes) == 1 else ""
        gids.append(gene)
        names.append(symbols.get(gene, ""))
    definitions["gene_id"] = gids
    definitions["gene_symbol"] = names
    definitions["gene_mapping_count"] = mapping_counts
    name = "primary_all26_event_definitions.csv.gz"
    definitions.to_csv(BASE / "inputs" / name, index=False, compression="gzip")
    provenance = pd.read_csv(BASE / "inputs/checksums.csv")
    provenance = provenance.loc[provenance.file != name]
    file = BASE / "inputs" / name
    row = dict(file=name, source="all 26 pooled animal counts; GRCm39.112 transcript introns", sha256=hashlib.sha256(file.read_bytes()).hexdigest(), bytes=file.stat().st_size, role="group-label-invariant prepared junction definitions; raw GTF SHA256 " + hashlib.sha256(gtf.read_bytes()).hexdigest())
    pd.concat([provenance, pd.DataFrame([row])], ignore_index=True).to_csv(BASE / "inputs/checksums.csv", index=False)
    print(dict(definitions=len(definitions), uniquely_mapped=int((definitions.gene_mapping_count == 1).sum()), measured_animals=26, feature_selection_uses_group_labels=False), flush=True)


if __name__ == "__main__":
    main()
