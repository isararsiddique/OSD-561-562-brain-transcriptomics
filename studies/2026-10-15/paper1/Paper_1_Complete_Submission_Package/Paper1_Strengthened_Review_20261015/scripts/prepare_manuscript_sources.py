#!/usr/bin/env python3
"""Retrieve exact reference metadata and extract ethics text from source snapshots."""
import concurrent.futures
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
DOIS = [
    "10.1016/j.cell.2020.10.050", "10.1093/nar/gkae1116", "10.1016/j.isci.2021.102361",
    "10.1038/nrg2825", "10.1038/s41593-022-01022-8", "10.1016/j.neuron.2018.04.001",
    "10.1074/jbc.M702044200", "10.1186/s13059-014-0550-8", "10.1093/bioinformatics/btad547",
    "10.1111/j.2517-6161.1995.tb02031.x", "10.1371/journal.pgen.0030161",
    "10.1038/nbt.2931", "10.1002/sim.8086", "10.1016/j.crmeth.2022.100325",
    "10.1038/s41526-022-00217-4", "10.1096/fj.202001754r",
    "10.1038/s41467-024-48916-8", "10.1093/biostatistics/kxv027",
    "10.1093/biostatistics/kxr034", "10.1016/j.neuroimage.2014.01.060",
    "10.1038/s41467-026-68737-1", "10.3390/cells14060417",
]


def walk(item):
    if isinstance(item, dict):
        yield item
        for value in item.values():
            yield from walk(value)
    elif isinstance(item, list):
        for value in item:
            yield from walk(value)


def reference(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    cache = HERE / "inputs/reference_metadata"
    cache.mkdir(exist_ok=True)
    destination = cache / (doi.replace("/", "_") + ".json")
    record = None
    candidates = [destination]
    candidates += list((HERE / "extensions/literature_novelty/cache").glob("crossref*.json"))
    # Existing verified bibliographic records may be reused; DOI identity is checked.
    candidates += list((HERE.parents[0] / "paper2_beta_20261014/config").glob("*references_verified.json"))
    for path in candidates:
        if not path.exists():
            continue
        obj = json.loads(path.read_text())
        entries = obj if isinstance(obj, list) else [obj.get("message", obj)]
        for item in entries:
            if str(item.get("DOI", item.get("doi", ""))).lower() == doi.lower():
                record = item
                break
        if record is not None:
            break
    if record is None:
        for attempt in range(4):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "OSDR-reference-review/1.0 (mailto:isarar@convolity.com)"})
                with urllib.request.urlopen(request, timeout=30) as response:
                    record = json.load(response)["message"]
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt+1))
    destination.write_text(json.dumps(record, indent=2) + "\n")
    if str(record.get("DOI", record.get("doi", ""))).lower() != doi.lower():
        raise ValueError(f"DOI mismatch: {doi}")
    date = record.get("published-print", record.get("published", record.get("published-online")))
    title = record["title"][0] if isinstance(record["title"], list) else record["title"]
    journal = record.get("container-title", [record.get("journal", "")])[0]
    return dict(doi=doi, title=title, authors=record.get("author", record.get("authors", [])),
                journal=journal, volume=record.get("volume", ""),
                pages=record.get("page", record.get("article-number", record.get("pages", ""))),
                year=date["date-parts"][0][0] if date else record["year"],
                source=url, verified_utc=datetime.now(timezone.utc).isoformat())


def main():
    output = HERE / "results/manuscript_sources.json"
    records = []
    for doi in DOIS:
        records.append(reference(doi))
        print("Verified", doi, flush=True)
    datasets = []
    for n in ["525", "561", "562", "563", "564", "612", "613"]:
        path = HERE / f"inputs/osd_{n}_public_metadata.json"
        snapshot = json.loads(path.read_text())
        objects = list(walk(snapshot))
        title = next(str(obj["title"]) for obj in objects if "title" in obj)
        ethics = [dict(protocol=obj.get("name", ""), statement=obj["description"])
                  for obj in objects if "IACUC" in str(obj.get("description", ""))]
        dois = sorted(set(re.findall(r'10\.\d{4,9}/[^\s"\\<>]+', path.read_text())))
        datasets.append(dict(accession=f"OSD-{n}", title=title, dois=dois,
                             url=f"https://osdr.nasa.gov/bio/repo/data/studies/OSD-{n}",
                             metadata_url=f"https://osdr.nasa.gov/osdr/data/osd/meta/{n}",
                             ethics=ethics, local_snapshot=str(path.relative_to(HERE))))
    output.write_text(json.dumps(dict(references=records, datasets=datasets), indent=2) + "\n")
    print(f"Verified {len(records)} reference DOIs; extracted {len(datasets)} dataset records.")


if __name__ == "__main__":
    main()
