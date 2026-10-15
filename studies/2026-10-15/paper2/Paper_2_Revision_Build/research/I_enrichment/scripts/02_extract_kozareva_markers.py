#!/usr/bin/env python3
"""Extract Kozareva et al. 2021 (Nature) Supplementary Tables 1-4 to CSV.

The published workbook is written in STRICT OOXML (purl.oclc.org namespace),
which openpyxl cannot open (it reports zero sheets), so the sheets are parsed
straight from the package XML. No marker gene is typed in by hand: everything
below comes out of the downloaded file.

Source : https://media.springernature.com/original/springer-static/esm/
         art%3A10.1038%2Fs41586-021-03220-z/MediaObjects/
         41586_2021_3220_MOESM4_ESM.xlsx
Article: Kozareva V. et al. A transcriptomic atlas of mouse cerebellar cortex
         comprehensively defines cell types. Nature 598, 214-219 (2021).
         DOI 10.1038/s41586-021-03220-z
"""
from __future__ import annotations
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
XLSX = HERE / "markers" / "Kozareva2021_SupplementaryTables_1-4.xlsx"
OUT = HERE / "markers"

NS_MAIN = "{http://purl.oclc.org/ooxml/spreadsheetml/main}"
NS_REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"
NS_ORL = "{http://purl.oclc.org/ooxml/officeDocument/relationships}"


def col_index(ref: str) -> int:
    letters = re.match(r"([A-Z]+)", ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def shared_strings(zf: zipfile.ZipFile) -> list[str]:
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall(f"{NS_MAIN}si"):
        out.append("".join(t.text or "" for t in si.iter(f"{NS_MAIN}t")))
    return out


def sheet_rows(zf: zipfile.ZipFile, path: str, sst: list[str]):
    root = ET.fromstring(zf.read(path))
    for row in root.iter(f"{NS_MAIN}row"):
        cells: dict[int, str] = {}
        for c in row.findall(f"{NS_MAIN}c"):
            ref = c.get("r") or ""
            if not ref:
                continue
            idx = col_index(ref)
            t = c.get("t")
            v = c.find(f"{NS_MAIN}v")
            if t == "s" and v is not None:
                val = sst[int(v.text)]
            elif t == "inlineStr":
                node = c.find(f"{NS_MAIN}is")
                val = "".join(x.text or "" for x in node.iter(f"{NS_MAIN}t")) if node is not None else ""
            else:
                val = v.text if v is not None else ""
            cells[idx] = val
        if cells:
            width = max(cells) + 1
            yield [cells.get(i, "") for i in range(width)]


def main() -> None:
    sha = hashlib.sha256(XLSX.read_bytes()).hexdigest()
    with zipfile.ZipFile(XLSX) as zf:
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {
            r.get("Id"): "xl/" + r.get("Target")
            for r in rels.findall(f"{NS_REL}Relationship")
        }
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        sst = shared_strings(zf)
        manifest = {"source_xlsx": str(XLSX), "sha256": sha, "sheets": {}}
        for sh in wb.find(f"{NS_MAIN}sheets").findall(f"{NS_MAIN}sheet"):
            name = sh.get("name")
            target = rid_to_target[sh.get(f"{NS_ORL}id")]
            rows = list(sheet_rows(zf, target, sst))
            dest = OUT / f"Kozareva2021_{name}.csv"
            with open(dest, "w", newline="") as fh:
                csv.writer(fh).writerows(rows)
            manifest["sheets"][name] = {
                "xml": target,
                "csv": dest.name,
                "n_rows": len(rows),
                "header": rows[0][:12] if rows else [],
            }
            print(name, "->", dest.name, len(rows), "rows; header:", rows[0][:10] if rows else [])
    (OUT / "kozareva_extraction_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
