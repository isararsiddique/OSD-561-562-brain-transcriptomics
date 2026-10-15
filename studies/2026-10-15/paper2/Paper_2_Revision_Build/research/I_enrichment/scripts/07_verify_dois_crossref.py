#!/usr/bin/env python3
"""Verify every DOI cited by the enrichment layer against Crossref."""
from __future__ import annotations
import json
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
RAW = HERE / "raw"

DOIS = {
    "kozareva2021": "10.1038/s41586-021-03220-z",
    "kolberg2023gprofiler": "10.1093/nar/gkad347",
    "raudvere2019gprofiler": "10.1093/nar/gkz369",
    "reimand2007gscs": "10.1093/nar/gkm226",
    "benjamini1995": "10.1111/j.2517-6161.1995.tb02031.x",
    "benjamini2001": "10.1214/aos/1013699998",
    "ashburner2000go": "10.1038/75556",
    "kanehisa2000kegg": "10.1093/nar/28.1.27",
    "milacic2024reactome": "10.1093/nar/gkad1025",
}
UA = "Paper2-enrichment-doi-check/1.0 (mailto:noreply@example.org)"


def fetch(doi: str) -> dict:
    req = urllib.request.Request(
        f"https://api.crossref.org/works/{doi}", headers={"User-Agent": UA}
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())["message"]


def main() -> None:
    out = {}
    for key, doi in DOIS.items():
        try:
            m = fetch(doi)
            date = m.get("issued", {}).get("date-parts", [[None]])[0]
            out[key] = {
                "doi": doi,
                "verified": True,
                "title": (m.get("title") or [""])[0],
                "container": (m.get("container-title") or [""])[0],
                "year": date[0] if date else None,
                "volume": m.get("volume"),
                "page": m.get("page"),
                "type": m.get("type"),
                "first_authors": [
                    f"{a.get('family','')} {a.get('given','')}".strip()
                    for a in (m.get("author") or [])[:3]
                ],
            }
            print(f"OK  {doi}  {out[key]['year']}  {out[key]['container']}  "
                  f"{out[key]['volume']}:{out[key]['page']}  {out[key]['title'][:70]}")
        except Exception as e:  # noqa: BLE001
            out[key] = {"doi": doi, "verified": False, "error": repr(e)}
            print(f"FAIL {doi}  {e}")
        time.sleep(1)
    (RAW / "crossref_doi_verification.json").write_text(json.dumps(out, indent=2))
    print("\nverified:", sum(v.get("verified", False) for v in out.values()), "/", len(out))


if __name__ == "__main__":
    main()
