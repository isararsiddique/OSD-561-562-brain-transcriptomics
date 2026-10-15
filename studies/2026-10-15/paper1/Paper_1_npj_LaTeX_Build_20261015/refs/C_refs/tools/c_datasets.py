#!/usr/bin/env python3
"""Worker C: verify the 11 NASA OSDR dataset entries (manuscript refs 23-33) and write ../C_datasets.bib.

For each OSD accession:
  * live OSDR metadata API  https://osdr.nasa.gov/osdr/data/osd/meta/<n>  (title, 'DOI' study comment)
  * deposited metadata snapshot in the source package inputs/osd_<n>_public_metadata.json (same fields)
  * DataCite REST API       https://api.datacite.org/dois/<doi>          (title, publisher, year, URL)
  * doi.org handle API      https://doi.org/api/handles/<doi>            (DOI resolves)
  * DataCite search for accessions without a DOI in OSDR metadata (query OSD-<n>)
"""
import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from c_verify import get_json, CACHE, OUT, old_reference_texts, tex  # noqa: E402

PKG = pathlib.Path(""
                   "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/inputs")
ACCESSED = "15 October 2026"
OSD = {23: 202, 24: 352, 25: 374, 26: 611, 27: 525, 28: 561, 29: 562, 30: 563, 31: 564, 32: 612, 33: 613}
KEYS = {202: "OSDR_OSD202", 352: "OSDR_OSD352", 374: "OSDR_OSD374", 611: "OSDR_OSD611", 525: "OSDR_OSD525",
        561: "OSDR_OSD561", 562: "OSDR_OSD562", 563: "OSDR_OSD563", 564: "OSDR_OSD564", 612: "OSDR_OSD612",
        613: "OSDR_OSD613"}


def study_fields(js):
    """Pull title and 'DOI' comment from an OSDR meta JSON (live or snapshot)."""
    if not js:
        return None, None, None
    st = js.get("study") or {}
    sid = next(iter(st), None)
    s = st.get(sid) or {}
    studies = s.get("studies") or []
    title, doi, rel = None, None, None
    if studies:
        title = studies[0].get("title")
        for c in studies[0].get("comments", []):
            if c.get("name") == "DOI" and c.get("value"):
                doi = c["value"].strip()
        rel = s.get("releaseDate") or studies[0].get("publicReleaseDate")
    if not doi:
        m = re.search(r'"name":\s*"DOI",\s*"value":\s*"([^"]+)"', json.dumps(js))
        doi = m.group(1) if m else None
    return title, doi, rel


def main():
    old = old_reference_texts()
    rows, bib = [], []
    for no, n in OSD.items():
        live = get_json("osdr_meta", f"OSD-{n}", f"https://osdr.nasa.gov/osdr/data/osd/meta/{n}", spacing=1.0)
        lt, ldoi, lrel = study_fields(live.get("json"))
        snap = PKG / f"osd_{n}_public_metadata.json"
        st, sdoi, _ = study_fields(json.loads(snap.read_text())) if snap.exists() else (None, None, None)
        doi = ldoi or sdoi
        doi_src = "OSDR live metadata" if ldoi else ("deposited snapshot" if sdoi else None)
        dc, resolves, search_hits = None, None, []
        if not doi:
            # OSDR study comment is empty: look the dataset up in DataCite (NASA GeneLab prefix 10.26030) by its exact
            # title, and accept a DOI only if DataCite's landing URL is this accession's OSDR page.
            sr = get_json("datacite_search", f"title_OSD-{n}", "https://api.datacite.org/dois",
                          params={"prefix": "10.26030", "query": f'titles.title:"{lt or st}"', "page[size]": 25}, spacing=1.0)
            for d in (sr.get("json") or {}).get("data", []):
                a = d.get("attributes") or {}
                search_hits.append({"doi": a.get("doi"), "title": ((a.get("titles") or [{}])[0]).get("title"), "url": a.get("url")})
            match = [h for h in search_hits if (h["url"] or "").rstrip("/").endswith(f"/OSD-{n}")
                     and re.sub(r"\W", "", (h["title"] or "").lower()) == re.sub(r"\W", "", (lt or st or "").lower())]
            if len(match) == 1:
                doi, doi_src = match[0]["doi"], "DataCite (prefix 10.26030, exact title, landing URL = this OSD page); absent from OSDR study metadata"
        if doi:
            dcr = get_json("datacite", doi, f"https://api.datacite.org/dois/{doi}", spacing=1.0)
            dc = ((dcr.get("json") or {}).get("data") or {}).get("attributes")
            h = get_json("doi_handle", doi, f"https://doi.org/api/handles/{doi}", spacing=0.5)
            resolves = bool(h.get("json") and h["json"].get("responseCode") == 1)

        title = lt or st or ((dc or {}).get("titles") or [{}])[0].get("title")
        ms = old.get(no, "")
        ms_title = re.sub(r"^NASA Open Science Data Repository\. ", "", ms).split(". OSD-")[0]
        pub = (dc or {}).get("publisher") if dc else None
        if isinstance(pub, dict):
            pub = pub.get("name")
        row = {
            "manuscript_ref_no": no, "key": KEYS[n], "accession": f"OSD-{n}",
            "url": f"https://osdr.nasa.gov/bio/repo/data/studies/OSD-{n}",
            "title_osdr_live": lt, "title_snapshot": st, "title_in_manuscript": ms_title,
            "title_matches_manuscript": bool(title) and re.sub(r"\W", "", (title or "").lower()) == re.sub(r"\W", "", ms_title.lower()),
            "doi": doi, "doi_source": doi_src, "doi_in_manuscript": bool(doi) and doi.lower() in ms.lower(),
            "doi_resolves": resolves,
            "datacite_title": ((dc or {}).get("titles") or [{}])[0].get("title") if dc else None,
            "datacite_publisher": pub,
            "datacite_year": (dc or {}).get("publicationYear") if dc else None,
            "datacite_url": (dc or {}).get("url") if dc else None,
            "datacite_url_matches_accession": bool(dc) and ((dc or {}).get("url") or "").rstrip("/").endswith(f"/OSD-{n}"),
            "datacite_search_hits_if_no_doi": search_hits,
            "osdr_live_status": live.get("status"), "release_date": lrel,
            "manuscript_string": ms,
        }
        rows.append(row)
        fields = [("author", "{{NASA Open Science Data Repository}}"), ("title", "{{" + tex(title or ms_title) + "}}"),
                  ("howpublished", f"NASA Open Science Data Repository (OSDR), accession OSD-{n}"),
                  ("publisher", tex(pub or "")),
                  ("year", str(row["datacite_year"] or "")), ("url", row["url"])]
        if doi and resolves and row["datacite_url_matches_accession"]:
            fields.append(("doi", doi))
        note = f"Accessed {ACCESSED}" + ("" if doi else "; no dataset DOI listed in OSDR metadata")
        fields.append(("note", note))
        body = ",\n".join(f"  {k} = {{{v}}}" if not v.startswith("{") else f"  {k} = {v}" for k, v in fields if v)
        bib.append(f"@misc{{{KEYS[n]},\n{body}\n}}\n")
    (OUT / "C_datasets.bib").write_text("% NASA OSDR dataset entries for Paper 1 (manuscript refs 23-33), Worker C, 15 October 2026.\n"
                                        "% Titles from the live OSDR metadata API; DOIs from OSDR study comments, checked in DataCite and doi.org.\n\n"
                                        + "\n".join(bib))
    (OUT / "cache" / "datasets_audit.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False))
    for r in rows:
        print(r["accession"], "| doi", r["doi"], r["doi_source"], "resolves", r["doi_resolves"], "| DC", r["datacite_publisher"],
              r["datacite_year"], "| url=acc", r["datacite_url_matches_accession"], "| in MS", r["doi_in_manuscript"],
              "| title match", r["title_matches_manuscript"], "| live", r["osdr_live_status"])
        if not r["title_matches_manuscript"]:
            print("    OSDR :", r["title_osdr_live"] or r["title_snapshot"])
            print("    MS   :", r["title_in_manuscript"])



if __name__ == "__main__":
    main()
