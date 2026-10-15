import json, urllib.request, urllib.parse, time, sys
Q = {
 "Q1_cerebellum_spaceflight_transcriptome": '(cerebellum OR cerebellar) AND (spaceflight OR "space flight" OR microgravity) AND (RNA-seq OR transcriptom* OR "gene expression")',
 "Q2_splicing_spaceflight_any": '("alternative splicing" OR "splice junction" OR "splicing") AND (spaceflight OR "space flight" OR "International Space Station")',
 "Q3_splicing_brain_spaceflight": '("alternative splicing" OR splicing OR "splice junction") AND (brain OR cerebell* OR hippocamp* OR neuron*) AND (spaceflight OR microgravity OR "hindlimb unloading")',
 "Q4_vivarium_vs_habitat_control": '("vivarium control" OR "vivarium controls") AND ("ground control" OR habitat) AND (spaceflight OR "Rodent Research")',
 "Q5_RR23_OSD525": '("Rodent Research-23" OR "RR-23" OR "OSD-525" OR "GLDS-525")',
 "Q6_RRRM2_OSD561": '("RRRM-2" OR "RRRM2" OR "Rodent Research Reference Mission-2" OR "OSD-561" OR "GLDS-561")',
 "Q7_RR10_OSD563": '("Rodent Research-10" OR "RR-10" OR "OSD-563" OR "GLDS-563")',
 "Q8_transportability_omics": '(transportability OR "external validity" OR "transport") AND (spaceflight OR "space omics" OR GeneLab)',
 "Q9_permutation_minimum_p_FDR_attainable": '("minimum achievable p" OR "minimum attainable p" OR "smallest possible p" OR "attainable p-value") AND ("false discovery rate" OR "multiple testing")',
 "Q10_permutation_small_sample_FDR_impossible": '(permutation OR "exact test") AND ("small sample" OR "few replicates") AND ("false discovery rate" OR Benjamini) AND (granularity OR discreteness OR "discrete p-values")',
 "Q11_cross_mission_replication_mouse": '(spaceflight OR "Rodent Research") AND (mouse OR mice) AND ("cross-mission" OR "multiple missions" OR "independent missions") AND (replicat* OR reproducib* OR consistency)',
 "Q12_GeneLab_reanalysis_confound_dissection_date": '(spaceflight OR GeneLab) AND ("dissection" OR "euthanasia" OR "collection date" OR "batch") AND (confound*) AND (transcriptom* OR RNA-seq)',
 "Q13_junction_counts_STAR_SJ_reuse_spaceflight": '("SJ.out.tab" OR "junction counts" OR "junction reads") AND (spaceflight OR GeneLab OR OSDR)',
 "Q14_tcaf1": '(Tcaf1 OR "TRPM8 channel-associated factor 1" OR FAM115A)',
 "Q15_female_mice_brain_spaceflight": '(female) AND (mice OR mouse) AND (brain OR cerebell*) AND (spaceflight OR "International Space Station") AND (transcriptom* OR RNA-seq)',
}
out = {}
for k, q in Q.items():
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 25, "resultType": "lite", "sort": "CITED desc"})
    try:
        d = json.load(urllib.request.urlopen(url, timeout=30))
        hits = [{"title": r.get("title"), "year": r.get("pubYear"), "journal": r.get("journalTitle"), "doi": r.get("doi"), "pmid": r.get("pmid"), "src": r.get("source")} for r in d["resultList"]["result"]]
        out[k] = {"query": q, "hitCount": d["hitCount"], "top": hits}
    except Exception as e:
        out[k] = {"query": q, "error": str(e)}
    time.sleep(0.4)
out["_meta"] = {"date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "service": "Europe PMC REST search, lite, sorted by citations, top 25"}
json.dump(out, open("epmc_novelty_search.json", "w"), indent=1)
for k, v in out.items():
    if k.startswith("_"): continue
    print(f"\n## {k}: {v.get('hitCount', v.get('error'))}")
    for h in v.get("top", [])[:12]:
        print(f"  {h['year']} | {(h['journal'] or '')[:28]} | {(h['title'] or '')[:120]}")
