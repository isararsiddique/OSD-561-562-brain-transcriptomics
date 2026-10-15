import json, urllib.request, urllib.parse, time
def TA(*terms):  # title-or-abstract OR group
    return "(" + " OR ".join(f'TITLE:{t} OR ABSTRACT:{t}' for t in terms) + ")"
FLIGHT = TA('"spaceflight"', '"space flight"', '"microgravity"', '"International Space Station"', '"Rodent Research"')
Q = {
 "A_splicing_x_flight": TA('"alternative splicing"','"splice junction"','"splicing"','"isoform"','"exon usage"') + " AND " + FLIGHT,
 "B_cerebellum_x_flight": TA('"cerebellum"','"cerebellar"') + " AND " + FLIGHT,
 "C_vivarium_x_flight_omics": TA('"vivarium"') + " AND " + FLIGHT + " AND " + TA('"transcriptome"','"transcriptomic"','"RNA-seq"','"gene expression"','"omics"'),
 "D_habitat_vs_vivarium_control": TA('"vivarium control"','"vivarium controls"','"basal control"') + " AND " + TA('"ground control"','"habitat control"','"habitat"'),
 "E_transport_x_flight": TA('"transportability"','"external validation"','"generalizability"','"cross-mission"','"independent mission"','"independent missions"') + " AND " + FLIGHT,
 "F_RR23_exact": '"Rodent Research-23" OR "Rodent Research 23" OR "RR-23 mission"',
 "G_RR10_exact": '"Rodent Research-10" OR "Rodent Research 10" OR "RR-10 mission"',
 "H_RRRM2_exact": '"Rodent Research Reference Mission" OR "RRRM-2" OR "RRRM2"',
 "I_discrete_pvalue_FDR": TA('"discrete p-values"','"discrete test statistics"','"discrete tests"','"minimum attainable"','"minimum achievable p"','"attainable p-value"') + " AND " + TA('"false discovery rate"','"FDR"','"multiple testing"'),
 "J_permutation_FDR_small_n_genomics": TA('"permutation"','"exact test"') + " AND " + TA('"false discovery rate"','"Benjamini"') + " AND " + TA('"small sample"','"small number of replicates"','"few replicates"','"low sample size"') + " AND " + TA('"RNA-seq"','"transcriptom*"','"gene expression"'),
 "K_split_reads_repeatability_discovery": TA('"downsampling"','"subsampling"','"read depth"','"sequencing depth"') + " AND " + TA('"differential splicing"','"splicing"') + " AND " + TA('"reproducibility"','"replicability"','"stability"'),
 "L_brain_spaceflight_mouse_RNAseq": TA('"brain"','"hippocampus"','"cortex"','"cerebellum"') + " AND " + FLIGHT + " AND " + TA('"RNA-seq"','"RNA sequencing"','"transcriptome"','"transcriptomic"') + " AND " + TA('"mice"','"mouse"'),
 "M_Tcaf1_any": TA('"Tcaf1"','"TCAF1"','"Fam115a"','"FAM115A"'),
 "N_comparator_control_choice_omics": TA('"choice of control"','"control group choice"','"comparator choice"','"choice of comparator"','"reference group"') + " AND " + TA('"transcriptom*"','"RNA-seq"','"gene expression"'),
}
out = {}
for k, q in Q.items():
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND (SRC:MED OR SRC:PMC)", "format": "json", "pageSize": 100, "resultType": "lite"})
    try:
        d = json.load(urllib.request.urlopen(url, timeout=40))
        hits = [{"title": r.get("title"), "year": r.get("pubYear"), "journal": r.get("journalTitle"), "doi": r.get("doi"), "pmid": r.get("pmid"), "authors": (r.get("authorString") or "")[:60]} for r in d["resultList"]["result"]]
        out[k] = {"query": q, "hitCount": d["hitCount"], "hits": hits}
    except Exception as e:
        out[k] = {"query": q, "error": str(e)}
    time.sleep(0.4)
out["_meta"] = {"date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "service": "Europe PMC REST, title/abstract fields, MED+PMC, first 100 by relevance"}
json.dump(out, open("epmc_novelty_search2.json", "w"), indent=1)
for k, v in out.items():
    if k.startswith("_"): continue
    print(f"## {k}: {v.get('hitCount', v.get('error'))}")
