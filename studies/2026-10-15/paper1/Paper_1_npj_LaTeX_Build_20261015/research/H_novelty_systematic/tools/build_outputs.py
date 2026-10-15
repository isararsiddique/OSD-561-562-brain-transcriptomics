"""Write H_screening.csv and H_search_log.json from the cached stages.

Decisions: records failing the registered 3a rule are EXCLUDE at stage 3a (the 120-record audit sample is marked
stage 3c). Records passing 3a were title-screened by the screener (stage 3b); abstracts were read for every title
judged possibly relevant. Everything not listed in FLAG below was excluded at 3b as off-topic for the claim it was
retrieved for. FLAG holds the screener's non-exclude decisions; stage '4' means the full text was read.
"""
import csv, json, os
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = os.path.join(OUT, 'cache')

R, N, A, U4 = 'RELEVANT, NOT ANTICIPATING', 'NARROWS', 'ANTICIPATES', 'UNRESOLVED (possible NARROWS/ANTICIPATES)'
FLAG = {
 # ---- C1 spaceflight
 'H00066': ('C1', '4', R, 'Lai Polo 2020: collection technique confounds liver/muscle RNA-seq; no brain, no IEG panel. Known prior art, already cited.'),
 'H00069': ('C1', '4', R, 'Overbey 2022: describes RRRM-2 ISS-terminal vs LAR design; scRNA of LAR non-brain tissues only; no protocol-effect size in brain.'),
 'H00712': ('C1', '3b', R, 'Choi 2020: validates on-orbit euthanasia/dissection hardware (RR-1); tissue quality, not brain IEG protocol effect.'),
 'H00129': ('C1', '4', R, 'Grandke 2026 (RRRM-1, BALB/c, small RNA, 13 tissues): TERM vs LAR contrast, both arms ketamine/xylazine; extraction explains 0.8% of pooled variance; no brain-specific size, no mRNA/IEG.'),
 'H02521': ('C1', '3b', R, 'Grandke 2025 preprint version of H00129.'),
 'H00135': ('C1', '4', R, 'Kremsky 2023 (RR-18, GeoMx): Arc is a top spaceflight DEG in cortex; flight effect only, no protocol contrast. Supports the paper\'s statement that IEGs appear among flight genes.'),
 'H00144': ('C1', '3b', R, 'Holley 2022 (35-day LAR, NanoString): flight effects only.'),
 'H00137': ('C1', '4', R, 'Masarapu 2024: RR-3 brains, ISS-terminal flight vs ground; no protocol contrast.'),
 'H00148': ('C1', '3b', R, 'Masarapu 2023 preprint of H00137.'),
 'H02189': ('C1', '4', R, 'Earth reference controls in spaceflight omics (2020): ground-control design; no brain, no IEG.'),
 'H04084': ('C1', '3b', R, 'SSRN version of H02189.'),
 'H00102': ('C1', '3b', R, 'GeneLab consistent sample/data processing pipeline abstract; names sample processing as a confounder; no data on brain.'),
 'H00104': ('C1', '3b', R, 'Tissue preservation for flight experiments (2016); RNA quality, no IEG.'),
 'H02470': ('C1', '3b', R, 'RR-1 biospecimen retrieval (2016); logistics.'),
 'H00180': ('C1', '3b', R, 'Neurolab (2007): Fos/FRA protein in rat forebrain during flight vs return; flight timing effect, not collection protocol.'),
 'H00183': ('C1', '3b', R, 'Neurolab Fos/FRA in precerebellar structures (2003); flight effect.'),
 'H02671': ('C1', '3b', R, 'Rat brain gene expression under space-flight conditions (2003); no abstract; title-level flight effect.'),
 'H04051': ('C3', '4', R, 'GeneLab batch-correction study (2023): seven liver datasets, Animal Return (LAR vs ISS-T) and preservation tabulated as batch variables; no held-out ground animals, no gain bound.'),
 'H02149': ('C3', '3b', R, 'Multi-mission liver harmonisation for ML (npj Microgravity 2024); no ground-only holdout design.'),
 'H02162': ('C3', '3b', R, 'Preprint of H02149.'),
 'H00058': ('C4', '4', R, 'Ganti 2026 OSDR OXPHOS census (Research Square): liver/muscle/retina only; none of the seven deposits analysed.'),
 'H02456': ('C4', '4', U4, 'Piekarski et al. 2026, Life Sci Space Res: single-nucleus brain transcriptomics of spaceflight mice with age effects; CLOSED access, no abstract in any source; dataset not identifiable. Could use OSD-612/613 snRNA component or RR-3 (Masarapu) data.'),
 'H02669': ('C4', '3b', R, 'Roy 2026 J Neuroinflammation: RR-1 brains.'),
 'H02393': ('C4', '3b', R, 'Cosmic kidney disease (2024): kidney; RR-10 kidney only.'),
 'H02101': ('C4', '3b', R, 'RRRM-2 kidney (OSD-771) preprint; different deposit and tissue.'),
 'H02681': ('C4', '3b', R, 'GeneLab dataset record: scRNA cerebellum, 3 and 8 month mice (RRRM-1 design of Grandke); not one of the seven deposits.'),
 'H02706': ('C4', '3b', R, 'KTH 2021 thesis: spatial transcriptomics of spaceflight mouse hippocampus; predates RRRM-2 flight (Dec 2021); RR-3 programme.'),
 'H02632': ('C4', '3b', R, 'GeneLab RR-10 cecum dataset record; not brain.'),
 'H02633': ('C4', '3b', R, 'GeneLab RR-10 lung dataset record; not brain.'),
 # ---- C1 non-spaceflight closest
 'H00412': ('C1', '3b', R, 'Staib-Lasarzik 2014: anaesthesia for euthanasia changes brain mRNA in mice (non-spaceflight). Already cited.'),
 'H02726': ('C1', '3b', R, 'Isoflurane alters IEG expression (Front Behav Neurosci 2016), non-spaceflight.'),
 'H02891': ('C1', '3b', R, 'Pentobarbital/isoflurane/propofol on IEG expression (2000), non-spaceflight.'),
 'H02892': ('C1', '3b', R, 'Pentobarbital/halothane on c-fos and junB (1992), non-spaceflight.'),
 'H04652': ('C1', '3b', R, 'Anaesthesia/euthanasia method for brain collection alters MAPK signalling (2019), non-spaceflight.'),
 'H00511': ('C1', '3b', R, 'CO2 vs focused-beam microwave euthanasia: hypothalamic mRNA (2011), non-spaceflight.'),
 'H00600': ('C1', '3b', R, 'Decapitation vs microwave: NPY mRNA (2001), non-spaceflight.'),
 'H00642': ('C1', '3b', R, 'mRNA integrity after microwave killing (1999), non-spaceflight.'),
 'H04180': ('C1', '3b', R, 'Ketamine IEG brain mapping (2022), non-spaceflight.'),
 'H04741': ('C1', '3b', R, 'Version of H04180.'),
 'H02822': ('C1', '3b', R, 'Hippocampal gene expression 2 days after general anaesthesia (2006), non-spaceflight.'),
 'H00372': ('C1', '3b', R, 'Post-mortem time-dependent activity/cell-specific gene changes (2021), non-spaceflight.'),
 'H00251': ('C1', '3b', R, 'Post-mortem interval effects on snRNA-seq (2025), non-spaceflight.'),
 'H04520': ('C1', '3b', R, 'Euthanasia, dissection and post-mortem delay on retina metabolome (2018), non-spaceflight.'),
 # ---- C2
 'H04759': ('C2', '4', R, 'Two-step batch correction (bioRxiv 2021/Biostatistics): derives the correlation operator (I-H)(I-H)^T induced by correction and corrects inference; orthogonal-type operator, no gain >1, no spectral-norm diagnostic. Closest operator-level prior art.'),
 'H05265': ('C2', '4', R, 'Batch Probing Score (bioRxiv 2025): batch stays classifier-detectable AFTER correction. Detectability used as a failure metric, not separated from loading precision or operator gain.'),
 'H00717': ('C2', '4', R, 'RiTex 2026 (radiomics harmonisation): theoretical stability bound uses the operator norm of the harmonisation map derivative; not omics, not proposed as a routine reported diagnostic.'),
 'H00718': ('C2', '4', R, 'Spectral deconfounding for high-dimensional regression (2025): spectral transform shrinks singular values; regression, not correction-operator gain reporting.'),
 'H03042': ('C2', '3b', R, 'Spectral deconfounding via perturbed sparse linear models (Cevid et al.).'),
 'H03435': ('C2', '3b', R, 'SVA overfitting in brain gene expression (bioRxiv Oct 2026): high false-positive rates not revealed by standard diagnostics; no operator-gain analysis.'),
 'H04758': ('C2', '3b', R, 'Batch correction distorts gene-gene correlation (bioRxiv 2026).'),
 'H00842': ('C2', '3b', R, 'Adjusting for PCs induces replicating false positives (2017).'),
 'H00847': ('C2', '3b', R, 'Nygaard 2016 (positive control).'),
 'H00790': ('C2', '3b', R, 'Zindler 2020 (positive control).'),
 'H00964': ('C2', '3b', R, 'Wang 2017 (positive control).'),
 'H03523': ('C2', '3b', R, 'Gerard & Stephens (positive control).'),
 'H01047': ('C2', '3b', R, 'Gagnon-Bartsch & Speed 2012 RUV-2.'),
 'H01051': ('C2', '3b', R, 'RUV when neither unwanted variation nor factor of interest is observed (2015).'),
 'H03513': ('C2', '3b', R, 'RUV with negative controls tech report (2013).'),
 'H04918': ('C2', '3b', R, 'RUV with PRPS (2022).'),
 'H00908': ('C2', '3b', R, 'RUVprps software (2026).'),
 'H00749': ('C2', '3b', R, 'RUV in pseudobulk scRNA-seq (2025).'),
 'H02994': ('C2', '3b', R, 'Post-integrated inference with surrogate-control outcomes (Biometrika 2026).'),
 'H04864': ('C2', '3b', R, 'CytofRUV author response (record type).'),
 'H03877': ('C3', '3b', R, 'Batch correction with sample remeasurement in confounded designs (2023).'),
 'H05504': ('C3', '3b', R, 'Practical limits of batch-effect correction (2019).'),
 'H03898': ('C3', '3b', R, 'Lipschitz-constrained neural networks (component prior art).'),
}


def main():
    U = [json.loads(l) for l in open(os.path.join(C, 'triaged.jsonl'))]
    audit = set(json.load(open(os.path.join(C, 'audit_sample_uids.json'))))
    rows = []
    for u in U:
        src = '+'.join(u['sources'])
        cl = ';'.join(u['claims'])
        if u['uid'] in FLAG:
            c, st, dec, why = FLAG[u['uid']]
        elif not u['pass3a']:
            c, st, dec = cl, ('3c' if u['uid'] in audit else '3a'), 'EXCLUDE'
            why = ('audit sample: manually confirmed off-topic' if u['uid'] in audit else
                   'fails registered 3a rule (<2 concept blocks of retrieving claim): ' + ','.join(u['blocks']))
        else:
            c, st, dec, why = ';'.join(u['pass3a']), '3b', 'EXCLUDE', 'title/abstract off-topic for claim (' + u['why3a'] + ')'
        rows.append([u['uid'], (u['title'] or '').replace('\n', ' ').strip(), u['year'], src, c, st, dec, why,
                     u['doi'], u['pmid'], u['pmcid']])
    # records found only in stage-4 follow-up searches (not part of the registered strings)
    rows.append(['H-FU1', 'A ferroptosis-related transcriptional signature in the brain after a 35-day spaceflight', 2026,
                 'OpenAlex (follow-up)', 'C4', '4', R,
                 'npj Microgravity 2026; analyses OSD-682/685/698/699 (RR-18 brain), not the seven deposits. Found via the related-articles panel, not by a registered string.',
                 '10.1038/s41526-026-00649-2', '', ''])
    with open(os.path.join(OUT, 'H_screening.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'title', 'year', 'source', 'claim', 'stage reached', 'decision', 'reason', 'doi', 'pmid', 'pmcid'])
        w.writerows(rows)
    log = [json.loads(l) for l in open(os.path.join(C, 'search_log.jsonl'))]
    calls = int(open(os.path.join(C, 'call_count.txt')).read())
    json.dump({'protocol_sha256_file': 'H_PROTOCOL.sha256', 'total_api_calls': calls, 'cap': 600,
               'n_logged_queries': len(log), 'queries': log}, open(os.path.join(OUT, 'H_search_log.json'), 'w'), indent=1)
    from collections import Counter
    print('rows', len(rows), Counter(r[5] for r in rows), Counter(r[6] for r in rows), 'calls', calls, 'queries', len(log))


if __name__ == '__main__':
    main()
