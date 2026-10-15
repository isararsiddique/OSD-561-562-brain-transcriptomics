#!/usr/bin/env python3
"""Apply the verified npj Microgravity structure to generated manuscript blocks.

Run after prepare_final_content.py and before build_final_documents.py.
Methods must reside in the main manuscript, including mathematical derivations,
simulation specifications and reproducibility instructions. SI contains results,
electronic table indexes and supplementary figures only.
"""
import json
import re
import shutil
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
R = HERE / 'results'
TITLE = 'Ground collection scores reveal limits of transcriptomic correction in mouse spaceflight brain'


def write_markdown(blocks, path):
    out = []
    for q in blocks:
        kind = q['type']
        if kind == 'title':
            out.append('# ' + q['text'])
        elif kind == 'heading':
            out.append('#' * (q['level'] + 1) + ' ' + q['text'])
        elif kind in ('paragraph', 'equation'):
            out.append(q['text'])
        elif kind == 'table':
            out += [q['title'], '| ' + ' | '.join(q['headers']) + ' |',
                    '| ' + ' | '.join(['---'] * len(q['headers'])) + ' |']
            out += ['| ' + ' | '.join(str(v).replace('|', ' / ').replace('\n', ' ')
                                     for v in row) + ' |' for row in q['rows']]
        elif kind == 'figure':
            out += [f"![{q['filename']}](publication_figures/{q['filename']}.png)", q['caption']]
    path.write_text('\n\n'.join(out) + '\n')


def main():
    main_blocks = json.loads((R / 'manuscript_content.json').read_text())
    supplement = json.loads((R / 'supplement_content.json').read_text())
    original_title = main_blocks[0]['text']
    for blocks in (main_blocks, supplement):
        for q in blocks:
            if 'text' in q:
                q['text'] = q['text'].replace(original_title, TITLE)
    abstract_index = next(i for i, q in enumerate(main_blocks) if q.get('text') == 'Abstract') + 1
    main_blocks[abstract_index]['text'] = (
        'Collection conditions can overlap biological responses in brain transcriptomics. '
        'We reanalysed seven NASA deposits containing 129 libraries and 83 animals using '
        'animal-restricted inference, paired regional tests and count models. Ground collection-score '
        'discrimination reproduced in both core regions; refitted label enumeration gave exact '
        'p=0.005 and BY q=0.039 across 18 tests. Across 362 overlapping ground holdouts '
        'with training-only preprocessing, original loadings worsened outside-sensor prediction '
        'by 9.1% in cerebellum and 3.5% in RR-23. Ground-derived subtraction was an oblique '
        'projection with spectral gain 8.13–15.31. A closest-loading constraint bounded '
        'amplification while retaining score removal. In 2,600 simulated datasets, bounding '
        'moderated instability under changed coupling but could leave substantial nuisance. '
        'No paired panel or cross-cohort gene test survived correction. These results distinguish '
        'reproducible collection signatures from transferable correction and motivate explicit stability validation.'
    )
    # Publisher requests funding within Acknowledgements rather than a separate
    # Funding section; the user-confirmed factual declaration is preserved.
    funding_index = next(i for i, q in enumerate(main_blocks) if q.get('text') == 'Funding')
    funding_text = main_blocks[funding_index + 1]['text']
    del main_blocks[funding_index:funding_index + 2]
    acknowledgements = next(i for i, q in enumerate(main_blocks) if q.get('text') == 'Acknowledgements')
    main_blocks[acknowledgements + 1]['text'] += ' ' + funding_text
    start = next(i for i, q in enumerate(supplement)
                 if q.get('text') == 'Supplementary Methods')
    end = next(i for i, q in enumerate(supplement)
               if q.get('text') == 'S7. Full table and software index')
    repro_start = next(i for i, q in enumerate(supplement)
                       if q.get('text') == 'S8. Reproduction order')
    repro_end = next(i for i in range(repro_start, len(supplement))
                     if supplement[i]['type'] == 'figure')
    methods = supplement[start + 1:end] + supplement[repro_start:repro_end]
    # The holdout/QC summaries are supplementary results, already documented in
    # the main Methods; retain these in the SI rather than duplicating results.
    ground_start = next(i for i, q in enumerate(methods)
                        if q.get('text') == 'Ground holdout and all deposit structure extensions')
    reproduction = next(i for i, q in enumerate(methods)
                        if q.get('text') == 'S8. Reproduction order')
    supplemental_results = methods[ground_start:reproduction]
    methods = methods[:ground_start] + methods[reproduction:]
    number_map = {'S1': '6', 'S2': '7', 'S3': '8'}
    for q in methods:
        if q['type'] == 'heading':
            q['text'] = re.sub(r'^S\d+\.\s*', '', q['text'])
        if q['type'] == 'equation':
            q['text'] = re.sub(r'\((S[123])\)', lambda m: '(' + number_map[m[1]] + ')', q['text'])
        if q['type'] == 'table':
            q['title'] = re.sub(r'^Supplementary Table S([123])\.',
                                lambda m: 'Table ' + str(int(m[1]) + 1) + '.', q['title'])
    for q in main_blocks:
        if 'text' in q:
            q['text'] = q['text'].replace('Supplementary Methods', 'Methods')
            q['text'] = q['text'].replace('Supplementary Tables', 'Supplementary Data')
        if q.get('text') == 'Reporting, software and AI assistance':
            q['text'] = 'Statistics and reproducibility'
        if q.get('text') == 'Author contributions':
            contribution_start = main_blocks.index(q) + 1
    initials = {'Isarar Siddique': 'I.S.', 'Raiyan Siddique': 'R.S.',
                'Arshi Niyaz': 'A.N.', 'Jamyang Ugyen Tshomo': 'J.U.T.', 'Atish Raza': 'A.R.'}
    for q in main_blocks[contribution_start:contribution_start + 5]:
        name = q['text'].split(':', 1)[0]
        q['text'] = q['text'].replace(name + ':', initials[name] + ' (' + name + '):', 1)
    insertion = next(i for i, q in enumerate(main_blocks)
                     if q.get('text') == 'Data availability')
    main_blocks[insertion:insertion] = methods
    supplement = (supplement[:start] +
                  [{'type': 'heading', 'text': 'Supplementary results', 'level': 1}] +
                  supplemental_results + supplement[end:repro_start] + supplement[repro_end:])
    for q in supplement:
        if q.get('text') == 'S7. Full table and software index':
            q['text'] = 'Electronic results and software index'
        if q['type'] == 'table':
            q['title'] = q['title'].replace('Supplementary Table S4.', 'Supplementary Table S1.')
            if 'Electronic table' in q['title']:
                for row in q['rows']:
                    if row[0] == 'Figure source and captions':
                        row[1] = ('publication_figures/*.svg; publication_figures/*.png; '
                                  'results/publication_figure_manifest.json. Intermediate sources: figures/.')
                q['rows'].append(['Journal structure and independent audit',
                                  'scripts/finalize_for_journal.py; extensions/independent_audit_20261015/; extensions/journal_audit_20261015/'])
    repro = next(q for q in methods if q['type'] == 'paragraph' and
                 q['text'].startswith('Create a Python 3.12 environment'))
    repro['text'] = (
        'Create a Python 3.12 environment and install requirements.txt. Limit BLAS threads '
        'with OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1. '
        'The README supplies the ordered scientific reconstruction commands, including '
        'forced count-model refitting and all retained robustness analyses. After analysis, '
        'run scripts/prepare_manuscript_sources.py, scripts/build_draft.py --figures-only, '
        'scripts/build_extended_figures.py, scripts/refine_panel_layouts.py, '
        'scripts/build_graphical_abstract.py, scripts/prepare_final_content.py, '
        'scripts/finalize_for_journal.py and scripts/build_final_documents.py, in that order. '
        'The graphical overview is optional. Regenerate equations with scripts/render_equations.cjs '
        'using Node, mathjax-full 3.2.2 and sharp when their sources change. Run '
        'scripts/refine_journal_artwork.cjs before building Word files. Render the Word '
        'files to PDF and page images using LibreOffice and render_docx.py. Finally run '
        'the unit tests, scripts/validate_revision.py and scripts/validate_fault_injections.py. '
        'Full reconstruction retains every test and failed method; cached count fits are '
        'distinct from forced refitting. A sealed file manifest must be rebuilt after reviewed changes.'
    )
    audit_dir = HERE / 'extensions/independent_audit_20261015'
    exact = pd.read_csv(audit_dir / 'exact_ground_sensor_association.csv')
    sensor = exact[exact.panel_or_comparator == 'Sensor31']
    assert len(exact) == 18 and len(sensor) == 2
    assert (sensor.p_exact_complete_label_enumeration == .005).all()
    summary = json.loads((audit_dir / 'audit_summary.json').read_text())
    assert summary['status'] == 'PASS' and summary['n_failed'] == 0, summary
    insertion = next(i for i, q in enumerate(main_blocks)
                     if q.get('text') == 'Real-data attenuation is heterogeneous and uncertain')
    main_blocks[insertion:insertion] = [dict(type='paragraph', text=(
        'An additional exploratory finite-sample test rebuilt supervised panel weights and '
        'training-only preprocessing under all 400 age-restricted ground-cohort label assignments '
        'per core region. The statistic aggregated all 81 balanced holdouts once for each '
        'assignment, avoiding split-level independence. Sensor31 retained AUC 1 in both '
        'regions (exact p=0.005, BH q=0.01125 and BY q=0.03932 across all 18 inherited '
        'panel/comparator-by-region tests; Supplementary Fig. 5). Shared regional animals '
        'and observational exchangeability still limit independent replication and causality. '
        'The result is not Sensor31-specific: rPRG also reached AUC 1 in both regions. '
        'The collection association survives this dependence sensitivity; it does not validate '
        'transfer of correction loadings.'))]
    insertion = next(i for i, q in enumerate(main_blocks)
                     if q.get('text') == 'Ground anchors, transfer and resampling')
    main_blocks[insertion:insertion] = [dict(type='heading', level=2,
        text='Independent finite-sample and age-model robustness'), dict(type='paragraph', text=(
        'The post hoc audit fixed the eight inherited panels plus the whole-genome comparator '
        'in each core region before executing 18 tests. All 400 assignments preserve three '
        'positive labels among six animals in each age block. Each assignment repeats all '
        '81 train-eight/test-four balanced splits; filtering, geometric-mean references and '
        'supervised contrast weights are refitted from training animals. Label-independent '
        'training-reference calculations are cached by animal training set, giving 225 distinct '
        'training sets per region. Mean within-age held-out pair discrimination is recomputed '
        'for every assignment, and the exact upper-tail proportion includes the observed '
        'assignment. BH and arbitrary-dependence BY adjustments include all 18 tests. '
        'No flight labels or held-out outcomes select a panel or tuning parameter. '
        'Separately, all eight panels in both regions are estimated from full age-by-flight-by-collection '
        'cell means and a saturated eight-column model. The 64 age-specific, equal-age '
        'and age-difference HC3 normal estimands form a separate exploratory BH family. '
        'Empty metadata-excluded cells remain unavailable. Saved-count reconstructions also '
        'decompose prediction error and report every training-expression tertile and animal, '
        'without treating overlapping holdouts as independent observations.'))]
    insertion = next(i for i, q in enumerate(main_blocks)
                     if q.get('text') == 'Ground-only prediction validation')
    # The saturated model is direction robustness, never a promotion of the
    # weaker metadata-sensitive biological significance result.
    result_position = next(i for i, q in enumerate(main_blocks)
                           if q.get('text') == 'Ground-derived correction has an oblique geometry') if any(
        q.get('text') == 'Ground-derived correction has an oblique geometry' for q in main_blocks) else None
    age_text = ('The independent saturated-age analysis retained negative hippocampal '
                'Circadian interactions at both ages (−0.510 log2 at 12 weeks and −0.307 '
                'at 29 weeks), with an equal-age estimate of −0.409. This agrees in direction '
                'with the additive-age estimate, but the equal-age BH q was 0.113 across '
                'all 64 saturated estimands and no result survived that family. Sparse '
                'metadata-excluded age cells remain a limitation. This is model-direction '
                'robustness, not additional corrected biological significance.')
    if result_position is None:
        result_position = next(i for i, q in enumerate(main_blocks)
                               if q['type'] == 'heading' and 'oblique' in q.get('text', '').lower())
    main_blocks.insert(result_position, dict(type='paragraph', text=age_text))
    supplement.insert(4, dict(type='paragraph', text=(
        'The independent audit supplies all 18 exact-label tests, both 400-assignment null '
        'tables, 64 saturated-age estimands, per-animal and expression-tertile prediction '
        'diagnostics, source checks and the fixed exploratory protocol under '
        'extensions/independent_audit_20261015/. No saturated-age test survives its BH '
        'family. Cerebellar original-score prediction worsens on average in each '
        'training-expression tertile (MSE ratios 1.096, 1.090 and 1.065), as does RR-23 '
        '(1.045, 1.008 and 1.025); individual animals remain heterogeneous.')))
    for block in main_blocks:
        if block['type'] == 'figure' and block['filename'].startswith('figure2_'):
            block['caption'] = block['caption'].replace('Unadjusted panel effects',
                'Panel effects before ground-derived transformation')
            block['caption'] += ' Age and, for pooled core contrasts, tissue nuisance adjustment are retained.'
        if block['type'] == 'figure' and block['filename'].startswith('figure10_'):
            block['caption'] += ' Translucent violins show the full 100-repetition distributions where variable; points are means.'
    statistics_index = next(i for i, q in enumerate(main_blocks)
                            if q.get('text') == 'Statistics and reproducibility')
    main_blocks.insert(statistics_index + 1, dict(type='paragraph', text=(
        'The significance threshold is 0.05. Primary panel contrasts, gradients, '
        'gene Wald tests, paired contrasts and wild interaction tests are two-sided; '
        'their exact probabilities and complete multiplicity families are reported. '
        'The independent collection-classification test uses an upper-tail AUC '
        'statistic with orientation learned exclusively from training labels. '
        'Its label assignments, rather than overlapping splits, define the null '
        'enumeration. Reanalysis sample sizes are fixed by the deposited animals. '
        'Error bars identify animal-bootstrap percentiles, approximate HC3 intervals '
        'or Monte Carlo mean uncertainty as appropriate; no post hoc significance '
        'threshold or sample-size extension was selected to obtain a favorable result.')))
    # Include the independent audit's complete multipanel result in the SI.
    figure_name = 'figureS5_independent_ground_and_age_validation'
    q = dict(type='figure', filename=figure_name, png_dpi=600, vector_format='SVG',
        source_filename='extensions/independent_audit_20261015/independent_audit_panel',
        caption=('Supplementary Figure 5. Independent finite-sample and prediction robustness. '
                 '(a,b) Refitted age-restricted label null distributions and observed Sensor31 '
                 'AUC in each core region. All 400 assignments enter each exact test; '
                 'all 18 panel/comparator tests enter BH and BY adjustment. '
                 '(c) Outside-sensor original-score prediction error by training-expression '
                 'tertile across all held-out animals. Points are means and bars are '
                 '2.5th–97.5th overlapping-split ranges, not independent-sample confidence limits. '
                 'Ratios compare prediction with the training-fitted intercept/age baseline. '
                 '(d) Saturated-age regional Circadian estimates with HC3 normal '
                 'intervals. No estimate passes BH across all 64 saturated estimands. '
                 'Exact collection association does not establish causal nuisance or biological recovery.'))
    for ext in ('png', 'svg'):
        shutil.copy2(audit_dir / ('independent_audit_panel.' + ext),
                     HERE / 'publication_figures' / (figure_name + '.' + ext))
    supplement.append(q)
    manifest = [x for x in main_blocks + supplement if x['type'] == 'figure']
    (R / 'publication_figure_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    # Tables follow the reference list in the main text file, before the figure
    # pages. Their cited numbering and scientific content are unchanged.
    tables = [q for q in main_blocks if q['type'] == 'table']
    main_blocks = [q for q in main_blocks if q['type'] != 'table']
    first_figure = next(i for i, q in enumerate(main_blocks) if q['type'] == 'figure')
    main_blocks[first_figure:first_figure] = tables
    # Extend the same auditable numeric registry used for the original text.
    registry = pd.read_csv(R / 'manuscript_numeric_lookups.csv')
    new_lookups = []
    source = 'extensions/independent_audit_20261015/exact_ground_sensor_association.csv'
    for _, row in sensor.iterrows():
        for field in ('mean_crossfitted_age_stratified_AUC',
                      'p_exact_complete_label_enumeration', 'q_BH_complete_audit_18_tests',
                      'q_BY_complete_audit_18_tests'):
            new_lookups.append(dict(source_file=source,
                filters=json.dumps({'accession': row.accession, 'panel_or_comparator': 'Sensor31'}, sort_keys=True),
                field=field, value=float(row[field])))
    age = pd.read_csv(audit_dir / 'saturated_age_panel_interactions.csv')
    for _, row in age[(age.accession == 'OSD-562') & (age.panel == 'Circadian') &
                      age.estimand.isin(['young', 'old', 'equal_age_mean'])].iterrows():
        for field in ('effect_log2', 'q_BH_all_64_saturated_estimands'):
            new_lookups.append(dict(source_file='extensions/independent_audit_20261015/saturated_age_panel_interactions.csv',
                filters=json.dumps({'accession': row.accession, 'panel': row.panel, 'estimand': row.estimand}, sort_keys=True),
                field=field, value=float(row[field])))
    pd.concat([registry, pd.DataFrame(new_lookups)]).drop_duplicates().to_csv(
        R / 'manuscript_numeric_lookups.csv', index=False)
    for name, blocks in [('manuscript_content.json', main_blocks),
                         ('supplement_content.json', supplement)]:
        (R / name).write_text(json.dumps(blocks, indent=2, ensure_ascii=False) + '\n')
    write_markdown(main_blocks, HERE / 'MANUSCRIPT_COMPLETE.md')
    write_markdown(supplement, HERE / 'SUPPLEMENTARY_INFORMATION.md')
    report = {'all_methods_in_main': True,
              'title_words': len(TITLE.split()),
              'abstract_words': len(main_blocks[abstract_index]['text'].split()),
              'main_equations': sum(q['type'] == 'equation' for q in main_blocks),
              'supplement_equations': sum(q['type'] == 'equation' for q in supplement),
              'publisher_url': 'https://www.nature.com/npjmgrav/submission-guidelines',
              'guidance_checked': '2026-10-15'}
    (R / 'journal_format_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
