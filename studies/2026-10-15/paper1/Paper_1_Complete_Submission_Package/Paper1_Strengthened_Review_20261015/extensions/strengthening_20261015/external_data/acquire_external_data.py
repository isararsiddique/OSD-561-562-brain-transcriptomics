"""Acquire design-selected external mouse-brain STAR counts without discovery selection.

Run with bundled Python and scientific PYTHONPATH. No normalized/VST/DE table is used.
All output is confined to this script's directory, or an explicitly given output path.
"""
from pathlib import Path
import argparse, csv, io, json, re, shutil, hashlib, datetime, zipfile
from urllib.parse import urljoin
from urllib.request import urlopen, Request
import numpy as np
import pandas as pd

DEFAULT_SOURCE = Path('/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/research/brain_integrative_2026')
SELECTED = (202, 352, 374, 611)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(1024 * 1024): h.update(b)
    return h.hexdigest()

def isa_table(z, name):
    rows = list(csv.reader(io.StringIO(z.read(name).decode('utf-8-sig')), delimiter='\t'))
    headers = rows[0]
    out = []
    for row in rows[1:]:
        if not any(row): continue
        row += [''] * (len(headers) - len(row))
        d = {}
        for i, h in enumerate(headers):
            if not h or h in ('Unit', 'Term Source REF', 'Term Accession Number', 'Protocol REF'): continue
            if h in d and d[h] != row[i]: raise ValueError(f'Conflicting duplicate field {h}')
            d[h] = row[i]
            if i + 1 < len(headers) and headers[i + 1] == 'Unit': d[h + '__unit'] = row[i + 1]
        out.append(d)
    return out

def canonical(s): return re.sub(r'\s+', '_', str(s).strip())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', type=Path, default=DEFAULT_SOURCE)
    ap.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument('--online', action='store_true')
    args = ap.parse_args(); out = args.output; out.mkdir(parents=True, exist_ok=True)
    provenance = []; inventory = []; all_meta = []; all_members = []; unmatched = []
    def copy_file(src, dst, accession, url='', expected=None):
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.exists():
            if src.resolve() != dst.resolve(): shutil.copy2(src, dst)
            mode = 'local NASA source file copied without numeric transformation'
        elif dst.exists(): mode = 'existing immutable acquisition copy'
        elif args.online and url:
            with urlopen(Request(url, headers={'User-Agent': 'brain-count-source-audit/1.0'}), timeout=90) as r:
                b = r.read()
            if b[:100].lstrip().lower().startswith((b'<html', b'<!doctype')): raise ValueError('HTML instead of data')
            dst.write_bytes(b); mode = 'NASA HTTP acquisition'
        else: raise FileNotFoundError(str(src))
        if expected is not None and dst.stat().st_size != expected: raise ValueError(f'NASA size mismatch: {dst}')
        provenance.append(dict(accession=accession, source_path=str(src), copied_path=str(dst), url=url,
                               sha256=sha(dst), bytes=dst.stat().st_size, expected_nasa_bytes=expected,
                               mode=mode, recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        return dst
    for n in SELECTED:
        accession = f'OSD-{n}'; raw = out / accession / 'raw_sources'; raw.mkdir(parents=True, exist_ok=True)
        for kind in ('meta', 'files'):
            src = args.source_root / 'metadata' / f'{accession}_{kind}.json'
            copy_file(src, raw / src.name, accession, f'https://osdr.nasa.gov/osdr/data/osd/{kind}/{n}')
        metaj = json.loads((raw / f'{accession}_meta.json').read_text())
        study = metaj['study'][accession]['studies'][0]
        filej = json.loads((raw / f'{accession}_files.json').read_text())
        files = filej['studies'][accession]['study_files']
        selected_files = [f for f in files if ('STAR_Unnormalized_Counts_GLbulkRNAseq.csv' in f['file_name'] or
                          'SampleTable_GLbulkRNAseq.csv' in f['file_name'] or f['file_name'].endswith('-ISA.zip'))]
        for f in selected_files:
            copy_file(args.source_root/'data'/accession/f['file_name'], raw/f['file_name'], accession,
                      urljoin('https://osdr.nasa.gov', f['remote_url']), f['file_size'])
        counts_path = next(raw.glob('*STAR_Unnormalized_Counts_GLbulkRNAseq.csv'))
        counts = pd.read_csv(counts_path, index_col=0)
        if not counts.index.is_unique or not counts.columns.is_unique: raise ValueError('Duplicate count identifiers')
        x = counts.to_numpy(dtype=float)
        if not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
            raise ValueError('Counts are not finite nonnegative integers')
        counts = counts.astype('int64'); counts.index.name = 'gene_id'
        counts.to_csv(out/accession/'counts.csv.gz', compression={'method':'gzip','mtime':0})
        z = zipfile.ZipFile(next(raw.glob('*-ISA.zip')))
        sname = next(f for f in z.namelist() if f.startswith('s_'))
        samples = isa_table(z, sname)
        aname = next(f for f in z.namelist() if 'rna-sequencing-(rna-seq)' in f and 'single-cell' not in f)
        assays = isa_table(z, aname)
        (out/accession/'isa_sample_records.json').write_text(json.dumps(samples, indent=2))
        (out/accession/'isa_assay_records.json').write_text(json.dumps(assays, indent=2))
        (out/accession/'study_protocols.json').write_text(json.dumps(study['protocols'], indent=2))
        sample_by_id = {canonical(s['Sample Name']):s for s in samples}
        sm = pd.read_csv(next(raw.glob('*SampleTable_GLbulkRNAseq.csv')),index_col=0)
        cohort = []
        for col in counts.columns:
            matches = [a for a in assays if canonical(a['Sample Name']) == col or
                       any(f'_{col}_GLbulkRNAseq' in v for k,v in a.items() if 'Merged Sequence Data' in k)]
            if not matches: raise ValueError(f'No assay mapping for {accession}/{col}')
            source_samples = [sample_by_id[canonical(a['Sample Name'])] for a in matches]
            sample = source_samples[0]; ass = matches[0]
            def get(*ks): return next((sample[k] for k in ks if sample.get(k)), '')
            for ss in source_samples:
                for k in ['Characteristics[Sex]','Characteristics[Strain]','Factor Value[Ionizing Radiation]']:
                    if ss.get(k) != sample.get(k): raise ValueError(f'Pool has inconsistent {k}')
            srcs = list(dict.fromkeys(s['Source Name'] for s in source_samples))
            if n == 611:
                # Original protocol says two hippocampi per library; ISA gives one pool label only.
                replicate_unit = 'two-animal hippocampal pool; member animal IDs unavailable'
                known_animals = ''; pool_size = 2; animal = f'{accession}:pool:{col}'
            else:
                pool_size = len(srcs); known_animals = pool_size
                replicate_unit = 'animal' if pool_size == 1 else 'disjoint two-animal hippocampal pool'
                animal = f'{accession}:{srcs[0]}' if pool_size == 1 else f'{accession}:pool:{col}'
            timekey = next((k for k in ['Factor Value[Time of Sample Collection After Treatment]',
                         'Parameter Value[Time of Sample Collection After Treatment]'] if sample.get(k)), '')
            agekey = next((k for k in ['Characteristics[Age at Launch]','Characteristics[age]'] if sample.get(k)), '')
            rad = get('Factor Value[Ionizing Radiation]'); dose = get('Factor Value[Absorbed Radiation Dose]',
                      'Parameter Value[Absorbed Radiation Dose]','Parameter Value[absorbed radiation dose]')
            dosekey=next((k for k in ['Factor Value[Absorbed Radiation Dose]', 'Parameter Value[Absorbed Radiation Dose]', 'Parameter Value[absorbed radiation dose]'] if sample.get(k)), '')
            doseunit=sample.get(dosekey+'__unit','')
            unitfactor={'gray':1.,'Gy':1.,'milligray':.001,'mGy':.001,'microgray':.000001}.get(doseunit)
            try: dosegy=float(dose)*unitfactor if unitfactor is not None else ''
            except ValueError: dosegy=''
            rawfiles = sorted(set(a.get('Raw Data File','') for a in matches))
            lanes = sorted(set(re.findall(r'_L(\d{3})_', ','.join(rawfiles))))
            row = dict(accession=accession, sample_id=col, animal_id=animal, source_name='|'.join(srcs),
                       source_names_json=json.dumps(srcs), independent_unit_id=animal, replicate_unit=replicate_unit,
                       pool_size_protocol=pool_size, identified_animal_members=known_animals,
                       organism=get('Characteristics[Organism]'), strain=get('Characteristics[Strain]'),
                       sex=get('Characteristics[Sex]'), age=get(agekey) if agekey else '',
                       age_unit=sample.get(agekey+'__unit',''), tissue=get('Characteristics[Material Type]'),
                       condition=(sm.loc[col,'condition'] if col in sm.index else f'{rad}|{dose} gray; absent from NASA SampleTable'),
                       present_in_nasa_sampletable=col in sm.index,
                       primary_analysis_included=col in sm.index,
                       spaceflight=get('Factor Value[Spaceflight]'),
                       radiation=rad, radiation_dose_value=dose, radiation_dose_unit=doseunit, radiation_dose_gy=dosegy,
                       hindlimb_unloading=get('Factor Value[Hindlimb Unloading]'),
                       collection_after_treatment=get(timekey) if timekey else '',
                       collection_after_treatment_unit=sample.get(timekey+'__unit',''),
                       duration=get('Parameter Value[duration]'), duration_unit=sample.get('Parameter Value[duration]__unit',''),
                       habitat=get('Parameter Value[Habitat]'), sacrifice_date=get('Comment[Sacrifice Date]'),
                       euthanasia=get('Parameter Value[Euthanasia Method]','Parameter Value[Euthanasia Chemical]'),
                       preservation=get('Parameter Value[Sample Preservation Method]'),
                       storage_temperature=get('Parameter Value[Sample Storage Temperature]'),
                       library_batch=ass.get('Parameter Value[Library Batch Number]',''), sequencing_lanes='|'.join(lanes),
                       sequencing_instrument=ass.get('Parameter Value[Sequencing Instrument]',''),
                       raw_data_files='|'.join(rawfiles), rnaseq_library_layout=ass.get('Parameter Value[Library Layout]',''),
                       rin=ass.get('Parameter Value[QA Score]',''), mission='RR-3' if n==352 else 'Earth radiation experiment',
                       independence_from_original='separate mission/experiment; no original RRRM-2/RR-10/RR-23 animals',
                       identity_warning='unknown members of two-animal pool' if n==611 else '',
                       anatomical_note='bulk pooled cryosections within one animal; same animals as RR-3 spatial/multiome' if n==352 else '')
            cohort.append(row); all_meta.append(row)
            for source in srcs:
                all_members.append(dict(accession=accession,sample_id=col,independent_unit_id=animal,
                                        source_name=source,actual_animal_id_known=(n!=611)))
        md = pd.DataFrame(cohort); md.to_csv(out/accession/'metadata.csv',index=False)
        used = {canonical(a['Sample Name']) for col in counts.columns for a in assays if canonical(a['Sample Name'])==col or
                any(f'_{col}_GLbulkRNAseq' in v for k,v in a.items() if 'Merged Sequence Data' in k)}
        for a in assays:
            if canonical(a['Sample Name']) not in used:
                unmatched.append(dict(accession=accession,sample_id=a['Sample Name'],reason='ISA RNAseq record not present in deposited STAR count matrix'))
        inventory.append(dict(accession=accession,title=study['title'],status='included_design_selected',
                              assay='bulk RNAseq STAR integer unnormalized counts', n_genes=len(counts),n_count_columns=counts.shape[1],
                              n_independent_units=md['independent_unit_id'].nunique(),n_identified_source_labels=len(set(md.source_name)),
                              experimental_units='animals' if n in (202,352) else 'two-animal pools',
                              conditions=json.dumps(md.groupby('condition',dropna=False).size().to_dict()),
                              tissue='|'.join(sorted(set(md.tissue))),sex='|'.join(sorted(set(md.sex))),
                              strain='|'.join(sorted(set(md.strain))),mission='RR-3' if n==352 else 'Earth analogue',
                              counts_path=str(out/accession/'counts.csv.gz'),metadata_path=str(out/accession/'metadata.csv'),
                              selection_reason='adult mouse brain bulk RNAseq counts and exposure/control design available; independent of original cohorts'))
        print('READY',accession,counts.shape,dict(md.groupby('condition').size()),flush=True)
    pd.DataFrame(all_meta).to_csv(out/'harmonized_metadata.csv',index=False)
    pd.DataFrame(all_members).to_csv(out/'sample_to_source_membership.csv',index=False)
    pd.DataFrame(inventory).to_csv(out/'cohort_inventory.csv',index=False)
    pd.DataFrame(unmatched).to_csv(out/'unavailable_count_samples.csv',index=False)
    (out/'acquisition_provenance.json').write_text(json.dumps(provenance,indent=2))
    output_manifest=[]
    for p in sorted(out.glob('OSD-*/*')):
        if p.is_file(): output_manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
    (out/'output_manifest.json').write_text(json.dumps(output_manifest,indent=2))

if __name__ == '__main__': main()
