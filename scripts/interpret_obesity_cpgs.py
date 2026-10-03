"""Describe frozen classifier CpGs without refitting or selecting on replication.

Platform annotation, published expression association and our observed group
differences are separate evidence layers. No expression data are reanalysed.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl

def sha(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f, 'sha256').hexdigest()

def main(root):
    root = Path(root); base = root / 'docs/chronic_20261003'
    out = base / 'obesity_annotation_v1'; out.mkdir(parents=True, exist_ok=True)
    inputs = [root / 'data/chronic_20261003/prepared_v3' / (s + '.npz') for s in ['obesity_discovery','obesity_replication']]
    model = root / 'models/chronic_20261003_obesity_discovery_supervised.joblib'
    hashes = {str(p.relative_to(root)): sha(p) for p in inputs + [model]}
    expected = json.loads((base / 'methylation_transfer_v3/obesity_discovery/supervised_summary.json').read_text())['model_artifact']['sha256']
    assert sha(model) == expected
    selected = pd.read_csv(base / 'methylation_transfer_v3/obesity_discovery/supervised_refit_features.tsv', sep='\t')
    assert selected.shape[0] == 50 and selected.feature_id.is_unique
    selected['source_coefficient_rank'] = selected.standardized_coefficient.abs().rank(method='first', ascending=False).astype(int)
    ann = pd.read_csv(out / 'selected_native_GPL13534.tsv', sep='\t', keep_default_na=False)
    result = selected.merge(ann, left_on='feature_id', right_on='ID', validate='one_to_one')
    arrays = [np.load(p, allow_pickle=False) for p in inputs]
    common = np.array(sorted(set(arrays[0]['probes']) & set(arrays[1]['probes'])))
    frequencies = Counter(); folds = []
    for seed in [20261003,20261020,20261103]:
        path = base / f'methylation_transfer_v3/obesity_discovery/supervised_{seed}_fold_audit.json'
        hashes[str(path.relative_to(root))] = sha(path)
        for fold in json.loads(path.read_text()):
            assert fold['all_learned_preprocessing_training_only']
            ids = common[np.array(fold['selected_indices'])]
            frequencies.update(ids.tolist()); folds.append({'seed':seed,'outer_fold':fold['outer_fold'],'selected_n':len(ids)})
    assert len(folds) == 15
    result['outer_selection_count'] = result.feature_id.map(lambda p: frequencies[p])
    result['outer_selection_denominator'] = 15
    long = []
    for cohort, a in zip(['discovery','replication'], arrays):
        lookup = {p:i for i,p in enumerate(a['probes'])}; y = a['y']; x = a['X'][:,[lookup[p] for p in result.feature_id]]
        assert x.shape == (48,50) and np.bincount(y).tolist() == [24,24]
        effects = []
        for j,p in enumerate(result.feature_id):
            control = x[y==0,j]; obese=x[y==1,j]
            control=control[np.isfinite(control)]; obese=obese[np.isfinite(obese)]
            rng=np.random.default_rng(20261003+j)
            lo=hi=delta=np.nan
            if len(control)>=2 and len(obese)>=2:
                samples = obese[rng.integers(len(obese),size=(2000,len(obese)))].mean(axis=1) - control[rng.integers(len(control),size=(2000,len(control)))].mean(axis=1)
                lo,hi = 100*np.quantile(samples,[.025,.975])
                delta=float(100*(obese.mean()-control.mean()))
            effects.append({'control_n':len(control),'obese_n':len(obese),'control_mean_beta':float(control.mean()) if len(control) else np.nan,'obese_mean_beta':float(obese.mean()) if len(obese) else np.nan,'delta_beta_pp':delta,'delta_beta_pp_ci_low':float(lo),'delta_beta_pp_ci_high':float(hi)})
            for sample,label,value in zip(a['samples'],y,x[:,j]):
                long.append({'cohort':cohort,'sample_id':sample,'obesity_label':int(label),'feature_id':p,'beta':float(value) if np.isfinite(value) else None})
        for key in effects[0]: result[cohort+'_'+key] = [v[key] for v in effects]
    workbook_path = root / 'data/chronic_20261003/obesity_annotation_v1/published_Supplementary_Data_1-26.xlsx'
    hashes[str(workbook_path.relative_to(root))] = sha(workbook_path)
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    published = []; matches = {}; nearest = {}
    for number in [2,7,9,11]:
        sheet = next(s for s in workbook if s.title.startswith(f'ST{number}_'))
        iterator = sheet.iter_rows(values_only=True); header = next(iterator)
        for excel_row, values in enumerate(iterator, start=2):
            r = dict(zip(header,values)); probe=r['sentinel']
            if probe not in set(result.feature_id): continue
            if number == 2:
                nearest[probe] = {'published_nearest_gene':r['Nearest_gene'] or '', 'published_nearest_region':r['Nearest_anno'] or ''}
            published.append({'source_doi':'10.1038/s41467-023-38439-z','sheet':sheet.title,'excel_row':excel_row,'evidence_kind':{2:'published obesity EWAS sentinel and nearest-location annotation',7:'published genic methylation-expression association',9:'published functional mapping and methylation-expression association',11:'published shared TAD and methylation-expression association'}[number],**{k:v for k,v in r.items() if k is not None}})
            matches.setdefault(number,set()).add(probe)
    # Serialize each published row without making RNA associations into our data.
    (out/'published_exact_probe_matches.json').write_text(json.dumps(published,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    for col in ['published_nearest_gene','published_nearest_region']:
        result[col]=result.feature_id.map(lambda p:nearest.get(p,{}).get(col,''))
    for number in [7,9,11]:
        subset = [r for r in published if r['sheet'].startswith(f'ST{number}_')]
        result[f'published_ST{number}_targets']=result.feature_id.map(lambda p:';'.join(sorted(set(r['target_gene'] for r in subset if r['sentinel']==p))))
    result['gene_symbols_native']=result.UCSC_RefGene_Name.map(lambda s:';'.join(dict.fromkeys(s.split(';'))) if s else '')
    result['gene_regions_native']=result.UCSC_RefGene_Group.map(lambda s:';'.join(dict.fromkeys(s.split(';'))) if s else '')
    result['cpg_island_context_native']=result.Relation_to_UCSC_CpG_Island.map(lambda s:s or 'Open sea (empty native field)')
    result['coordinate_build']=result.Genome_Build
    pairs=[]
    for _,r in result.iterrows():
        names=r.UCSC_RefGene_Name.split(';') if r.UCSC_RefGene_Name else []
        regions=r.UCSC_RefGene_Group.split(';') if r.UCSC_RefGene_Group else []
        assert len(names)==len(regions)
        for g,context in sorted(set(zip(names,regions))):
            pairs.append({'feature_id':r.feature_id,'native_gene_symbol':g,'native_refgene_context':context,'CHR':r.CHR,'MAPINFO':r.MAPINFO,'Genome_Build':r.Genome_Build,'evidence_kind':'platform genomic-location annotation; not proven regulatory targeting'})
    result=result.sort_values('source_coefficient_rank')
    # Trim oligonucleotide/hardware columns from the interpretation-facing table.
    cols=['source_coefficient_rank','feature_id','standardized_coefficient','CHR','MAPINFO','coordinate_build','gene_symbols_native','gene_regions_native','cpg_island_context_native','UCSC_RefGene_Accession','Probe_SNPs','Probe_SNPs_10','outer_selection_count','outer_selection_denominator']+[cohort+'_'+key for cohort in ['discovery','replication'] for key in effects[0]]+['published_nearest_gene','published_nearest_region','published_ST7_targets','published_ST9_targets','published_ST11_targets']
    result[cols].to_csv(out/'cpg_gene_interpretation_all50.tsv',sep='\t',index=False)
    result[cols].head(12).to_csv(out/'cpg_gene_interpretation_top12.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(out/'native_gene_context_pairs.tsv',sep='\t',index=False)
    pd.DataFrame(long).to_csv(out/'selected_cpg_beta_source_data.tsv',sep='\t',index=False)
    annotated = result.gene_symbols_native.ne('')
    summary={'actual_final_selected_probes':50,'native_platform_coordinate_matches':50,'native_gene_annotated_probes':int(annotated.sum()),'without_native_gene_assignment':int((~annotated).sum()),'unique_native_gene_symbols':sorted(set(g for s in result.gene_symbols_native if s for g in s.split(';'))),'same_delta_direction_probes':int((np.sign(result.discovery_delta_beta_pp)==np.sign(result.replication_delta_beta_pp)).sum()),'native_build_values':sorted(set(result.coordinate_build)),'selection_stability_denominator':15,'top12_rule':'Descending absolute coefficient of frozen source-only final refit; independent of replication effect/association results','published_exact_probe_match_counts':{str(n):len(matches.get(n,set())) for n in [2,7,9,11]},'mean_difference_interval':'Unadjusted obese minus control beta percentage points; separate group bootstrap 2000; conditional descriptive intervals, not selection-adjusted or age/sex-adjusted EWAS','posthoc_replication_interpretation':True,'RNA_reanalysed':False,'model_refitted':False,'causal_target_assignment_claimed':False,'input_hashes':hashes,'fold_counts':folds,'script_sha256':sha(__file__)}
    comparable=result.discovery_delta_beta_pp.notna() & result.replication_delta_beta_pp.notna()
    summary['two_cohort_observed_effect_comparable_probes']=int(comparable.sum())
    summary['insufficient_replication_observations']=result.loc[~comparable,['feature_id','replication_control_n','replication_obese_n']].to_dict('records')
    assert all(sha(root/p)==digest for p,digest in hashes.items())
    summary['all_inputs_and_frozen_model_unchanged']=True
    (out/'interpretation_receipt.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['input_hashes','fold_counts','script_sha256','unique_native_gene_symbols']},ensure_ascii=False))
    print(result[['feature_id','gene_symbols_native','published_nearest_gene','published_ST7_targets','published_ST9_targets','published_ST11_targets','discovery_delta_beta_pp','replication_delta_beta_pp','outer_selection_count']].head(12).to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root');main(p.parse_args().root)
