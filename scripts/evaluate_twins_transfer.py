"""Frozen isolated-adipocyte model applied to BMI-discordant twin bulk SAT.

The endpoint is higher versus lower BMI within a discordant monozygotic twin
pair. It is not an obese-versus-normoweight clinical endpoint. Entire families
are the uncertainty resampling unit; the 11 concordant pairs remain descriptive.
Only the frozen primary model's 50 measured probes are extracted from complete
official sample tables. No target fit/selection/calibration is performed.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import gzip
import json
import os
from pathlib import Path
import re
import sys
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(v,'1')
import joblib
import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.evaluation import performance,exact_binomial_interval
from nutriomics_methylation.io import sha256,versions,write_json

def extract(path,selected):
    records={};values={};current=None;in_table=False;columns=None;row_counts={};definitions={};series={};detection_failed={};table_columns={}
    with gzip.open(path,'rt',encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('!Series_'):
                key,_,value=line.strip().partition(' = ');series.setdefault(key,[]).append(value)
            if line.startswith('^SAMPLE = '):
                current=line.strip().split(' = ',1)[1]
                if current in records:raise ValueError('Duplicate native sample accession')
                records[current]={'sample_id':current};values[current]={};row_counts[current]=0;detection_failed[current]=0;columns=None;continue
            if current is None:continue
            if line.startswith('!sample_table_begin'):in_table=True;columns=None;continue
            if line.startswith('!sample_table_end'):in_table=False;continue
            if line.startswith('!Sample_'):
                field,_,value=line.strip().partition(' = ')
                if field=='!Sample_characteristics_ch1':
                    name,separator,val=value.partition(':')
                    if separator:records[current]['characteristic_'+re.sub(r'\W+','_',name.lower().strip())]=val.strip()
                else:records[current].setdefault(field.removeprefix('!Sample_'),[]).append(value)
            if line.startswith('#'):definitions.setdefault(current,[]).append(line.strip())
            if in_table:
                fields=next(csv.reader([line],delimiter='\t'))
                if columns is None:
                    columns=fields
                    if 'ID_REF' not in columns or 'VALUE' not in columns:raise ValueError(f'Native sample table has no CpG/beta column: {columns}')
                    idcol,valcol=columns.index('ID_REF'),columns.index('VALUE')
                    detection=[i for i,c in enumerate(columns) if re.sub(r'\W+','',c.lower()).startswith('detectionp')]
                    if len(detection)!=1:raise ValueError('Expected exactly one deposited detection-P column')
                    pcol=detection[0];table_columns[current]=columns;continue
                if len(fields)!=len(columns):raise ValueError('Numeric sample table width mismatch')
                row_counts[current]+=1;probe=fields[idcol]
                if probe in selected:
                    if probe in values[current]:raise ValueError('Duplicate selected CpG in native sample table')
                    v=float(fields[valcol]) if fields[valcol] not in ['null','NA','NaN',''] else np.nan
                    if np.isfinite(v) and not 0<=v<=1:raise ValueError('Expected source beta in [0,1]')
                    detection_p=float(fields[pcol]) if fields[pcol] not in ['null','NA','NaN',''] else np.nan
                    if np.isfinite(detection_p) and not 0<=detection_p<=1:raise ValueError('Detection P outside [0,1]')
                    if not np.isfinite(detection_p) or detection_p>.01:
                        v=np.nan;detection_failed[current]+=1
                    values[current][probe]=v
    # A complete gzip read verifies the trailer CRC, including all unselected
    # numeric rows, without calculating an outcome-dependent target filter.
    meta=pd.DataFrame.from_dict(records,orient='index');meta.index.name='sample_id'
    if len(meta)!=70:raise ValueError('Expected 70 twin samples')
    native=meta['characteristic_twin_id'].str.extract(r'^Twin(\d+)_(\d+)$')
    if native.isna().any().any():raise ValueError('Unresolved native family IDs')
    meta['family_id']=native[0].astype(int);meta['twin_member']=native[1].astype(int)
    counts=meta.groupby('family_id').size()
    if len(counts)!=35 or not (counts==2).all():raise ValueError('Expected 35 intact twin pairs')
    meta['declared_bmi_group']=meta['characteristic_bmi'].str.lower()
    meta['declared_disease_state']=meta['characteristic_disease_state']
    meta['discordant']=meta.declared_disease_state.str.startswith('BMI_Disc_')
    families=meta.groupby('family_id').discordant.agg(['all','any'])
    if not (families['all']==families['any']).all() or int(families['all'].sum())!=24:raise ValueError('Expected 24 discordant and 11 concordant native families')
    primary=meta.loc[meta.discordant]
    if set(primary.declared_bmi_group)!={'high','low'} or not (primary.groupby(['family_id','declared_bmi_group']).size()==1).all():raise ValueError('Discordant high/low endpoint does not pair exactly')
    target=pd.DataFrame.from_dict(values,orient='index').reindex(index=meta.index,columns=selected)
    return target,meta,{'all_sample_table_rows':row_counts,'numeric_value_definitions':definitions,'series_metadata':series,'gzip_full_read_crc_passed':True,
      'native_table_columns':table_columns,'selected_beta_failed_detection_p_by_sample':detection_failed,'detection_p_threshold':.01,
      'all_50_selected_native_beta_range_validated':True,'other_numeric_rows_width_checked_not_all_beta_ranges':True}

def bh(p):
    p=np.asarray(p,float);result=np.full(len(p),np.nan);valid=np.flatnonzero(np.isfinite(p));order=valid[np.argsort(p[valid],kind='stable')]
    if len(order):result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return result

def run(root,source,out):
    if out.exists():raise FileExistsError('Fresh twins transfer directory required')
    out.mkdir(parents=True);modelpath=root/'models/chronic_20261003_obesity_discovery_supervised.joblib';digest=sha256(modelpath);bundle=joblib.load(modelpath);model=bundle['model'];ix=model['probes'].selected_indices_;selected=bundle['probe_ids'][ix]
    if len(selected)!=50:raise ValueError('Unchanged primary50 required')
    write_json(out/'protocol.json',{'frozen_utc':datetime.now(timezone.utc).isoformat(),'source':'GSE68336','frozen_model_sha256':digest,'selected_probes':selected.tolist(),
      'unit':'Twin family','expected_samples':70,'expected_pairs':35,'primary_discordant_pairs':24,'concordant_pairs_descriptive':11,
      'endpoint':'Relative higher versus lower BMI within discordant monozygotic twins','absolute_obesity_threshold_available':False,
      'tissue':'Bulk subcutaneous adipose, transferred from isolated adipocyte model','target_training_feature_selection_tuning_calibration':False,'threshold':.5,
      'clinical_obesity_external_validation':False,'source_processed_beta_raw_IDAT_reprocessing':False,
      'detection_p_threshold':.01,'QC_repair_after_V1_values_inspected':True,'V1_kept_as_no_detectionP_sensitivity':True})
    beta,meta,audit=extract(source,selected)
    # Reconstruct only the model's required input positions; the frozen selector
    # discards every other column and imputes absent selected values from source.
    X=np.full((len(meta),len(bundle['probe_ids'])),np.nan,np.float32);X[:,ix]=beta.to_numpy(np.float32)
    with threadpool_limits(limits=1):p=model.predict_proba(X)[:,1]
    meta['probability']=p;meta.to_csv(out/'native_metadata.tsv',sep='\t')
    pd.DataFrame({'sample_id':meta.index,'family_id':meta.family_id,'discordant':meta.discordant,'bmi_group':meta.declared_bmi_group,'probability':p}).to_csv(out/'all_predictions.tsv',sep='\t',index=False)
    primary=meta.loc[meta.discordant];y=(primary.declared_bmi_group=='high').to_numpy(int);metrics=performance(y,primary.probability.to_numpy(),seed=20261003,groups=primary.family_id.to_numpy())
    metrics['interval_method']['discrimination_brier']='Discordant twin-family cluster bootstrap, conditional on fixed predictions'
    pairs=primary.pivot(index='family_id',columns='declared_bmi_group',values='probability');paired_delta=pairs.high-pairs.low;wins=int((paired_delta>0).sum());ties=int((paired_delta==0).sum())
    paired={'pairs':len(pairs),'higher_probability_for_higher_BMI':wins,'ties':ties,'paired_ordering_fraction':wins/len(pairs),
      'ordering_ci95':exact_binomial_interval(wins,len(pairs)),'mean_probability_difference':float(paired_delta.mean())}
    effect=[]
    for probe in selected:
        mapping=primary[['family_id','declared_bmi_group']].copy();mapping['beta']=beta.loc[primary.index,probe]
        cells=mapping.pivot(index='family_id',columns='declared_bmi_group',values='beta');delta=(cells.high-cells.low).dropna();n=len(delta)
        if n>=2 and delta.std(ddof=1)>0:
            se=delta.std(ddof=1)/np.sqrt(n);mean=delta.mean();pvalue=2*stats.t.sf(abs(mean/se),n-1);ci=stats.t.ppf(.975,n-1)*se
            effect.append({'feature_id':probe,'paired_families':n,'delta_beta_high_minus_low':float(mean),'standard_error':float(se),'confidence_low':float(mean-ci),'confidence_high':float(mean+ci),'p_value':float(pvalue)})
        else:effect.append({'feature_id':probe,'paired_families':n,'delta_beta_high_minus_low':float(delta.mean()) if n else None,'standard_error':None,'confidence_low':None,'confidence_high':None,'p_value':None})
    effects=pd.DataFrame(effect);effects['q_value_selected_50_family']=bh(effects.p_value);effects.to_csv(out/'paired_cpg_effects.tsv',sep='\t',index=False)
    annotation=root/'docs/chronic_20261003/obesity_annotation_v1/cpg_gene_interpretation_all50.tsv'
    if annotation.exists():
        context=pd.read_csv(annotation,sep='\t')
        columns=['feature_id','CHR','MAPINFO','coordinate_build','gene_symbols_native','gene_regions_native','cpg_island_context_native','Probe_SNPs','Probe_SNPs_10',
                 'outer_selection_count','outer_selection_denominator','discovery_delta_beta_pp','replication_delta_beta_pp','published_ST7_targets','published_ST9_targets','published_ST11_targets']
        annotated=effects.merge(context[columns],on='feature_id',validate='one_to_one')
        if len(annotated)!=50:raise ValueError('Official primary50 annotation not complete')
        annotated['twin_delta_beta_pp']=annotated.delta_beta_high_minus_low*100
        annotated['direction_matches_source']=np.sign(annotated.twin_delta_beta_pp)==np.sign(annotated.discovery_delta_beta_pp)
        annotated.to_csv(out/'paired_cpg_gene_context.tsv',sep='\t',index=False)
    coverage=pd.DataFrame({'feature_id':selected,'measured_target_samples':beta.notna().sum().to_numpy()});coverage.to_csv(out/'probe_coverage.tsv',sep='\t',index=False)
    if sha256(modelpath)!=digest:raise ValueError('Frozen model changed')
    write_json(out/'source_audit.json',audit)
    write_json(out/'summary.json',{'primary_relative_BMI_transfer_metrics':metrics,'paired_probability_ordering':paired,'target_samples':70,'primary_samples':48,'families':35,'primary_families':24,
      'selected_CpGs_all_missing':int((coverage.measured_target_samples==0).sum()),'frozen_model_unchanged':True,'model_sha256':digest,'source_sha256':sha256(source),
      'clinical_obesity_external_validation':False,'no_target_fitting_or_calibration':True,'versions':versions(),'script_sha256':sha256(Path(__file__)),
      'selected_beta_failed_detection_p_total':int(sum(audit['selected_beta_failed_detection_p_by_sample'].values())),
      'QC_repair_after_V1_values_inspected':True,
      'limits':'Relative BMI phenotype and bulk SAT differ from isolated-adipocyte obese/lean source; public BMI is ordinal; upstream processed beta/QC scope is not rebuilt',
      'generated_utc':datetime.now(timezone.utc).isoformat()})
    print(json.dumps({'relative_BMI_AUROC':metrics['auroc'],'paired_ordering':paired,'missing_selected':int((coverage.measured_target_samples==0).sum())}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.root,a.source,a.output)
