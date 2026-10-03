"""Center-held-out endocrine versus primary HTN from repository raw beta.

All four endpoints are declared before modeling. No healthy controls exist in
this study. Repository raw beta, not global quantile/ComBat normalized data, is
used. Age/sex-only predictions are a prespecified confounding comparator.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import sys
import re
import shutil
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(v,'1')
import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.extension_grid import FAMILIES,FEATURE_COUNTS,REGULARIZATION,nested_all,tune_all,calibrated_fit,splits
from nutriomics_methylation.evaluation import performance
from nutriomics_methylation.io import sha256,versions,write_json

EXPECTED={'PHT':42,'CS':57,'PA':101,'PPGL':55}

def prepare(raw,sheet,dest):
    meta=pd.read_excel(sheet).set_index('Sample')
    if meta.index.has_duplicates or meta.Hypertension_type.value_counts().to_dict()!=EXPECTED:raise ValueError('Native endocrine phenotype counts differ')
    if meta['Enrolling center'].isna().any() or not set(meta.Sex)<={'M','F'} or meta.Age.isna().any():raise ValueError('Missing required published center/age/sex metadata')
    meta.index=meta.index.astype(str)
    with raw.open(encoding='utf-8') as stream:
        header=next(csv.reader([stream.readline()],delimiter='\t'));first=next(csv.reader([stream.readline()],delimiter='\t'))
    # Official BetaValue_Raw.txt is an R write.table export: the 255-column
    # header omits the CpG row-name heading, while each data row has 256 fields.
    implicit_row_names=len(first)==len(header)+1
    columns=header if implicit_row_names else header[1:]
    if len(first)!=256 or len(columns)!=255:raise ValueError('Unexpected official raw-beta header width')
    resolved=[re.sub(r'^S_(\d+)$',r'S\1',c) for c in columns]
    if set(resolved)!=set(meta.index) or len(set(resolved))!=255:raise ValueError(f'Raw-beta columns cannot be mapped exactly to supplementary participant IDs: {columns[:5]}')
    if len(set(columns))!=len(columns):raise ValueError('Repeated raw-beta participant')
    records=[];blocks=[];total=0
    read_options={'names':['feature_id',*columns],'header':None,'skiprows':1,'index_col':0,'sep':'\t','chunksize':4000,
      'dtype':{c:np.float32 for c in columns},'na_values':['NA','NaN','null','NULL']}
    for chunk in pd.read_csv(raw,**read_options):
        if chunk.shape[1]!=255 or chunk.columns.tolist()!=columns:raise ValueError('Raw-beta width or sample order mismatch')
        probes=chunk.index.astype(str);keep=probes.str.startswith('cg');values=chunk.loc[keep].to_numpy(np.float32)
        finite=values[np.isfinite(values)]
        if len(finite) and (finite.min()<0 or finite.max()>1):raise ValueError('Raw methylation beta outside [0,1]')
        records.extend(probes[keep]);blocks.append(values);total+=len(chunk)
        if total%80000==0:print(json.dumps({'raw_beta_rows_read':total}),flush=True)
    if len(set(records))!=len(records):raise ValueError('Repeated raw CpG IDs')
    X=np.concatenate(blocks,axis=0).T.copy();meta=meta.loc[resolved];samples=np.array(resolved,str);probes=np.array(records,str)
    np.savez(dest/'raw_prepared.npz',X=X,samples=samples,probes=probes)
    meta.to_csv(dest/'native_metadata.tsv',sep='\t')
    write_json(dest/'preparation.json',{'n':255,'raw_rows':total,'cg_rows':len(probes),'raw_beta_sha256':sha256(raw),'metadata_sha256':sha256(sheet),
      'all_source_beta_rows_read':True,'global_normalization_or_ComBat':False,'native_label_counts':EXPECTED,
      'center_subtype_table':pd.crosstab(meta['Enrolling center'],meta.Hypertension_type).to_dict(),
      'raw_header_implicit_R_row_names':implicit_row_names,'native_to_paper_sample_ID_rule':'Exact S_<positive integer> to S<same integer>; all255 IDs form a one-to-one set match',
      'native_header_sample_ids':columns,'paper_sample_ids_in_numeric_order':resolved,
      'raw_IDAT_or_detectionP_available_here':False,'normalization_scope':'Repository unnormalized beta; no new cohort-wide learned normalization applied'})
    return X,meta,samples,probes

def clustered_metrics(y,p,groups):
    m=performance(y,p,seed=20261003,groups=groups)
    m['interval_method']['discrimination_brier']='Recruitment-center cluster bootstrap, conditional on fixed OOF predictions'
    m['interval_method']['sensitivity_specificity']='Clopper-Pearson participant intervals; center clustering is not represented by these intervals'
    return m

def run(raw,sheet,out):
    if out.exists():raise FileExistsError('Fresh endocrine output required')
    if raw.name!='BetaValue_Raw.txt':raise ValueError('Only repository raw beta admitted; normalized/ComBat data excluded')
    out.mkdir(parents=True)
    code=out/'executed_code';code.mkdir()
    module=Path(__file__).resolve().parents[1]/'src/nutriomics_methylation/extension_grid.py'
    shutil.copy2(__file__,code/'run_endocrine_extension.py');shutil.copy2(module,code/'extension_grid.py')
    write_json(out/'protocol.json',{'frozen_utc':datetime.now(timezone.utc).isoformat(),'source':'S-BSST831','paper_pmid':'36329530',
      'contrasts':['CS_vs_PHT','PA_vs_PHT','PPGL_vs_PHT','EHT_vs_PHT'],'unit':'Participant, recruitment center held out',
      'outer_folds':'3 (only 3 PHT recruitment centers)','inner_and_calibration_folds':'2 within outer training centers',
      'group_partition':'Exhaustive feasible minimum label-count imbalance for <=10 centers; no features/predictions enter fold selection',
      'families':FAMILIES,'k':FEATURE_COUNTS,'C':REGULARIZATION,'elasticnet_l1_ratio':.5,'seed':20261003,'threshold':.5,
      'healthy_controls':0,'all_models_and_endpoints_reported':True,'primary_vs_normotensive_HTN_claim':False,
      'global_quantile_normalized_or_ComBat_beta_used':False,'confounding_comparator':'Fixed C1 ridge, age and sex only, identical outer center folds',
      'predeclared_source_raw':str(raw),'metadata':str(sheet),'script_sha256':sha256(Path(__file__)),
      'extension_grid_sha256':sha256(module),'metadata_sha256':sha256(sheet)})
    X,meta,samples,probes=prepare(raw,sheet,out);rows=[]
    with threadpool_limits(limits=1):
        for subtype in ('CS','PA','PPGL','EHT'):
            keep=(meta.Hypertension_type=='PHT')|((meta.Hypertension_type!='PHT') if subtype=='EHT' else (meta.Hypertension_type==subtype))
            m=meta.loc[keep];xx=X[keep.to_numpy()];y=(m.Hypertension_type!='PHT').to_numpy(int);ids=samples[keep.to_numpy()];groups=m['Enrolling center'].to_numpy(int)
            dest=out/(subtype+'_vs_PHT');dest.mkdir();results=nested_all(xx,y,ids,20261003,groups)
            choices=tune_all(xx,y,20261003,groups);write_json(dest/'source_only_selection.json',choices)
            for family,item in results.items():
                pd.DataFrame({'sample_id':ids,'center':groups,'y':y,'probability':item['probability'],'fold':item['fold']}).to_csv(dest/(family+'_oof.tsv'),sep='\t',index=False)
                write_json(dest/(family+'_fold_audit.json'),item['audits']);metrics=clustered_metrics(y,item['probability'],groups)
                fitted=calibrated_fit(xx,y,choices[family],20261003,groups);joblib.dump({'model':fitted,'probe_ids':probes,'training_samples':ids},dest/(family+'_source_only.joblib'))
                pd.DataFrame({'feature_id':probes[fitted.base['probes'].selected_indices_],'standardized_coefficient':fitted.base['model'].coef_[0]}).to_csv(dest/(family+'_features.tsv'),sep='\t',index=False)
                result={'endpoint':subtype+'_vs_PHT','family':family,'nested_center_metrics':metrics,'selected_k':choices[family]['k'],'selected_C':choices[family]['C'],
                  'healthy_control_classification':False,'independent_external_study':False,'individual_age_sex_adjustment':False,
                  'limitations':'Center/subtype strongly confounded; raw beta has no IDAT/detectionP reprocessing; OOF probabilities are research scores'}
                write_json(dest/(family+'_summary.json'),result);rows.append(result);print(json.dumps({'endpoint':result['endpoint'],'family':family,'AUROC':metrics['auroc']}),flush=True)
            cov=np.column_stack([m.Age.to_numpy(float),(m.Sex=='M').to_numpy(float)]);p=np.full(len(y),np.nan);fold=np.full(len(y),-1)
            for index,(tr,te) in enumerate(splits(y,20261003,5,groups)):
                baseline=Pipeline([('scale',StandardScaler()),('model',LogisticRegression(C=1.,solver='liblinear',class_weight='balanced',random_state=20261003))]).fit(cov[tr],y[tr])
                p[te]=baseline.predict_proba(cov[te])[:,1];fold[te]=index
            pd.DataFrame({'sample_id':ids,'center':groups,'y':y,'probability':p,'fold':fold}).to_csv(dest/'age_sex_only_oof.tsv',sep='\t',index=False)
            metrics=clustered_metrics(y,p,groups);result={'endpoint':subtype+'_vs_PHT','family':'age_sex_only','nested_center_metrics':metrics,'selected_k':2,'selected_C':1.,'not_methylation_model':True}
            write_json(dest/'age_sex_only_summary.json',result);rows.append(result)
    write_json(out/'summary.json',{'analyses':rows,'versions':versions(),'generated_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha256(Path(__file__)),
      'extension_grid_sha256':sha256(module),'all_four_endpoints_reported':True,'healthy_control_count':0,'independent_external_cohort':False})
    pd.DataFrame([{'endpoint':r['endpoint'],'family':r['family'],'n':r['nested_center_metrics']['n'],'auroc':r['nested_center_metrics']['auroc'],'ci95_low':r['nested_center_metrics']['ci95']['auroc'][0],'ci95_high':r['nested_center_metrics']['ci95']['auroc'][1],'brier':r['nested_center_metrics']['brier']} for r in rows]).to_csv(out/'comparison.tsv',sep='\t',index=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.raw,a.metadata,a.output)
