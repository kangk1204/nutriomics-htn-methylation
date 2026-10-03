"""Prespecified nested screening on three verified public disease cohorts."""
from __future__ import annotations
import argparse,json,os,sys
from datetime import datetime,timezone
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.chronic import nested,tune
from nutriomics_methylation.io import read_matrix,sha256,versions,write_json
from nutriomics_methylation.evaluation import performance,point_metrics

COHORTS={
 'GSE89632':{'disease':'NAFLD','modality':'expression','tissue':'liver','field':'characteristic_diagnosis','labels':{'HC':0,'SS':1,'NASH':1},'expected':[24,39]},
 'GSE38642':{'disease':'T2D','modality':'expression','tissue':'pancreatic islets','field':'characteristic_status','labels':{'non-diabetic donors':0,'T2D donors':1},'expected':[54,9]},
 'GSE21232':{'disease':'T2D','modality':'DNAm','tissue':'pancreatic islets','field':'characteristic_disease_state','labels':{'non-diabetic':0,'type 2 diabetes':1},'expected':[11,5]},
}
SEEDS=(20261003,20261020,20261103)

def one(acc,raw,out):
    config=COHORTS[acc];path=raw/(acc+'_series_matrix.txt.gz');digest=sha256(path)
    table,meta=read_matrix(path)
    mapped=meta[config['field']].map(config['labels'])
    if mapped.isna().any():raise ValueError('Unknown disease labels; adapter review required')
    y=mapped.to_numpy(dtype=int);samples=meta.index.to_numpy(dtype=str);X=table.to_numpy(dtype=np.float32).T
    if np.bincount(y).tolist()!=config['expected']:raise ValueError('Declared disease denominator changed')
    if len(set(samples))!=len(samples):raise ValueError('Duplicate participant accession')
    if config['modality']=='DNAm':
        finite=X[np.isfinite(X)]
        if not len(finite) or finite.min()<0 or finite.max()>1:raise ValueError('Expected processed beta values in [0,1]')
    if len(np.unique(np.nan_to_num(X,nan=-999),axis=0))!=len(X):raise ValueError('Exact duplicate profiles require grouped processing')
    target=out/acc;target.mkdir()
    meta['y']=y;meta.to_csv(target/'sample_audit.tsv',sep='\t')
    source={'accession':acc,'source_file':str(path),'sha256':digest,'n':len(y),'n_positive':int(y.sum()),'features':X.shape[1],**config,
            'public_processed_data':True,'normalization_scope':'Author-supplied processed cohort; original cohort normalization cannot be undone here.',
            'individual_covariate_adjustment_performed':False,'independent_external_evaluation':False,'tissue_representative_of_blood':False}
    write_json(target/'source.json',source)
    summaries=[]
    with threadpool_limits(limits=1):
        for ranking in ('supervised','variance'):
            probabilities=[];individual=[]
            for seed in SEEDS:
                probability,fold,audits=nested(X,y,samples,seed,ranking)
                probabilities.append(probability)
                pd.DataFrame({'sample_id':samples,'y':y,'probability':probability,'fold':fold}).to_csv(target/f'{ranking}_{seed}_oof.tsv',sep='\t',index=False)
                write_json(target/f'{ranking}_{seed}_fold_audit.json',audits)
                m=point_metrics(y,probability);individual.append({'seed':seed,**m})
                print(json.dumps({'accession':acc,'ranking':ranking,'seed':seed,'nested_auroc':m['auroc']}),flush=True)
            averaged=np.mean(probabilities,axis=0)
            pd.DataFrame({'sample_id':samples,'y':y,'probability':averaged}).to_csv(target/f'{ranking}_repeated_oof.tsv',sep='\t',index=False)
            m=performance(y,averaged,seed=SEEDS[0],bootstraps=1000)
            fitted,audit=tune(X,y,SEEDS[0],ranking)
            apparent=point_metrics(y,fitted.predict_proba(X)[:,1])
            selection=fitted['probes'].selected_indices_
            pd.DataFrame({'feature_id':table.index.to_numpy()[selection],'standardized_coefficient':fitted['model'].coef_[0],
                          'status':'exploratory full-source refit; not independently validated marker'}).to_csv(target/f'{ranking}_refit_features.tsv',sep='\t',index=False)
            write_json(target/f'{ranking}_final_tuning.json',audit)
            summary={'accession':acc,'disease':config['disease'],'modality':config['modality'],'tissue':config['tissue'],'ranking':ranking,
                     'repeated_nested_oof_metrics':m,'seed_metrics':individual,'training_resubstitution_metrics':apparent,
                     'independent_external_evaluation':False,'interpretation':'Exploratory tissue-specific discrimination; not dietary effects, causal markers, clinical risk or externally validated improvement.'}
            write_json(target/f'{ranking}_summary.json',summary);summaries.append(summary)
    if sha256(path)!=digest:raise ValueError('Source changed during analysis')
    return summaries

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cores',type=int,default=3);a=p.parse_args()
    if a.out.exists():raise FileExistsError('Preserve all earlier analyses; new output required')
    if not 1<=a.cores<=3:raise ValueError('At most three cohort processes')
    a.out.mkdir(parents=True)
    plan={'fixed_before_chronic_model_evaluation':True,'cohorts':COHORTS,'seeds':list(SEEDS),'outer':'Stratified min(5,minority n)','inner':'Stratified min(3,training minority n)',
          'selection':['supervised F statistic inside folds','variance inside folds'],'k':[10,50,200],'C':[0.1,1.0],'threshold':0.5,
          'class_weight':'balanced','primary':'Repeated nested OOF; training score separately','negative_results_reported':True,'multiple_cohort_screening':'Exploratory; no confirmatory selected-positive claim'}
    write_json(a.out/'plan.json',plan)
    results=Parallel(n_jobs=a.cores)(delayed(one)(acc,a.raw,a.out) for acc in COHORTS)
    combined=[r for group in results for r in group]
    write_json(a.out/'summary.json',combined)
    rows=[{'accession':r['accession'],'disease':r['disease'],'modality':r['modality'],'ranking':r['ranking'],
           'n':r['repeated_nested_oof_metrics']['n'],'n_positive':r['repeated_nested_oof_metrics']['n_positive'],
           'nested_oof_auroc':r['repeated_nested_oof_metrics']['auroc'],'ci_lower':r['repeated_nested_oof_metrics']['ci95']['auroc'][0],
           'ci_upper':r['repeated_nested_oof_metrics']['ci95']['auroc'][1],'brier':r['repeated_nested_oof_metrics']['brier'],
           'training_auroc':r['training_resubstitution_metrics']['auroc'],'external':'not evaluated'} for r in combined]
    pd.DataFrame(rows).to_csv(a.out/'comparison.tsv',sep='\t',index=False)
    write_json(a.out/'receipt.json',{'generated_utc':datetime.now(timezone.utc).isoformat(),'versions':versions(),'script_sha256':sha256(Path(__file__)),
               'learning_code_sha256':sha256(Path(__file__).resolve().parents[1]/'src/nutriomics_methylation/chronic.py'),
               'models_evaluated':len(combined),'nested_prediction_runs':len(combined)*len(SEEDS),'raw_data_not_redistributed':True,
               'all_evaluated_models_reported':True,'no_person_matched_multiomics_merge':True})

if __name__=='__main__':main()
