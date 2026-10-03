"""Frozen source-only learning and cohort transfer for obesity and fatty liver."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,os,re,sys
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.chronic import nested,tune
from nutriomics_methylation.evaluation import performance,point_metrics
from nutriomics_methylation.io import sha256,versions,write_json

SEEDS=(20261003,20261020,20261103)

def liver_matrix(path):
    headers={}
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for number,line in enumerate(f):
            if line.startswith('!series_matrix_table_begin'):
                skip=number+1;break
            if line.startswith('!Sample_'):
                values=next(csv.reader([line],delimiter='\t'))
                headers.setdefault(values[0],[]).append(values[1:])
        else:raise ValueError('Missing GEO table')
    samples=headers['!Sample_geo_accession'][0]
    meta=pd.DataFrame(index=samples);meta.index.name='sample_id'
    for key,entries in headers.items():
        for j,values in enumerate(entries):
            if len(values)!=len(samples):raise ValueError('Metadata width mismatch')
            column=key.removeprefix('!Sample_')+(f'_{j+1}' if len(entries)>1 else '')
            meta[column]=values
            if key.startswith('!Sample_characteristics'):
                for sample,value in zip(samples,values):
                    if ':' in value:
                        field,val=value.split(':',1)
                        meta.loc[sample,'characteristic_'+re.sub(r'\W+','_',field.strip().lower())]=val.strip()
    types=defaultdict(lambda:np.float32);types['ID_REF']=str
    table=pd.read_csv(path,sep='\t',skiprows=skip,comment='!',index_col=0,dtype=types,na_values=['null','NULL','NA','NaN'])
    if table.columns.tolist()!=samples or table.index.has_duplicates:raise ValueError('GEO matrix alignment or probe duplication')
    return table,meta

def checked_beta(X):
    finite=X[np.isfinite(X)]
    if not finite.size or finite.min()<0 or finite.max()>1:raise ValueError('Expected methylation beta in [0,1]')

def save_prepared(name,table,meta,labels,prepared,out,source):
    selected=labels.notna();m=meta.loc[selected].copy();y=labels[selected].to_numpy(dtype=int)
    if len(set(m.index))!=len(m):raise ValueError('Duplicate sample IDs')
    probes=table.index.to_numpy(dtype=str)
    cg=np.array([p.startswith('cg') for p in probes])
    X=table.loc[cg,m.index].to_numpy(dtype=np.float32).T.copy();probes=probes[cg]
    checked_beta(X)
    if len(np.unique(np.nan_to_num(X,nan=-999),axis=0))!=len(X):raise ValueError('Duplicate within-cohort profiles require group review')
    m['y']=y;m.to_csv(out/(name+'_samples.tsv'),sep='\t')
    np.savez(prepared/(name+'.npz'),X=X,y=y,samples=m.index.to_numpy(dtype=str),probes=probes)
    receipt={'cohort':name,'n':len(y),'controls':int((y==0).sum()),'cases':int(y.sum()),'probes':len(probes),
             'original_samples':len(meta),'excluded_sample_ids':meta.index[~selected].tolist(),
             'sample_id_sha256':hashlib.sha256('\n'.join(m.index).encode()).hexdigest(),
             'source':source,'individual_age_sex_adjustment':False,'raw_IDAT_reprocessing':False,
             'normalization_scope':'Author-provided processed data; upstream normalization scope is not independently recreated.'}
    write_json(out/(name+'_source.json'),receipt)
    print(json.dumps({'prepared':name,'n':len(y),'counts':np.bincount(y).tolist(),'probes':len(probes)}),flush=True)
    return receipt

def prepare_obesity(assets,metadata,prepared,out):
    frames=[pd.read_csv(metadata/(f'GSE222595-{platform}_series_matrix_metadata.tsv'),sep='\t',index_col=0)
            for platform in ['GPL13534','GPL21145']]
    meta=pd.concat(frames)
    ha=meta['description_1'].astype(str)
    if not ha.str.fullmatch(r'HA\d+').all() or ha.duplicated().any():raise ValueError('Native HA-to-GSM map ambiguous')
    mapping=pd.Series(meta.index,index=ha)
    labels=meta['title'].str.extract(r'^adipocyte-(control|obese)-rep\d+$',expand=False).map({'control':0,'obese':1})
    if labels.isna().any():raise ValueError('Obesity phenotype title adapter requires review')
    if set(meta['characteristic_cohort'])!={'discovery','replication'}:raise ValueError('Unknown obesity cohort')
    path=assets/'GSE222595_processed_beta_comb.txt.gz';digest=sha256(path)
    with gzip.open(path,'rt') as f:header=next(csv.reader([next(f)],delimiter='\t'))
    names=header[1::2]
    if len(header)!=385 or header[2::2]!=['Detection Pval']*192 or set(names)!=set(mapping.index):
        raise ValueError('Processed beta/detection P header differs from native mapping')
    pieces=[];probe_chunks=[];detected=np.zeros(192,dtype=np.int64);rows=0
    types=defaultdict(lambda:np.float32);types['ID_REF']=str
    for chunk in pd.read_csv(path,sep='\t',index_col=0,dtype=types,na_values=['NA','null','NULL','NaN'],chunksize=4000):
        if chunk.shape[1]!=384:raise ValueError('Processed numeric width mismatch')
        numeric=chunk.to_numpy();beta=numeric[:,0::2].copy();p=numeric[:,1::2]
        checked_beta(beta)
        if ((p[np.isfinite(p)]<0)|(p[np.isfinite(p)]>1)).any():raise ValueError('Detection P outside [0,1]')
        failed=(~np.isfinite(p))|(p>0.01)|(~np.isfinite(beta));beta[failed]=np.nan
        detected+=(~failed).sum(axis=0);pieces.append(beta);probe_chunks.append(chunk.index.to_numpy(dtype=str));rows+=len(chunk)
        if rows%80000==0:print(json.dumps({'obesity_beta_rows_read':rows}),flush=True)
    probes=np.concatenate(probe_chunks)
    if len(set(probes))!=len(probes):raise ValueError('Duplicate processed probe IDs')
    table=pd.DataFrame(np.concatenate(pieces),index=probes,columns=mapping.loc[names].tolist())
    if sha256(path)!=digest:raise ValueError('Obesity source changed')
    # Availability is a measurement QC property, not an outcome-based feature filter.
    # EPIC-only probes remain missing in discovery and are removed within training folds.
    meta['finite_passed_probes']=pd.Series(detected,index=mapping.loc[names])
    meta['total_union_probes']=rows
    meta.to_csv(out/'obesity_native_sample_mapping.tsv',sep='\t')
    receipts=[]
    for cohort in ['discovery','replication']:
        keep=(meta['characteristic_cohort']==cohort)&(meta['source_name_ch1']=='subcutaneous')
        if np.bincount(labels[keep].to_numpy(dtype=int)).tolist()!=[24,24]:raise ValueError('Declared SAT cohort denominator changed')
        # All deposited SAT samples are retained; original article's 190/192 QC
        # is not reproduced or silently substituted for this declared 48/48 design.
        lab=labels.where(keep)
        receipts.append(save_prepared('obesity_'+cohort,table,meta,lab,prepared,out,
            {'path':str(path),'sha256':digest,'mapping':'Native description_1 HA code; exact header mapping',
             'detection_P_threshold':0.01,'only_subcutaneous':True,'upstream_qc_sample_exclusions_reproduced':False}))
    return receipts

def prepare_liver(assets,prepared,out):
    receipts=[]
    for acc in ['GSE48325','GSE61258']:
        path=assets/(acc+'_series_matrix.txt.gz');digest=sha256(path);table,meta=liver_matrix(path)
        if acc=='GSE48325':
            groups=meta['characteristic_group']
            surgery=[c for c in meta if c.startswith('characteristic_bariatric_surgery')]
            if len(surgery)!=1:raise ValueError('Surgery timing metadata ambiguous')
            timings=meta[surgery[0]].str.strip()
            if not set(timings)<= {'NA','1','2'}:raise ValueError('Unknown surgery timing')
            labels=groups.map({'NORMAL':0,'Healthy Obese':0,'NAFLD':1,'NASH':1}).where(timings!='2')
            exclusions='Exclude all post-bariatric-surgery samples and unknown disease groups; NORMAL/Healthy Obese controls versus NAFLD/NASH.'
        else:
            groups=meta['characteristic_diseasestatus']
            labels=groups.map({'Control':0,'HealthyObese':0,'NAFLD':1,'NASH':1})
            exclusions='Exclude PSC/PBC; Control/Healthy Obese versus NAFLD/NASH. Epigenetic age-derived fields are never predictors.'
        meta['included_primary_nafld']=labels.notna();meta.to_csv(out/(acc+'_all_sample_audit.tsv'),sep='\t')
        receipts.append(save_prepared('nafld_'+acc,table,meta,labels,prepared,out,
                       {'path':str(path),'sha256':digest,'selection':exclusions,'detection_P_data_available_in_this_matrix':False}))
        if sha256(path)!=digest:raise ValueError('Liver source changed')
    return receipts

def align(source,target):
    common=np.array(sorted(set(source['probes'])&set(target['probes'])),dtype=str)
    if not len(common):raise ValueError('No common CpG IDs')
    a={v:i for i,v in enumerate(source['probes'])};b={v:i for i,v in enumerate(target['probes'])}
    return source['X'][:,[a[p] for p in common]],target['X'][:,[b[p] for p in common]],common

def transfer_audit(X,Z,source_ids,target_ids):
    overlap=sorted(set(source_ids)&set(target_ids))
    if overlap:raise ValueError('Cross-cohort sample accessions overlap')
    def signatures(data):return [hashlib.sha256(np.nan_to_num(row,nan=-999).tobytes()).hexdigest() for row in data]
    a=signatures(X);b=signatures(Z)
    duplicate=[str(target_ids[i]) for i,s in enumerate(b) if s in set(a)]
    # Cross-cohort biological donor identity can remain undisclosed even when
    # accessions and processed profiles differ. Do not convert this audit to a
    # claim of independence for liver cohorts with possible related recruitment.
    return {'shared_sample_accessions':overlap,'exact_duplicate_target_profiles':duplicate,
            'all_common_CpGs_used_for_profile_hash':X.shape[1],'unknown_donor_linkage_excluded':False}

def evaluate(source_name,target_name,prepared,out,external_scope):
    source=dict(np.load(prepared/(source_name+'.npz')));target=dict(np.load(prepared/(target_name+'.npz')))
    X,Z,probes=align(source,target);y=source['y'];samples=source['samples'];z=target['y']
    duplicate=transfer_audit(X,Z,samples,target['samples'])
    target_keep=~np.isin(target['samples'],duplicate['exact_duplicate_target_profiles'])
    Z=Z[target_keep];z=z[target_keep];target_samples=target['samples'][target_keep]
    destination=out/source_name;destination.mkdir();write_json(destination/'transfer_overlap_audit.json',duplicate)
    if min(np.bincount(z))<2:raise ValueError('Insufficient target classes after duplicate exclusion')
    combined=[]
    with threadpool_limits(limits=1):
        for ranking in ['supervised','variance']:
            probabilities=[];seed_metrics=[]
            for seed in SEEDS:
                pred,fold,audits=nested(X,y,samples,seed,ranking)
                pd.DataFrame({'sample_id':samples,'y':y,'probability':pred,'fold':fold}).to_csv(destination/f'{ranking}_{seed}_oof.tsv',sep='\t',index=False)
                write_json(destination/f'{ranking}_{seed}_fold_audit.json',audits)
                probabilities.append(pred);m=point_metrics(y,pred);seed_metrics.append({'seed':seed,**m})
                print(json.dumps({'source':source_name,'ranking':ranking,'seed':seed,'nested_auroc':m['auroc']}),flush=True)
            p=np.mean(probabilities,axis=0)
            pd.DataFrame({'sample_id':samples,'y':y,'probability':p}).to_csv(destination/f'{ranking}_repeated_oof.tsv',sep='\t',index=False)
            fitted,tuning=tune(X,y,SEEDS[0],ranking)
            # No fit/tuning/threshold selection/calibration on target values or outcomes.
            external=fitted.predict_proba(Z)[:,1]
            pd.DataFrame({'sample_id':target_samples,'y':z,'probability':external}).to_csv(destination/f'{ranking}_external.tsv',sep='\t',index=False)
            write_json(destination/f'{ranking}_source_only_tuning.json',tuning)
            ix=fitted['probes'].selected_indices_
            pd.DataFrame({'feature_id':probes[ix],'standardized_coefficient':fitted['model'].coef_[0]}).to_csv(destination/f'{ranking}_refit_features.tsv',sep='\t',index=False)
            npz=Path(__file__).resolve().parents[1]/'models'/('chronic_20261003_'+source_name+'_'+ranking+'.joblib')
            npz.parent.mkdir(exist_ok=True)
            if npz.exists():raise FileExistsError('Model artifact already exists')
            import joblib
            joblib.dump({'model':fitted,'probe_ids':probes,'training_samples':samples,'source_cohort':source_name},npz)
            result={'source':source_name,'target':target_name,'modality':'DNAm','ranking':ranking,
                    'common_probe_ids':len(probes),'repeated_nested_oof_metrics':performance(y,p,seed=SEEDS[0]),
                    'seed_metrics':seed_metrics,'training_resubstitution_metrics':point_metrics(y,fitted.predict_proba(X)[:,1]),
                    'source_only_frozen_external_metrics':performance(z,external,seed=SEEDS[0]),
                    'external_scope':external_scope,'target_excluded_exact_duplicates':duplicate['exact_duplicate_target_profiles'],
                    'target_fitting_tuning_calibration_or_threshold_selection':False,
                    'external_cohort_raw_processing_independence_established':False,
                    'individual_age_sex_adjustment_performed':False,'clinical_validation':False,
                    'model_artifact':{'path':str(npz),'sha256':sha256(npz)},
                    'interpretation':'Exploratory tissue-specific discrimination; source cohort CV and frozen cohort transfer are reported separately.'}
            write_json(destination/(ranking+'_summary.json'),result);combined.append(result)
            print(json.dumps({'source':source_name,'ranking':ranking,'external_auroc':result['source_only_frozen_external_metrics']['auroc']}),flush=True)
    return combined

def main():
    p=argparse.ArgumentParser();p.add_argument('--assets',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--prepared',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists() or a.prepared.exists():raise FileExistsError('New output/prepared directories required')
    a.out.mkdir(parents=True);a.prepared.mkdir(parents=True)
    write_json(a.out/'plan.json',{'frozen_before_these_cohort_model_results':True,'seeds':list(SEEDS),
      'source_target_pairs':[['obesity_discovery','obesity_replication'],['nafld_GSE48325','nafld_GSE61258']],
      'obesity_design':'SAT only; discovery 450K 24 controls/24 obese, held-out replication EPIC 24 controls/24 obese; exact native HA mapping; mask Detection P>0.01; no target training.',
      'liver_design':'Exclude all postoperative and unknown GSE48325 samples; healthy obese/normal controls versus NAFLD/NASH. GSE61258 Control/Healthy Obese versus NAFLD/NASH; PSC/PBC excluded.',
      'liver_external_claim':'Exploratory cohort transfer; donor independence not established; related recruitment and upstream normalization may limit inference.',
      'source_CV':'Nested 5-fold/minority-limited, inner 3-fold/minority-limited; every learned feature, imputation and scaler training-only.',
      'hyperparameters':{'k':[10,50,200],'C':[0.1,1.0],'ranking':['supervised','variance'],'threshold':0.5},
      'all_results_reported':True,'raw_source_matrices_not_redistributed':True})
    receipts=prepare_obesity(a.assets,a.metadata,a.prepared,a.out)+prepare_liver(a.assets,a.prepared,a.out)
    write_json(a.out/'preparation_receipt.json',receipts)
    plans=[('obesity_discovery','obesity_replication','Separate discovery/replication cohorts declared by original study; author-processed combined beta matrix upstream pipeline not independently reconstructed.'),
           ('nafld_GSE48325','nafld_GSE61258','Exploratory liver cohort transfer; source papers and related recruitment do not establish independent donor identities.')]
    groups=Parallel(n_jobs=2)(delayed(evaluate)(src,tgt,a.prepared,a.out,scope) for src,tgt,scope in plans)
    results=[r for group in groups for r in group];write_json(a.out/'summary.json',results)
    rows=[{'source':r['source'],'target':r['target'],'ranking':r['ranking'],'train_n':r['repeated_nested_oof_metrics']['n'],
      'source_nested_auroc':r['repeated_nested_oof_metrics']['auroc'],'external_n':r['source_only_frozen_external_metrics']['n'],
      'external_auroc':r['source_only_frozen_external_metrics']['auroc'],
      'external_ci_lower':r['source_only_frozen_external_metrics']['ci95']['auroc'][0],
      'external_ci_upper':r['source_only_frozen_external_metrics']['ci95']['auroc'][1],
      'training_auroc':r['training_resubstitution_metrics']['auroc'],'claim':r['external_scope']} for r in results]
    pd.DataFrame(rows).to_csv(a.out/'comparison.tsv',sep='\t',index=False)
    write_json(a.out/'receipt.json',{'generated_utc':datetime.now(timezone.utc).isoformat(),'versions':versions(),
      'script_sha256':sha256(Path(__file__)),'learning_code_sha256':sha256(Path(__file__).resolve().parents[1]/'src/nutriomics_methylation/chronic.py'),
      'models_reported':len(results),'nested_prediction_runs':len(results)*len(SEEDS),'all_evaluated_models_reported':True})

if __name__=='__main__':main()
