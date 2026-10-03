"""Independently regenerate source OOF and apply the unchanged primary50 model.

This is reproducibility verification on already inspected cohorts. It does not
create a new confirmation cohort, retune the frozen model or refit its target
preprocessing. Source OOF fits are new independent fits using the original
prespecified grid; the original serialized model is only used for prediction.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import os
from pathlib import Path
import sys
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(v,'1')
import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.chronic import nested
from nutriomics_methylation.evaluation import performance
from nutriomics_methylation.io import sha256,versions,write_json
from run_chronic_methylation import liver_matrix,checked_beta,SEEDS

def align_bundle(data,probes):
    ix=pd.Index(data['probes']).get_indexer(probes);X=np.full((len(data['y']),len(probes)),np.nan,np.float32)
    present=ix>=0;X[:,present]=data['X'][:,ix[present]]
    return X

def compare_predictions(path,ids,p):
    if not path.exists():return {'previous_prediction_file_present':False,'path':str(path)}
    old=pd.read_csv(path,sep='\t').set_index('sample_id')
    if set(old.index)!=set(ids):raise ValueError('Replay participant set differs')
    delta=float(np.max(np.abs(old.loc[ids,'probability'].to_numpy()-p)))
    return {'previous_prediction_file_present':True,'path':str(path),'sha256':sha256(path),'max_absolute_probability_difference':delta,'within_1e_minus_10':delta<1e-10}

def run(root,out):
    if out.exists():raise FileExistsError('Fresh replay directory required')
    out.mkdir(parents=True);modelpath=root/'models/chronic_20261003_obesity_discovery_supervised.joblib';mh=sha256(modelpath)
    frozen=joblib.load(modelpath);probes=frozen['probe_ids'];model=frozen['model']
    if len(model['probes'].selected_indices_)!=50:raise ValueError('Primary50 immutable artifact differs')
    prior=root/'docs/chronic_20261003/methylation_transfer_v3/obesity_discovery'
    write_json(out/'protocol.json',{'frozen_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Independent computation reproduction of previously inspected results',
      'source':'GSE222595 SAT discovery','old_replication':'GSE222595 SAT EPIC replication','third':'GSE67024',
      'primary_model_sha256':mh,'frozen_model_fit_called':False,'target_fitting_tuning_calibration':False,'threshold':.5,'source_nested_seeds':SEEDS,
      'source_grid':{'k':[10,50,200],'C':[.1,1.]},'independent_raw_IDAT_reprocessing':False})
    source=dict(np.load(root/'data/chronic_20261003/prepared_v3/obesity_discovery.npz',allow_pickle=False))
    replicate=dict(np.load(root/'data/chronic_20261003/prepared_v3/obesity_replication.npz',allow_pickle=False))
    X=align_bundle(source,probes);Z=align_bundle(replicate,probes);rows=[];comparison={};preds=[]
    with threadpool_limits(limits=1):
        for seed in SEEDS:
            pred,fold,audit=nested(X,source['y'],source['samples'],seed,'supervised');preds.append(pred)
            pd.DataFrame({'sample_id':source['samples'],'y':source['y'],'probability':pred,'fold':fold}).to_csv(out/f'source_{seed}_oof.tsv',sep='\t',index=False)
            write_json(out/f'source_{seed}_audit.json',audit)
            comparison[f'source_seed_{seed}']=compare_predictions(prior/f'supervised_{seed}_oof.tsv',source['samples'],pred)
        p=np.mean(preds,axis=0)
        pd.DataFrame({'sample_id':source['samples'],'y':source['y'],'probability':p}).to_csv(out/'source_repeated_oof.tsv',sep='\t',index=False)
        rows.append({'cohort':'source_nested','metrics':performance(source['y'],p,seed=SEEDS[0])})
        comparison['source_repeated']=compare_predictions(prior/'supervised_repeated_oof.tsv',source['samples'],p)
        for name,data,matrix,old in [('source_resubstitution',source,X,None),('old_epic_replication',replicate,Z,prior/'supervised_external.tsv')]:
            p=model.predict_proba(matrix)[:,1];pd.DataFrame({'sample_id':data['samples'],'y':data['y'],'probability':p}).to_csv(out/(name+'_predictions.tsv'),sep='\t',index=False)
            rows.append({'cohort':name,'metrics':performance(data['y'],p,seed=SEEDS[0])})
            if old:comparison[name]=compare_predictions(old,data['samples'],p)
        table,meta=liver_matrix(root/'data/obesity_third_20261003/GSE67024_series_matrix.txt.gz');checked_beta(table.to_numpy(np.float32))
        title=meta['title'];label=title.str.contains(r'_PO1(?:\.|$)',regex=True).astype(int);control=title.str.contains(r'_Control(?:\.|$)',regex=True)
        if not ((label==1)^control).all() or np.bincount(label).tolist()!=[14,15]:raise ValueError('GSE67024 native phenotype mapping changed')
        third={'probes':table.index.to_numpy(str),'X':table.to_numpy(np.float32).T,'y':label.to_numpy(int),'samples':meta.index.to_numpy(str)}
        thirdX=align_bundle(third,probes);p=model.predict_proba(thirdX)[:,1]
        pd.DataFrame({'sample_id':third['samples'],'y':third['y'],'probability':p}).to_csv(out/'gse67024_predictions.tsv',sep='\t',index=False)
        rows.append({'cohort':'GSE67024_seen_replay','metrics':performance(third['y'],p,seed=SEEDS[0])})
        comparison['GSE67024']=compare_predictions(root/'docs/obesity_third_20261003_v1/original_supervised_50_predictions.tsv',third['samples'],p)
    if sha256(modelpath)!=mh:raise ValueError('Immutable model changed during replay')
    write_json(out/'summary.json',{'results':rows,'comparisons':comparison,'all_existing_predictions_reproduced':all(v.get('within_1e_minus_10',False) for v in comparison.values()),
      'primary_model_sha256_before_after':mh,'frozen_model_unchanged':True,'generated_utc':datetime.now(timezone.utc).isoformat(),'versions':versions(),'script_sha256':sha256(Path(__file__)),
      'scope':'Reexecution on existing processed beta; not raw IDAT reconstruction or unused-cohort validation'})
    pd.DataFrame([{'cohort':r['cohort'],'n':r['metrics']['n'],'auroc':r['metrics']['auroc'],'auroc_ci_low':r['metrics']['ci95']['auroc'][0],'auroc_ci_high':r['metrics']['ci95']['auroc'][1],'brier':r['metrics']['brier']} for r in rows]).to_csv(out/'comparison.tsv',sep='\t',index=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.root,a.output)
