"""Fresh HTN grid; all models and every source OOF prediction are retained."""
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
from nutriomics_methylation.extension_grid import FAMILIES,FEATURE_COUNTS,REGULARIZATION,nested_all,tune_all,calibrated_fit
from nutriomics_methylation.evaluation import performance,match_groups
from nutriomics_methylation.io import sha256,versions,write_json

def run(root,out):
    if out.exists():raise FileExistsError('Fresh output required')
    out.mkdir(parents=True);path=root/'data/prepared/matrices.npz';digest=sha256(path);data=np.load(path,allow_pickle=False)
    keep=data['primary_phenotypes']!='preHT';y=(data['primary_phenotypes'][keep]=='HTN').astype(int);samples=data['primary_samples'][keep]
    if np.bincount(y).tolist()!=[44,44] or (~keep).sum()!=44:raise ValueError('Primary HTN denominator changed')
    seeds=[20261003,20261020,20261103]
    write_json(out/'protocol.json',{'frozen_utc':datetime.now(timezone.utc).isoformat(),'source':'GSE193795','n':88,'preHT_excluded':44,
      'primary_endpoint':'Depositor HTN versus control','families':FAMILIES,'feature_counts':FEATURE_COUNTS,'C':REGULARIZATION,'elasticnet_l1_ratio':.5,
      'nested_outer_folds':5,'inner_folds':3,'calibration_folds':3,'seeds':seeds,'selection_metric':'Internal AUROC; no target selection',
      'spaces':['common','full'],'threshold':.5,'external':'GSE42774; previously seen diagnostic only',
      'all_models_reported':True,'prepared_sha256':digest,'individual_age_sex_adjustment':False})
    emeta=pd.read_csv(root/'data/prepared/external_metadata.tsv',sep='\t',index_col=0);egroups=match_groups(emeta);rows=[]
    with threadpool_limits(limits=1):
        for space in ('common','full'):
            dest=out/space;dest.mkdir();probes=data['common_probes'] if space=='common' else data['probe_ids']
            ix=pd.Index(data['probe_ids']).get_indexer(probes)
            if (ix<0).any():raise ValueError('Probe alignment failed')
            X=data['primary'][keep][:,ix].astype(np.float32);allpred={f:[] for f in FAMILIES}
            for seed in seeds:
                results=nested_all(X,y,samples,seed)
                for family,item in results.items():
                    pd.DataFrame({'sample_id':samples,'y':y,'probability':item['probability'],'fold':item['fold']}).to_csv(dest/f'{family}_{seed}_oof.tsv',sep='\t',index=False)
                    write_json(dest/f'{family}_{seed}_fold_audit.json',item['audits']);allpred[family].append(item['probability'])
            choices=tune_all(X,y,seeds[0]);write_json(dest/'source_only_selection.json',choices)
            for family in FAMILIES:
                probability=np.mean(allpred[family],axis=0);metrics=performance(y,probability,seed=seeds[0])
                fitted=calibrated_fit(X,y,choices[family],seeds[0]);joblib.dump({'model':fitted,'probe_ids':probes,'training_samples':samples},dest/f'{family}_source_only.joblib')
                pd.DataFrame({'sample_id':samples,'y':y,'probability':probability}).to_csv(dest/f'{family}_repeated_oof.tsv',sep='\t',index=False)
                write_json(dest/f'{family}_calibration_audit.json',fitted.audit)
                pd.DataFrame({'feature_id':probes[fitted.base['probes'].selected_indices_],'standardized_coefficient':fitted.base['model'].coef_[0]}).to_csv(dest/f'{family}_features.tsv',sep='\t',index=False)
                external=None
                if space=='common':
                    p=fitted.predict_proba(data['external'])[:,1]
                    pd.DataFrame({'sample_id':data['external_samples'],'y':data['external_y'],'probability':p}).to_csv(dest/f'{family}_seen_external.tsv',sep='\t',index=False)
                    external=performance(data['external_y'],p,seed=seeds[0],groups=egroups)
                result={'space':space,'family':family,'source_nested':metrics,'seen_external_diagnostic':external,'selected_k':choices[family]['k'],'selected_C':choices[family]['C'],
                  'external_model_selection':False,'external_calibration':False,'threshold':.5,'clinical_validation':False}
                write_json(dest/f'{family}_summary.json',result);rows.append(result);print(f'{space}/{family}: source AUROC={metrics["auroc"]:.6f}',flush=True)
    if sha256(path)!=digest:raise ValueError('Original input changed')
    write_json(out/'summary.json',{'analyses':rows,'generated_utc':datetime.now(timezone.utc).isoformat(),'versions':versions(),'script_sha256':sha256(Path(__file__)),'code_sha256':sha256(Path(__file__).resolve().parents[1]/'src/nutriomics_methylation/extension_grid.py'),'original_matrix_unchanged':True})
    pd.DataFrame([{'space':r['space'],'family':r['family'],'source_auroc':r['source_nested']['auroc'],'source_ci_low':r['source_nested']['ci95']['auroc'][0],'source_ci_high':r['source_nested']['ci95']['auroc'][1],'source_brier':r['source_nested']['brier'],'diagnostic_external_auroc':None if r['seen_external_diagnostic'] is None else r['seen_external_diagnostic']['auroc']} for r in rows]).to_csv(out/'comparison.tsv',sep='\t',index=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.root,a.output)
