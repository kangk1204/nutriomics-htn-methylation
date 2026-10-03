"""Frozen-model evaluation of a new, separately published SAT-adipocyte cohort.

GSE67024: 15 obese and 14 never-obese women, 450K, PMID26351548.
All three already fitted models are fixed before accessing target beta values.
No feature selection, tuning, calibration or threshold fitting on the target.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(variable,'1')
import joblib
import numpy as np
import pandas as pd
import requests
from threadpoolctl import threadpool_limits
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from nutriomics_methylation.evaluation import performance
from nutriomics_methylation.io import sha256, versions, write_json
from run_chronic_methylation import liver_matrix, checked_beta

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67024/matrix/GSE67024_series_matrix.txt.gz'


def run(root, out):
    if out.exists():
        raise FileExistsError('Use fresh third-cohort output')
    out.mkdir(parents=True)
    models = {
        'original_supervised_50': root / 'models/chronic_20261003_obesity_discovery_supervised.joblib',
        'original_variance': root / 'models/chronic_20261003_obesity_discovery_variance.joblib',
        'snp_excluded_supervised_10': root / 'models/obesity_snp_sensitivity_v1_source_only.joblib'}
    hashes={name:sha256(path) for name,path in models.items()}
    write_json(out/'protocol.json',{'frozen_utc':datetime.now(timezone.utc).isoformat(),
        'primary_model':'original_supervised_50','comparators':list(models)[1:], 'model_sha256':hashes,
        'target':'GSE67024','source':'GSE222595 discovery SAT 48 participants',
        'expected_target_class_counts':[14,15], 'tissue':'isolated abdominal subcutaneous fat cells',
        'scope':'Separately published 2015 women-only cohort; cross-study exploratory transfer',
        'target_beta_or_predictions_previously_inspected':False,
        'source_replication_results_previously_known':True,
        'source_only_frozen_preprocessing':True,'threshold':.5,
        'target_training_feature_selection_tuning_calibration':False,
        'individual_covariate_adjustment':False,'raw_IDAT_reprocessing':False,
        'data_url':URL,'metadata_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67024',
        'paper_pmid':'26351548','primary_model_must_not_be_selected_by_this_test':True})
    data = root/'data/obesity_third_20261003';data.mkdir(exist_ok=True)
    matrix = data/'GSE67024_series_matrix.txt.gz'
    if not matrix.exists():
        temporary = matrix.with_suffix('.partial')
        response=requests.get(URL,stream=True,timeout=(15,90));response.raise_for_status()
        with temporary.open('wb') as stream:
            for block in response.iter_content(1024*1024):stream.write(block)
        temporary.replace(matrix)
    table,meta=liver_matrix(matrix)
    checked_beta(table.to_numpy(dtype=np.float32))
    # Native GEO titles explicitly distinguish PO1 obesity and Control samples.
    title=meta['title']
    label=title.str.contains(r'_PO1(?:\.|$)',regex=True).astype(int)
    control=title.str.contains(r'_Control(?:\.|$)',regex=True)
    if not ((label==1)^control).all() or np.bincount(label).tolist()!=[14,15]:
        raise ValueError('Native obesity/control title mapping differs from prespecified cohort')
    meta['y']=label
    meta.to_csv(out/'native_metadata.tsv',sep='\t')
    y=label.to_numpy(int);rows=[]
    originals=np.load(root/'data/chronic_20261003/prepared_v3/obesity_discovery.npz',allow_pickle=False)
    replication=np.load(root/'data/chronic_20261003/prepared_v3/obesity_replication.npz',allow_pickle=False)
    if set(meta.index)&(set(originals['samples'])|set(replication['samples'])):
        raise ValueError('Cross-study sample accession overlap')
    beta=table.to_numpy(dtype=np.float32).T
    target_lookup={p:i for i,p in enumerate(table.index)}
    with threadpool_limits(limits=1):
        for name,path in models.items():
            bundle=joblib.load(path);fitted=bundle['model'];probes=bundle['probe_ids']
            X=np.full((len(y),len(probes)),np.nan,dtype=np.float32)
            common=[(j,target_lookup[p]) for j,p in enumerate(probes) if p in target_lookup]
            X[:,[j for j,_ in common]]=beta[:,[k for _,k in common]]
            selected=fitted['probes'].selected_indices_
            metrics=performance(y,fitted.predict_proba(X)[:,1],seed=20261003)
            pd.DataFrame({'sample_id':meta.index,'y':y,'probability':fitted.predict_proba(X)[:,1]}).to_csv(out/(name+'_predictions.tsv'),sep='\t',index=False)
            selected_ids=probes[selected]
            coverage=pd.DataFrame({'feature_id':selected_ids,'measured_target_participants':np.isfinite(X[:,selected]).sum(axis=0)})
            coverage.to_csv(out/(name+'_probe_coverage.tsv'),sep='\t',index=False)
            result={'model':name,'target':'GSE67024','metrics':metrics,
                    'selected_probes':len(selected),'selected_probes_all_target_missing':int((coverage.measured_target_participants==0).sum()),
                    'original_model_sha256':hashes[name],'source_only_frozen_model':True,
                    'external_label_or_value_fitting':False,'clinical_validation':False,
                    'interpretation':'Cross-study 29-woman isolated SAT-adipocyte discrimination; raw processing/covariates differ and are not harmonized with test data'}
            write_json(out/(name+'_summary.json'),result);rows.append(result)
            print(json.dumps({'model':name,'auroc':metrics['auroc'],'brier':metrics['brier'],'missing_selected':result['selected_probes_all_target_missing']}),flush=True)
    if hashes!={name:sha256(path) for name,path in models.items()}:
        raise ValueError('A frozen model changed')
    write_json(out/'summary.json',{'results':rows,'source_models_unchanged':True,'target_matrix_sha256':sha256(matrix),
                'target_samples':len(meta),'target_cases':int(y.sum()),'target_controls':int((y==0).sum()),
                'target_probe_rows':len(table),'gzip_full_read_crc_passed':True,'versions':versions(),
                'script_sha256':sha256(Path(__file__)),'generated_utc':datetime.now(timezone.utc).isoformat()})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();run(a.root,a.out)
