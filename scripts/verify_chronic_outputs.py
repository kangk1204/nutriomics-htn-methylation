"""Recalculate exported scores, denominators and inner/outer leakage audits."""
import argparse,json,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.evaluation import point_metrics
from nutriomics_methylation.io import sha256,write_json

SEEDS=(20261003,20261020,20261103)

def verify_metric(path,expected,probability_dtype=None):
    d=pd.read_csv(path,sep='\t')
    assert d.sample_id.is_unique
    probability=d.probability.to_numpy(dtype=probability_dtype)
    m=point_metrics(d.y.to_numpy(),probability)
    for key in ['n','n_positive','auroc','average_precision','brier','sensitivity','specificity','tp','tn','fp','fn']:
        assert np.isclose(m[key],expected[key],atol=1e-12,rtol=0),(path,key,m[key],expected[key])
    return {'path':str(path),'sha256':sha256(path),'n':len(d),'auroc':m['auroc'],
            'probability_serialization_dtype':'float32 native saga output round-trip' if probability_dtype==np.float32 else 'float64'}

def verify_fold(path,prediction):
    audits=json.loads(path.read_text());seen=[]
    for a in audits:
        train=set(a['train_sample_ids']);test=set(a['test_sample_ids'])
        assert not train&test and a['preprocessing_fit_n']==len(train)
        seen+=list(test)
        for inner in a['inner_splits']:
            tr=set(inner['train_indices']);va=set(inner['validation_indices'])
            assert not tr&va and tr|va==set(range(len(train)))
            assert inner['preprocessing_fit_n']==len(tr)
    assert len(seen)==len(set(seen)) and set(seen)==set(prediction.sample_id)
    return len(audits)

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.repo.resolve()
    checked=[];outer=0
    training=r/'docs/trainset_20261003'
    receipt=json.loads((training/'receipt.json').read_text())
    assert receipt['historical_sources_unchanged'] and receipt['same_training_and_evaluation_participants']
    for model in ['common_ridge','common_elasticnet','full_ridge','full_elasticnet']:
        result=json.loads((training/(model+'_metrics.json')).read_text())
        checked.append(verify_metric(training/(model+'_training_predictions.tsv'),result['metrics'],
                                    np.float32 if 'elasticnet' in model else None))
        assert result['fresh_fit_n']==88 and result['original_refit_max_probability_difference']<=1e-12
    root=r/'docs/chronic_20261003'
    for experiment in ['screen_v1','methylation_transfer_v3']:
        directory=root/experiment;summary=json.loads((directory/'summary.json').read_text());run=json.loads((directory/'receipt.json').read_text())
        script='run_chronic_screen.py' if experiment=='screen_v1' else 'run_chronic_methylation.py'
        assert run['script_sha256']==sha256(r/'scripts'/script)
        assert run['learning_code_sha256']==sha256(r/'src/nutriomics_methylation/chronic.py')
        for item in summary:
            sub=directory/item.get('accession',item.get('source'));ranking=item['ranking']
            m=item['repeated_nested_oof_metrics'];repeated=sub/(ranking+'_repeated_oof.tsv')
            checked.append(verify_metric(repeated,m));matrices=[]
            for seed in SEEDS:
                prediction=sub/f'{ranking}_{seed}_oof.tsv';d=pd.read_csv(prediction,sep='\t')
                expected=next(s for s in item['seed_metrics'] if s['seed']==seed)
                checked.append(verify_metric(prediction,expected));matrices.append(d.probability.to_numpy())
                outer+=verify_fold(sub/f'{ranking}_{seed}_fold_audit.json',d)
            np.testing.assert_allclose(np.mean(matrices,axis=0),pd.read_csv(repeated,sep='\t').probability.to_numpy(),atol=1e-12,rtol=0)
            if 'source_only_frozen_external_metrics' in item:
                external=sub/(ranking+'_external.tsv');checked.append(verify_metric(external,item['source_only_frozen_external_metrics']))
                ext=pd.read_csv(external,sep='\t');assert not set(ext.sample_id)&set(pd.read_csv(repeated,sep='\t').sample_id)
                assert not item['target_fitting_tuning_calibration_or_threshold_selection']
    write_json(a.out,{'generated_utc':datetime.now(timezone.utc).isoformat(),'checked_prediction_tables':len(checked),
        'outer_fold_audits':outer,'metric_denominators_and_recalculation':'passed','inner_outer_membership_checks':'passed',
        'training_refit_matches_frozen_final_predictions':'passed','learned_source_only_transfer_metadata':'passed',
        'prediction_tables':checked,'script_sha256':sha256(Path(__file__)),
        'limits':'Checks exported results/audits. Does not establish donor independence, clinical performance, raw normalization independence or causality.'})
    print(json.dumps({'prediction_tables_verified':len(checked),'outer_fold_audits':outer,'status':'passed'}))

if __name__=='__main__':main()
