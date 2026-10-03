"""Independently check final fixed predictions, center folds and code receipts.

This audit never fits or selects a model. It requires all declared endpoints
and model families, so a partial result cannot be passed as final completion.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ENDPOINT_N = {'CS_vs_PHT': 99, 'PA_vs_PHT': 143,
              'PPGL_vs_PHT': 97, 'EHT_vs_PHT': 255}
FAMILIES = ('ridge', 'elasticnet', 'linear_svm', 'age_sex_only')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_auroc(y, probability):
    y = np.asarray(y, dtype=int)
    n1, n0 = int(sum(y == 1)), int(sum(y == 0))
    if not n1 or not n0:
        raise ValueError('Both classes required')
    return float((rankdata(probability)[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def check_partitions(records, groups, valid_key):
    groups = np.asarray(groups)
    coverage = np.zeros(len(groups), dtype=int)
    for record in records:
        train = np.asarray(record['train_indices'], dtype=int)
        valid = np.asarray(record[valid_key], dtype=int)
        if len(set(train)) != len(train) or len(set(valid)) != len(valid):
            raise ValueError('Repeated fold indices')
        if set(train) & set(valid) or set(train) | set(valid) != set(range(len(groups))):
            raise ValueError('Fold overlap or incomplete partition')
        if set(groups[train]) & set(groups[valid]):
            raise ValueError('Recruitment-center leakage')
        if record['preprocessing_fit_n'] != len(train):
            raise ValueError('Preprocessing fit size is not training-only')
        coverage[valid] += 1
    if not np.all(coverage == 1):
        raise ValueError('Incomplete or duplicated held-out coverage')


def audit(directory, preflight_path, acquisition_path):
    summary = json.loads((directory / 'summary.json').read_text(encoding='utf-8-sig'))
    protocol = json.loads((directory / 'protocol.json').read_text(encoding='utf-8-sig'))
    preparation = json.loads((directory / 'preparation.json').read_text(encoding='utf-8-sig'))
    preflight = json.loads(preflight_path.read_text(encoding='utf-8-sig'))
    acquisition = json.loads(acquisition_path.read_text(encoding='utf-8-sig'))
    if acquisition['state'] != 'complete' or acquisition['bytes'] != acquisition['official_bytes']:
        raise ValueError('Official raw transfer has not passed whole-file verification')
    if preparation['raw_beta_sha256'] != acquisition['sha256']:
        raise ValueError('Prepared input differs from verified full raw file')
    if protocol['metadata_sha256'] != preflight['metadata_sha256']:
        raise ValueError('Metadata version changed from fold preflight')
    if not summary['all_four_endpoints_reported'] or summary['healthy_control_count'] != 0:
        raise ValueError('Endpoint completeness or control definition changed')
    for name, key in [('run_endocrine_extension.py', 'script_sha256'), ('extension_grid.py', 'extension_grid_sha256')]:
        if sha(directory / 'executed_code' / name) != protocol[key] or summary[key] != protocol[key]:
            raise ValueError('Executed code hash changed')
    expected = {(endpoint, family) for endpoint in ENDPOINT_N for family in FAMILIES}
    actual = {(r['endpoint'], r['family']) for r in summary['analyses']}
    if len(summary['analyses']) != len(expected) or actual != expected:
        raise ValueError('All 16 declared model/endpoint results required')
    rows = []
    for endpoint, family in sorted(expected):
        path = directory / endpoint / (family + '_oof.tsv')
        pred = pd.read_csv(path, sep='\t', float_precision='round_trip')
        y = pred.y.to_numpy(int)
        probability = pred.probability.to_numpy(float)
        groups = pred.center.to_numpy(int)
        ids = pred.sample_id.astype(str).to_numpy()
        if pred.sample_id.duplicated().any() or len(pred) != ENDPOINT_N[endpoint] or sum(y == 0) != 42:
            raise ValueError('Participant count, uniqueness or shared PHT control count changed')
        if set(y) != {0, 1} or not np.isfinite(probability).all() or not ((probability >= 0) & (probability <= 1)).all():
            raise ValueError('Invalid outcome or probability')
        planned = next(r for r in preflight['analyses'] if r['endpoint'] == endpoint)
        if set(pred.fold) != {0, 1, 2}:
            raise ValueError('Three center-held-out folds required')
        for fold in planned['outer_folds']:
            held = pred.fold.eq(fold['fold']).to_numpy()
            if set(ids[held]) != set(fold['test_participant_ids']):
                raise ValueError('Actual outer participant membership changed from preflight')
            if set(groups[held]) & set(groups[~held]):
                raise ValueError('Outer center leakage')
        if family != 'age_sex_only':
            audits = json.loads((directory / endpoint / (family + '_fold_audit.json')).read_text())
            if len(audits) != 3:
                raise ValueError('Missing nested fold audit')
            for record in audits:
                held = pred.fold.eq(record['outer_fold']).to_numpy()
                if record['outer_train_samples'] != list(ids[~held]) or record['outer_test_samples'] != list(ids[held]):
                    raise ValueError('Nested audit patient ordering differs from fixed OOF table')
                check_partitions(record['choice']['inner_splits'], groups[~held], 'validation_indices')
                check_partitions(record['calibration']['calibration_splits'], groups[~held], 'calibration_score_indices')
                if record['calibration']['base_fit_n'] != sum(~held) or record['calibration']['calibrator_fit_n'] != sum(~held):
                    raise ValueError('Outer test data reached model or calibration')
        result = next(r for r in summary['analyses'] if r['endpoint'] == endpoint and r['family'] == family)
        observed_auc = rank_auroc(y, probability)
        observed_brier = float(np.mean((y - probability) ** 2))
        reported = result['nested_center_metrics']
        if abs(observed_auc - reported['auroc']) > 1e-12 or abs(observed_brier - reported['brier']) > 1e-12:
            raise ValueError('Reported AUROC/Brier differs from fixed predictions')
        rows.append({'endpoint': endpoint, 'family': family, 'n': len(pred), 'controls_primary_HTN': int(sum(y == 0)),
                     'rank_auroc': observed_auc, 'brier': observed_brier, 'predictions_sha256': sha(path),
                     'center_leakage_detected': False})
    return {'audited_utc': datetime.now(timezone.utc).isoformat(), 'all_16_results_passed': True,
            'no_model_fit_or_selection': True, 'whole_raw_sha256': acquisition['sha256'],
            'official_source_hash_available': False, 'metadata_sha256': protocol['metadata_sha256'],
            'summary_sha256': sha(directory / 'summary.json'),
            'protocol_sha256': sha(directory / 'protocol.json'),
            'preparation_sha256': sha(directory / 'preparation.json'),
            'metadata_preflight_sha256': sha(preflight_path),
            'acquisition_receipt_sha256': sha(acquisition_path),
            'audit_script_sha256': sha(Path(__file__)),
            'unique_source_participants': 255, 'healthy_control_count': 0,
            'shared_42_controls_across_contrasts_not_independent_people': True, 'analyses': rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--metadata-preflight', type=Path, required=True)
    p.add_argument('--acquisition', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = audit(args.directory, args.metadata_preflight, args.acquisition)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({'all_16_results_passed': True, 'output': str(args.output)}))
