"""Meaningful fail-closed checks for final center-held-out receipts."""
import importlib.util
from pathlib import Path
import pytest

PATH = Path(__file__).resolve().parents[1] / 'scripts/audit_endocrine_outputs.py'
SPEC = importlib.util.spec_from_file_location('endocrine_output_audit', PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def folds():
    return [dict(train_indices=[2, 3], validation_indices=[0, 1], preprocessing_fit_n=2),
            dict(train_indices=[0, 1], validation_indices=[2, 3], preprocessing_fit_n=2)]


def test_rank_auc_preserves_ties():
    assert AUDIT.rank_auroc([0, 1, 0, 1], [.1, .9, .5, .5]) == .875


def test_center_disjoint_partition_passes():
    AUDIT.check_partitions(folds(), [1, 1, 2, 2], 'validation_indices')


def test_center_leak_is_rejected_even_without_patient_overlap():
    with pytest.raises(ValueError, match='center leakage'):
        AUDIT.check_partitions(folds(), [1, 2, 1, 2], 'validation_indices')


def test_repeated_heldout_coverage_is_rejected():
    with pytest.raises(ValueError, match='coverage'):
        AUDIT.check_partitions([folds()[0], folds()[0]], [1, 1, 2, 2], 'validation_indices')


def test_global_preprocessing_size_is_rejected():
    records = folds()
    records[0]['preprocessing_fit_n'] = 4
    with pytest.raises(ValueError, match='training-only'):
        AUDIT.check_partitions(records, [1, 1, 2, 2], 'validation_indices')
