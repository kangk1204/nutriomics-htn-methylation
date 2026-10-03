"""Validate full-source coverage rather than accepting an incomplete receipt."""
import importlib.util
from pathlib import Path
import pytest

PATH = Path(__file__).resolve().parents[1] / 'scripts/collect_endocrine_acquisition.py'
SPEC = importlib.util.spec_from_file_location('endocrine_collection', PATH)
COLLECT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COLLECT)


def receipts():
    status = dict(state='complete', bytes=8, official_bytes=8, verified_bytes=8,
                  completed_chunks=2, official_hash_available=False,
                  sha256='1' * 64, crc32='0' * 8, etag='fixed')
    chunks = [dict(number=i, start=i * 4, stop=i * 4 + 3, bytes=4,
                   etag='fixed', sha256='2' * 64, crc32='3' * 8) for i in range(2)]
    return status, chunks


def test_exact_full_coverage_passes():
    COLLECT.validate(*receipts())


def test_incomplete_whole_source_rejected():
    status, chunks = receipts()
    status['state'] = 'running'
    with pytest.raises(ValueError, match='incomplete'):
        COLLECT.validate(status, chunks)


def test_range_gap_rejected():
    status, chunks = receipts()
    chunks[1]['start'] = 5
    with pytest.raises(ValueError, match='range'):
        COLLECT.validate(status, chunks)


def test_mixed_etag_rejected():
    status, chunks = receipts()
    chunks[1]['etag'] = 'changed'
    with pytest.raises(ValueError, match='Mixed upstream'):
        COLLECT.validate(status, chunks)


def test_unverified_extra_bytes_rejected():
    status, chunks = receipts()
    status['verified_bytes'] = 7
    with pytest.raises(ValueError, match='mismatch'):
        COLLECT.validate(status, chunks)
