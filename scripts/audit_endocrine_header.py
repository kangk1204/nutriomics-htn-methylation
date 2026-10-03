"""Validate the already CRC-verified first official range without fitting."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from datetime import datetime, timezone

import pandas as pd


def audit(chunk: Path, metadata: Path) -> dict:
    with chunk.open(encoding='utf-8') as stream:
        header = next(csv.reader([stream.readline()], delimiter='\t'))
        rows = [next(csv.reader([stream.readline()], delimiter='\t')) for _ in range(10)]
    implicit = len(rows[0]) == len(header) + 1
    columns = header if implicit else header[1:]
    ids = [re.sub(r'^S_(\d+)$', r'S\1', x) for x in columns]
    meta = pd.read_excel(metadata)
    sample_ids = meta.Sample.astype(str)
    assert len(columns) == 255 and all(len(x) == 256 for x in rows)
    assert len(set(ids)) == 255 and set(ids) == set(sample_ids)
    assert not sample_ids.duplicated().any()
    finite = []
    for row in rows:
        assert row[0].startswith('cg')
        for value in row[1:]:
            if value not in ('NA', 'NaN', 'null', 'NULL'):
                number = float(value)
                assert 0 <= number <= 1
                finite.append(number)
    return {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'first_verified_chunk': str(chunk),
        'chunk_sha256': hashlib.sha256(chunk.read_bytes()).hexdigest(),
        'metadata_sha256': hashlib.sha256(metadata.read_bytes()).hexdigest(),
        'native_participant_columns': len(columns),
        'exact_native_to_paper_one_to_one_mapping': True,
        'mapping_rule': 'S_<integer> -> S<same integer>',
        'implicit_R_row_name_header': implicit,
        'beta_rows_inspected': len(rows),
        'observed_finite_beta_min': min(finite),
        'observed_finite_beta_max': max(finite),
        'whole_file_integrity_checked': False,
        'model_fitted': False,
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--chunk', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = audit(args.chunk, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
