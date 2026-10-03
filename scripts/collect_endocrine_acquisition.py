"""Collect bounded whole-source receipts; never transfer raw beta or weights."""
from __future__ import annotations
import argparse
import base64
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import shlex
import subprocess
import tarfile

REMOTE_ROOT = '/home/keunsoo/projects/nutriomics-2026-final/nutriomics-htn-methylation'


def validate(status, chunks):
    if status['state'] != 'complete' or status['bytes'] != status['official_bytes']:
        raise ValueError('Whole source is incomplete')
    if status['verified_bytes'] != status['bytes'] or status['completed_chunks'] != len(chunks):
        raise ValueError('Whole/range count mismatch')
    if status['official_hash_available']:
        raise ValueError('Unexpected published source checksum claim')
    if not re.fullmatch(r'[0-9a-f]{64}', status['sha256']) or not re.fullmatch(r'[0-9a-f]{8}', status['crc32']):
        raise ValueError('Invalid whole-source local digest')
    offset = 0
    for number, chunk in enumerate(chunks):
        if chunk['number'] != number or chunk['start'] != offset:
            raise ValueError('Missing, overlapping or reordered verified range')
        if chunk['stop'] < chunk['start'] or chunk['bytes'] != chunk['stop'] - chunk['start'] + 1:
            raise ValueError('Range length mismatch')
        if chunk['etag'] != status['etag']:
            raise ValueError('Mixed upstream source versions')
        if not re.fullmatch(r'[0-9a-f]{64}', chunk['sha256']) or not re.fullmatch(r'[0-9a-f]{8}', chunk['crc32']):
            raise ValueError('Range local checksums missing')
        offset = chunk['stop'] + 1
    if offset != status['bytes']:
        raise ValueError('Verified ranges do not exactly cover whole raw source')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError('Fresh immutable collection directory required')
    code = f'''import json,pathlib,sys,tarfile
root=pathlib.Path({REMOTE_ROOT!r});raw=root/'data/endocrine_ranges_20261003_v2'
status=json.loads((raw/'range_acquisition_status.json').read_text())
if status['state']!='complete':raise RuntimeError('Raw whole-file verification incomplete')
archive=tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz')
for name in ('range_acquisition_status.json','range_chunk_manifest.json','S-BSST831_study.json'):
 archive.add(raw/name,arcname=name,recursive=False)
archive.add(root/'data/endocrine_ranges_20261003_v2_job.json',arcname='chain_job_receipt.json',recursive=False)
archive.close()
'''
    encoded = base64.b64encode(code.encode()).decode()
    command = 'python3 -c ' + shlex.quote("import base64;exec(base64.b64decode('" + encoded + "'))")
    run = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8',
                          'keunsoo@100.82.52.41', command], capture_output=True, timeout=120)
    if run.returncode:
        raise RuntimeError(run.stderr.decode(errors='replace'))
    with tarfile.open(fileobj=io.BytesIO(run.stdout)) as archive:
        admitted = {'range_acquisition_status.json', 'range_chunk_manifest.json',
                    'S-BSST831_study.json', 'chain_job_receipt.json'}
        if {m.name for m in archive.getmembers()} != admitted or any(not m.isfile() for m in archive.getmembers()):
            raise ValueError('Unadmitted acquisition collection entry')
        args.output.mkdir(parents=True)
        archive.extractall(args.output, filter='data')
    status = json.loads((args.output / 'range_acquisition_status.json').read_text())
    chunks = json.loads((args.output / 'range_chunk_manifest.json').read_text())
    validate(status, chunks)
    records = [{'name': path.name, 'bytes': path.stat().st_size,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sorted(args.output.iterdir())]
    receipt = {'collected_utc': datetime.now(timezone.utc).isoformat(),
               'whole_raw_verified_complete': True, 'bytes': status['bytes'],
               'whole_raw_sha256': status['sha256'], 'whole_raw_crc32': status['crc32'],
               'source_completed_utc': status['completed_utc'], 'etag': status['etag'],
               'all_verified_ranges_cover_whole_file_once': True, 'ranges': len(chunks),
               'official_published_checksum_available': False,
               'checksum_scope': 'Locally generated downloaded-byte integrity hashes; no official source checksum claimed',
               'raw_or_model_weights_transferred': False, 'model_results_in_this_collection': False,
               'files': records, 'collection_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output / 'collection_receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
