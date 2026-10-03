"""Acquire official probe annotations and the published study supplement.

The GPL self-only query deliberately avoids the very large platform-family
archive, which also contains unrelated samples. Raw downloads stay in data/.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import requests

GPL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534&targ=self&form=text&view=full'
SUPP = 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-38439-z/MediaObjects/41467_2023_38439_MOESM3_ESM.xlsx'

def main(root):
    root = Path(root)
    data = root / 'data/chronic_20261003/obesity_annotation_v1'
    out = root / 'docs/chronic_20261003/obesity_annotation_v1'
    data.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    selected = set(row['feature_id'] for row in csv.DictReader((root / 'docs/chronic_20261003/methylation_transfer_v3/obesity_discovery/supervised_refit_features.tsv').open(), delimiter='\t'))
    receipts = []
    for url, name in [(GPL, 'GPL13534_self.soft.gz'), (SUPP, 'published_Supplementary_Data_1-26.xlsx')]:
        dest = data / name
        if not dest.exists():
            response = requests.get(url, stream=True, timeout=(15, 90)); response.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + '.partial')
            n = 0
            with (gzip.open(tmp, 'wb') if name.endswith('.gz') else tmp.open('wb')) as f:
                for block in response.iter_content(1024 * 1024):
                    n += len(block)
                    if n > 600_000_000: raise ValueError('Self-only annotation response unexpectedly large')
                    f.write(block)
            tmp.replace(dest)
        receipts.append({'url': url, 'filename': name, 'bytes': dest.stat().st_size, 'sha256': hashlib.file_digest(dest.open('rb'), 'sha256').hexdigest(), 'license': 'NCBI GEO public platform annotation' if name.endswith('.gz') else 'CC-BY-4.0; McAllan et al., Nat Commun 14:2784, DOI 10.1038/s41467-023-38439-z; subset and reformatted, no original workbook redistributed'})
    rows = []; metadata = []; columns = None; nrows = 0; active = False
    with gzip.open(data / 'GPL13534_self.soft.gz', 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!platform_table_begin'): active = True; continue
            if line.startswith('!platform_table_end'): active = False; continue
            if not active:
                if line.startswith(('^PLATFORM', '!Platform_', '#')): metadata.append(line.rstrip())
                continue
            parts = line.rstrip('\r\n').split('\t')
            if columns is None:
                columns = parts
                if columns[0] != 'ID': raise ValueError('Expected native platform table')
                continue
            nrows += 1
            if parts[0] in selected:
                if len(parts) != len(columns): raise ValueError('Annotation column mismatch')
                rows.append(dict(zip(columns, parts)))
    if {r['ID'] for r in rows} != selected or len(rows) != 50: raise ValueError('All 50 actual selected probes must be annotated exactly once')
    with (out / 'selected_native_GPL13534.tsv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t'); w.writeheader(); w.writerows(rows)
    (out / 'GPL13534_metadata.txt').write_text('\n'.join(metadata) + '\n', encoding='utf-8')
    record = {'sources': receipts, 'native_platform_rows': nrows, 'selected_probes': len(selected), 'matched_probes': len(rows), 'annotation_columns': columns, 'gzip_full_read_crc_passed': True, 'self_only_query': True, 'coordinate_convention': 'Native Illumina MAPINFO; Genome_Build retained verbatim; do not treat as BED coordinates', 'script_sha256': hashlib.file_digest(Path(__file__).open('rb'), 'sha256').hexdigest()}
    (out / 'acquisition.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: record[k] for k in ['native_platform_rows', 'selected_probes', 'matched_probes']}))

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('root'); main(p.parse_args().root)
