"""Acquire official design metadata and optionally the untouched endocrine beta file.

No outcome/feature selection is performed by this acquisition utility. Published
normalized/ComBat matrices are retained as provenance only and are not eligible
for a fresh nested prediction analysis.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import requests

API = 'https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BSST831'
SUPPLEMENT = 'https://static-content.springer.com/esm/art%3A10.1186%2Fs13148-022-01347-y/MediaObjects/13148_2022_1347_MOESM5_ESM.xlsx'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def fetch(url,path):
    if path.exists():return {'url':url,'path':str(path),'sha256':digest(path),'bytes':path.stat().st_size,'existing_cache':True}
    tmp=path.with_name(path.name+'.partial')
    r=requests.get(url,stream=True,timeout=(30,180));r.raise_for_status()
    expected=r.headers.get('Content-Length');count=0
    with tmp.open('wb') as f:
        for block in r.iter_content(1024*1024):
            f.write(block);count+=len(block)
            if count%(128*1024*1024)<1024*1024:print(json.dumps({'download':path.name,'bytes':count}),flush=True)
    # Content-Length describes transport bytes. requests decompresses an HTTP
    # content-encoded directory page; comparing that header with decoded bytes
    # falsely declares a complete page truncated.
    encoding=r.headers.get('Content-Encoding','identity').lower()
    if expected is not None and encoding=='identity' and count!=int(expected):raise ValueError(f'Incomplete HTTP body: expected {expected}, got {count}')
    tmp.replace(path)
    return {'url':url,'path':str(path),'sha256':digest(path),'bytes':count,'transport_content_length':expected,'transport_encoding':encoding,'retrieved_utc':datetime.now(timezone.utc).isoformat()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--download-beta',action='store_true');p.add_argument('--endocrine-only',action='store_true');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);receipts=[]
    receipts.append(fetch(API,a.output/'S-BSST831_study.json'))
    obj=json.loads((a.output/'S-BSST831_study.json').read_text());print(json.dumps(obj)[:16000],flush=True)
    receipts.append(fetch(SUPPLEMENT,a.output/'sample_characteristics.xlsx'))
    for acc in ([] if a.endocrine_only else ['GSE68336','GSE314180','GSE262817']):
        prefix=acc[:-3]+'nnn'
        receipts.append(fetch(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/soft/{acc}_family.soft.gz',a.output/(acc+'_family.soft.gz')))
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/'
        receipts.append(fetch(url,a.output/(acc+'_matrix_listing.html')))
    if a.download_beta:
        # Select a repository-declared unnormalized/raw beta file, never an
        # adjusted matrix. Ambiguous names require inspection rather than guessing.
        def files(node):
            if isinstance(node,dict):
                if 'path' in node and 'size' in node:yield node
                for value in node.values():yield from files(value)
            elif isinstance(node,list):
                for value in node:yield from files(value)
        candidates=[f for f in files(obj) if 'raw' in f['path'].lower() and ('beta' in f['path'].lower() or 'methyl' in f['path'].lower())]
        if len(candidates)!=1:raise ValueError(f'Expected exactly one raw-beta file, got {candidates}')
        name=candidates[0]['path'];link=f'https://www.ebi.ac.uk/biostudies/files/S-BSST831/{name}'
        record=fetch(link,a.output/Path(name).name)
        if record['bytes']!=candidates[0]['size']:raise ValueError('Raw beta size differs from official BioStudies file listing')
        receipts.append(record)
    (a.output/'acquisition_receipt.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
if __name__=='__main__':main()
