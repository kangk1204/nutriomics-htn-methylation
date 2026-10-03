"""Acquire declared chronic-disease GEO matrices and audit their actual metadata."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

ACCESSIONS = ("GSE61258", "GSE48325", "GSE38642", "GSE222595", "GSE21232", "GSE89632", "GSE135251")

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def acquire(accession, filename, base, raw, docs):
    target=raw/filename
    url=base+filename
    if target.exists(): raise FileExistsError('New acquisition directory required')
    partial=target.with_suffix(target.suffix+'.partial')
    with urllib.request.urlopen(url,timeout=120) as response, partial.open('xb') as f:
        while block:=response.read(1024*1024): f.write(block)
    partial.rename(target)
    headers={};series={};nrows=0;columns=None;in_table=False
    with gzip.open(target,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):in_table=True;continue
            if line.startswith('!series_matrix_table_end'):in_table=False;continue
            values=next(csv.reader([line],delimiter='\t'))
            if in_table:
                if columns is None:columns=values[1:]
                else:
                    if len(values)!=len(columns)+1:raise ValueError('Matrix width mismatch')
                    nrows+=1
            elif line.startswith('!Sample_'):headers.setdefault(values[0],[]).append(values[1:])
            elif line.startswith('!Series_'):series.setdefault(values[0],[]).append(values[1:])
    samples=headers.get('!Sample_geo_accession',[[]])[0]
    if columns is not None and columns!=samples:raise ValueError('Sample metadata/order mismatch')
    meta=pd.DataFrame(index=samples);meta.index.name='sample_id'
    for key,entries in headers.items():
        for j,vals in enumerate(entries):
            if len(vals)!=len(samples):raise ValueError('Sample metadata width mismatch')
            name=key.removeprefix('!Sample_')+(f'_{j+1}' if len(entries)>1 else '')
            meta[name]=vals
            if key.startswith('!Sample_characteristics'):
                for sample,val in zip(samples,vals):
                    if ':' in val:
                        field,value=val.split(':',1)
                        meta.loc[sample,'characteristic_'+re.sub(r'\W+','_',field.strip().lower())]=value.strip()
    metadata_path=docs/(filename.removesuffix('.txt.gz')+'_metadata.tsv')
    meta.to_csv(metadata_path,sep='\t')
    characteristic={c:meta[c].dropna().value_counts().to_dict() for c in meta if c.startswith('characteristic_') or c in ['source_name_ch1','platform_id']}
    result={'accession':accession,'url':url,'file':str(target),'sha256':sha(target),'bytes':target.stat().st_size,
            'gzip_crc':'passed_to_eof','schema_width':'passed','numeric_rows':nrows,'samples':len(samples),
            'sample_ids_unique':len(samples)==len(set(samples)),'matrix_contains_values':nrows>0,
            'series_metadata':series,'metadata_file':str(metadata_path),'metadata_sha256':sha(metadata_path),
            'characteristics':characteristic,'sample_titles':meta.get('title',pd.Series(dtype=str)).tolist(),
            'rights':'Official public GEO research download; cite accession and original study. Source data are excluded from code redistribution.'}
    (docs/(filename.removesuffix('.txt.gz')+'_source.json')).write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'accession':accession,'file':filename,'n':len(samples),'features':nrows,'bytes':result['bytes'],'characteristics':characteristic}),flush=True)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--docs',type=Path,required=True)
    a=p.parse_args();a.raw.mkdir(parents=True,exist_ok=True);a.docs.mkdir(parents=True,exist_ok=True)
    jobs=[];listing=[]
    for accession in ACCESSIONS:
        base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/matrix/'
        text=urllib.request.urlopen(base,timeout=30).read().decode()
        names=sorted(set(re.findall(r'href="([^\"]+\.txt\.gz)"',text)))
        listing.append({'accession':accession,'url':base,'matrix_files':names})
        jobs.extend((accession,name,base) for name in names)
    (a.docs/'acquisition_plan.json').write_text(json.dumps({'fixed_before_outcome_evaluation':True,'accessions':list(ACCESSIONS),'listings':listing,'purpose':'Exploratory NAFLD, T2D and obesity screening; report every evaluable cohort and failure; no selected-positive generalization claim.'},indent=2)+'\n')
    results=[];errors=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(acquire,acc,name,base,a.raw,a.docs):(acc,name) for acc,name,base in jobs}
        for f in as_completed(futures):
            try:results.append(f.result())
            except Exception as e:
                errors.append({'accession':futures[f][0],'file':futures[f][1],'error':type(e).__name__,'detail':str(e)})
                print(json.dumps(errors[-1]),flush=True)
    (a.docs/'acquisition_receipt.json').write_text(json.dumps({'acquired_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),'sources':results,'errors':errors,'all_planned_downloads_finished':len(results)+len(errors)==len(jobs)},ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
