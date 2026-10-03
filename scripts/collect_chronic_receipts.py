"""Separate copied raw GEO metadata from versionable derived receipts."""
import argparse,hashlib,json,shutil
from datetime import datetime,timezone
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args();repo=a.repo.resolve()
    acquisition=repo/'docs/chronic_20261003/acquisition_v1';archive=repo/'data/chronic_20261003/collected_native_metadata'
    archive.mkdir(parents=True,exist_ok=True);operations=[]
    for source in sorted(acquisition.glob('*_series_matrix_metadata.tsv')):
        target=archive/source.name
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        if target.exists():
            assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
            source.unlink()
        else:shutil.move(source,target)
        operations.append({'source_filename':source.name,'local_file':str(target),'sha256':digest,'versioned':False,'reason':'Native GEO metadata retained under ignored data; derived sample audits remain versionable.'})
    for source in sorted(acquisition.glob('*_series_matrix_source.json')):
        target=source.with_name(source.name.replace('_series_matrix_source.json','_source_receipt.json'))
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        if target.exists():
            assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
            source.unlink()
        else:source.rename(target)
        operations.append({'source_filename':source.name,'local_file':str(target),'sha256':digest,'versioned':True,'reason':'Copied derived source receipt renamed; original content and all source hashes unchanged.'})
    out=repo/'docs/chronic_20261003/collection_receipt.json'
    if out.exists():raise FileExistsError('A collection receipt already exists')
    out.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'operations':operations,'raw_matrices_and_models_not_collected_or_committed':True},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'native_metadata_retained_ignored':sum(not x['versioned'] for x in operations),'receipts_copied_unchanged':sum(x['versioned'] for x in operations)}))

if __name__=='__main__':main()
