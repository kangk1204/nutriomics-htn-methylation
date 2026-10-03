"""Preserve the original receipt and explicitly resolve historical source bytes."""
import argparse,hashlib,json
from pathlib import Path

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args();r=a.repo.resolve()
    receipt=r/'docs/trainset_20261003/receipt.json';original=receipt.with_name('receipt_original.txt')
    if original.exists():raise FileExistsError('Original receipt already preserved; refuse second annotation')
    raw=receipt.read_bytes();value=json.loads(raw)
    maps={};preserved=json.loads((r/'docs/validation/artifact_manifest.json').read_text())['historical_source_preserved_mapping']
    for item in value['source_files']:
        relative=item['path'].split('/nutriomics-htn-methylation/',1)[-1]
        if relative in preserved:
            annotation=preserved[relative]
            assert annotation['executed_sha256']==item['sha256']==digest(r/annotation['archived_path'])
            maps[item['path']]=annotation
    assert len(maps)==2
    original.write_bytes(raw)
    value['historical_source_preserved_mapping']=maps
    value['receipt_original']={'path':'docs/trainset_20261003/receipt_original.txt','sha256':digest(original)}
    value['annotation_scope']='Original source/output hashes and numeric results unchanged. Earlier executed CLI/io bytes are preserved in exact archives; original unannotated receipt is separately retained.'
    receipt.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'annotated':str(receipt),'preserved_sources':len(maps),'original_sha256':digest(original)}))

if __name__=='__main__':main()
