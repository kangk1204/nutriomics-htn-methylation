"""Size/CRC-verified supplementary downloads with curl transfer retries."""
from __future__ import annotations
import argparse,gzip,hashlib,json,subprocess,urllib.request
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path

FILES=(
 ('GSE48325','matrix','GSE48325_series_matrix.txt.gz'),
 ('GSE61258','matrix','GSE61258_series_matrix.txt.gz'),
 ('GSE222595','suppl','GSE222595_processed_beta_comb.txt.gz'),
)

def fetch(item,dest):
 acc,folder,name=item
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/{folder}/{name}'
 with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=45) as response:
  expected=int(response.headers['Content-Length']);modified=response.headers.get('Last-Modified')
 target=dest/name;partial=dest/(name+'.partial')
 if target.exists():raise FileExistsError('Verified downloads must not be overwritten')
 commands=[]
 for attempt in range(1,9):
  command=['curl','--fail','--location','--silent','--show-error','--connect-timeout','20','--max-time','1200','--retry','3','--retry-all-errors','--continue-at','-','--output',str(partial),url]
  result=subprocess.run(command,capture_output=True,timeout=1500)
  actual=partial.stat().st_size if partial.exists() else 0
  commands.append({'attempt':attempt,'exit_code':result.returncode,'bytes':actual})
  print(json.dumps({'download':name,'attempt':attempt,'bytes':actual,'expected':expected,'exit':result.returncode}),flush=True)
  if actual==expected:break
  if actual>expected:raise ValueError('Download exceeds official content length')
 if not partial.exists() or partial.stat().st_size!=expected:raise ValueError('Download size mismatch after retries')
 with gzip.open(partial,'rb') as f:
  while f.read(1024*1024):pass
 h=hashlib.sha256()
 with partial.open('rb') as f:
  while block:=f.read(1024*1024):h.update(block)
 partial.rename(target)
 return {'accession':acc,'url':url,'file':str(target),'bytes':expected,'last_modified':modified,'sha256':h.hexdigest(),'gzip_crc':'passed_to_eof','attempts':commands}

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise FileExistsError('New supplemental acquisition directory required')
 a.out.mkdir(parents=True)
 results=[];errors=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures={pool.submit(fetch,item,a.out):item for item in FILES}
  for f in as_completed(futures):
   try:results.append(f.result())
   except Exception as e:errors.append({'file':futures[f][2],'type':type(e).__name__,'detail':str(e)});print(json.dumps(errors[-1]),flush=True)
   (a.out/'receipt.json').write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'sources':results,'errors':errors,'scope':'Supplemental beta and corrected complete matrix downloads; original truncated transfer receipts preserved.'},indent=2)+'\n')

if __name__=='__main__':main()
