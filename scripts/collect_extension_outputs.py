"""Collect only derived outputs from the authorized isolated healthy server."""
import argparse
import base64
import io
from pathlib import Path
import shlex
import subprocess
import tarfile

def main():
    p=argparse.ArgumentParser();p.add_argument('--remote',required=True);p.add_argument('--destination',type=Path,required=True);p.add_argument('--server',choices=['gpu','cpu'],default='gpu');a=p.parse_args()
    prefix='/home/keunsoo4070/projects/nutriomics-phase1-20261003/' if a.server=='gpu' else '/home/keunsoo/projects/nutriomics-2026-final/'
    if not a.remote.startswith(prefix) or '..' in Path(a.remote).parts:raise ValueError('Only this task healthy server workspace is admitted')
    code=f"import sys,tarfile; t=tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz'); t.add({a.remote!r},arcname={Path(a.remote).name!r},filter=lambda v:None if v.issym() or v.islnk() or '/staged-inputs' in v.name or v.name.endswith(('.joblib','.npz')) or v.name.endswith('native_metadata.tsv') else v); t.close()"
    encoded=base64.b64encode(code.encode()).decode();command='python3 -c '+shlex.quote("import base64;exec(base64.b64decode('"+encoded+"'))")
    transport=['tailscale','ssh','keunsoo4070@keunsoo4070-z790-aorus-elite-ax'] if a.server=='gpu' else ['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','keunsoo@100.82.52.41']
    result=subprocess.run([*transport,command],capture_output=True,timeout=300)
    if result.returncode:raise RuntimeError(result.stderr.decode(errors='replace'))
    a.destination.mkdir(parents=True,exist_ok=True)
    tarfile.open(fileobj=io.BytesIO(result.stdout)).extractall(a.destination,filter='data')
    print('Collected derived output bytes',len(result.stdout))
if __name__=='__main__':main()
