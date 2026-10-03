import csv
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pandas as pd

script=Path(__file__).resolve().parents[1]/'scripts/run_endocrine_extension.py'
spec=importlib.util.spec_from_file_location('endocrine',script);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def test_official_R_header_and_paper_ID_adapter(tmp_path):
    labels=np.concatenate([np.repeat(k,v) for k,v in module.EXPECTED.items()]);n=len(labels)
    meta=pd.DataFrame({'Sample':[f'S{i+1}' for i in range(n)],'Hypertension_type':labels,'Sex':'F','Age':50,'Enrolling center':1})
    sheet=tmp_path/'metadata.xlsx';meta.to_excel(sheet,index=False);raw=tmp_path/'BetaValue_Raw.txt'
    with raw.open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow([f'S_{i+1}' for i in range(n)]);writer.writerow(['cg00000001',*np.linspace(.1,.8,n)]);writer.writerow(['cg00000002',*np.linspace(.8,.1,n)])
    X,m,samples,probes=module.prepare(raw,sheet,tmp_path)
    assert X.shape==(255,2);assert samples[0]=='S1';assert np.isclose(X[0,0],.1);assert probes.tolist()==['cg00000001','cg00000002']
