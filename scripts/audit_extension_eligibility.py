"""Report why two additional GEO studies are not binary model test matrices."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re

def record(path,url):
    return {'path':str(path),'url':url,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def audit(source,out):
    out.mkdir(parents=True,exist_ok=True);rows=[]
    for acc in ('GSE314180','GSE262817'):
        metadata=source/(acc+'_series_metadata.txt');brief=source/(acc+'_all_sample_brief.txt');listing=source/(acc+'_supplementary_listing.html')
        series=metadata.read_text(encoding='utf-8');samples=brief.read_text(encoding='utf-8');files=listing.read_text(encoding='utf-8')
        ids=re.findall(r'^\^SAMPLE = (GSM\d+)',samples,re.M);titles=re.findall(r'^!Sample_title = (.+)',samples,re.M)
        if len(ids)!=len(set(ids)) or len(ids)!=len(titles):raise ValueError('Sample-title identity mismatch')
        links=re.findall(r'href="([^"]+)"',files);numerical=[v for v in links if v.endswith(('.csv.gz','.txt.gz','.tsv.gz'))]
        item={'accession':acc,'sample_records':len(ids),'unique_native_titles':len(set(titles)),'official_processed_numeric_files':numerical,
          'retrieval_utc':datetime.now(timezone.utc).isoformat(),'sources':[record(metadata,f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),
          record(brief,f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief'),
          record(listing,f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/')]}
        if acc=='GSE314180':
            hypertensive=sum('hypertensive patients' in t.lower() for t in titles)
            item.update({'native_titles_explicit_hypertensive':hypertensive,'verified_normotensive_controls':0,
              'paper_analysis_patients':176,'participant_equivalence_of_218_GSM_and_paper176_resolved':False,
              'eligible_binary_HTN_vs_normal_test':False,'permitted_role':'Mechanistic leukocyte CpG-inflammatory biomarker support; no disease discrimination test',
              'gate_reason':'All 218 deposited sample titles identify hypertensive patients; no verified normal control outcome. Paper reports176 analyzed patients, exclusion crosswalk unavailable.'})
            if hypertensive!=len(titles):raise ValueError('A sample title requires renewed control-label review')
        else:
            labels=[]
            for title in titles:
                match=re.fullmatch(r'PBMC DNA, (Low Blood Pressure|Elevated/Mid Blood Pressure|High Blood Pressure), ([NMH]\d+)',title)
                if not match:raise ValueError('Native youth BP title requires review: '+title)
                labels.append(match[1])
            counts={label:labels.count(label) for label in sorted(set(labels))}
            item.update({'native_BP_title_counts':counts,'processed_per_person_CpG_matrix_available':bool(numerical),
              'eligible_binary_HTN_vs_normal_test':False,'future_proxy_comparison':'High BP10 versus low BP10; elevated/mid BP10 excluded from fit/tuning',
              'permitted_role':'Group-level DMR/gene mechanism support; raw SRA processing requires a separately completed participant matrix',
              'gate_reason':'Official supplementary files supply a group tDMR/gene workbook and no per-person numerical CpG matrix. DMR group labels cannot be used as sample classifier data.'})
            if counts!={'Elevated/Mid Blood Pressure':10,'High Blood Pressure':10,'Low Blood Pressure':10} or numerical:raise ValueError('Revisit youth source eligibility: design or processed-file availability changed')
        rows.append(item)
    result={'studies':rows,'generated_utc':datetime.now(timezone.utc).isoformat(),'no_classifier_or_AUROC_fabricated':True}
    (out/'eligibility.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();audit(a.source,a.output)
