"""Descriptive measurement availability; no fitting or outcome-based filtering."""
import argparse,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepared',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('New audit output required')
    reports=[]
    for name in ['obesity_discovery','obesity_replication','nafld_GSE48325','nafld_GSE61258']:
        with np.load(a.prepared/(name+'.npz')) as d:
            X=d['X'];samples=d['samples'];y=d['y'];finite=np.isfinite(X)
            measured=finite.any(axis=0);within_missing=(~finite[:,measured]).mean(axis=1)
            allnan=samples[~finite.any(axis=1)]
            reports.append({'cohort':name,'n':len(samples),'controls':int((y==0).sum()),'cases':int(y.sum()),
                'stored_union_probe_ids':X.shape[1],'any_finite_in_cohort_probes':int(measured.sum()),
                'cohort_probes_with_at_least_90percent_finite':int((finite.mean(axis=0)>=.9).sum()),
                'all_measurements_missing_sample_ids':allnan.tolist(),
                'sample_missing_fraction_among_cohort_observed_probes_quantiles':np.quantile(within_missing,[0,.25,.5,.75,1]).tolist(),
                'sample_filtering_performed_by_this_audit':False,'diagnostic_uses_outcome_labels_for_probe_selection':False,
                'interpretation':'Union rows are not all measured in each platform. Availability counts are descriptive; each model learns its own missingness filter inside training folds.'})
    a.out.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'cohorts':reports,'model_results_modified':False},indent=2)+'\n')
    print(json.dumps(reports))

if __name__=='__main__':main()
