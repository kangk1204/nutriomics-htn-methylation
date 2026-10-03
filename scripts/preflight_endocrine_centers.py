"""Audit actual supplementary participant/center folds without methylation fitting."""
import argparse
import json
from pathlib import Path
import sys
import hashlib
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from nutriomics_methylation.extension_grid import splits


def main():
    p=argparse.ArgumentParser();p.add_argument('--metadata',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    metadata=pd.read_excel(a.metadata).set_index('Sample')
    if metadata.index.duplicated().any():raise ValueError('Repeated supplementary participant IDs')
    if metadata.Hypertension_type.value_counts().to_dict()!={'PHT':42,'CS':57,'PA':101,'PPGL':55}:raise ValueError('Unexpected native counts')
    results=[]
    for subtype in ['CS','PA','PPGL','EHT']:
        select=metadata.Hypertension_type.eq('PHT')|(metadata.Hypertension_type.ne('PHT') if subtype=='EHT' else metadata.Hypertension_type.eq(subtype))
        subset=metadata.loc[select];y=subset.Hypertension_type.ne('PHT').to_numpy(int);groups=subset['Enrolling center'].to_numpy(int)
        coverage=np.zeros(len(y),int);folds=[]
        for number,(train,test) in enumerate(splits(y,20261003,5,groups)):
            if set(groups[train])&set(groups[test]):raise ValueError('Recruitment-center leakage')
            coverage[test]+=1;inner=splits(y[train],20261003+number,3,groups[train])
            calibration=splits(y[train],20261003+number+107,3,groups[train])
            for development,heldout in [*inner,*calibration]:
                if set(groups[train][development])&set(groups[train][heldout]):raise ValueError('Inner/calibration center leakage')
            folds.append({'fold':number,'train_n':len(train),'test_n':len(test),'train_centers':np.unique(groups[train]).tolist(),
                'test_centers':np.unique(groups[test]).tolist(),'test_PHT':int(sum(y[test]==0)),'test_endocrine':int(sum(y[test]==1)),
                'inner_fold_count':len(inner),'calibration_fold_count':len(calibration),
                'test_participant_ids':subset.index[test].tolist(),'center_disjoint_all_levels':True})
        if not np.all(coverage==1):raise ValueError('Incomplete or repeated outer participant coverage')
        results.append({'endpoint':subtype+'_vs_PHT','participants':len(subset),'controls_are_primary_HTN_not_healthy':True,'outer_folds':folds})
    result={'metadata_sha256':hashlib.sha256(a.metadata.read_bytes()).hexdigest(),'participants':len(metadata),
        'all_four_center_partitions_feasible':True,'no_methylation_fit_or_performance_computed':True,'analyses':results}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':main()
