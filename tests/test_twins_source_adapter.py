import gzip
import importlib.util
from pathlib import Path
import numpy as np
import pytest

script=Path(__file__).resolve().parents[1]/'scripts/evaluate_twins_transfer.py'
spec=importlib.util.spec_from_file_location('twins',script);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def source(path):
    with gzip.open(path,'wt') as f:
        for family in range(1,36):
            for member in (1,2):
                bmi='high' if member==1 else 'low';state=('BMI_Disc_HighBMI' if member==1 else 'BMI_Disc_LowBMI') if family<=24 else 'BMI_Concordant'
                f.write(f'^SAMPLE = GSM{family*10+member}\n!Sample_characteristics_ch1 = twin id: Twin{family}_{member}\n!Sample_characteristics_ch1 = bmi: {bmi}\n!Sample_characteristics_ch1 = disease state: {state}\n#VALUE = Average Beta\n!sample_table_begin\nID_REF\tVALUE\tDetection Pval\n')
                p=.05 if family==1 and member==1 else .001
                f.write(f'cg00000001\t0.6\t{p}\n!sample_table_end\n')

def test_native_detection_P_mask_and_family_endpoint(tmp_path):
    path=tmp_path/'cohort.soft.gz';source(path);beta,meta,audit=module.extract(path,['cg00000001'])
    assert meta.discordant.sum()==48;assert meta.family_id.nunique()==35
    assert np.isnan(beta.loc['GSM11','cg00000001']);assert audit['selected_beta_failed_detection_p_by_sample']['GSM11']==1
    assert audit['gzip_full_read_crc_passed']

def test_truncated_gzip_is_rejected(tmp_path):
    path=tmp_path/'cohort.soft.gz';source(path);data=path.read_bytes();path.write_bytes(data[:-8])
    with pytest.raises((EOFError,OSError)):module.extract(path,['cg00000001'])
