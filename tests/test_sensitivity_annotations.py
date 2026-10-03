import gzip
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location("sensitivity", Path(__file__).parents[1] / "scripts/run_obesity_sensitivity.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def platform(tmp_path, rows, end=True):
    path = tmp_path / "native.soft.gz"
    with gzip.open(path, "wt") as stream:
        stream.write("^PLATFORM = GPL13534\n!platform_table_begin\nID\tProbe_SNPs\tProbe_SNPs_10\n")
        stream.write(rows)
        if end:
            stream.write("!platform_table_end\n")
    return path


def test_all_platform_flags_not_only_selected_features(tmp_path):
    flagged, seen = module.flagged_platform(platform(tmp_path, "cg1\t\t\ncg2\trs1\t\ncg3\t\trs2\n"))
    assert flagged == {"cg2", "cg3"}
    assert seen == {"cg1", "cg2", "cg3"}


def test_duplicate_and_truncated_annotation_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="Duplicate"):
        module.flagged_platform(platform(tmp_path, "cg1\t\t\ncg1\t\t\n"))
    with pytest.raises(ValueError, match="Incomplete"):
        module.flagged_platform(platform(tmp_path, "cg1\t\t\n", end=False))


def test_column_shift_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="width"):
        module.flagged_platform(platform(tmp_path, "cg1\tmissing-column\n"))
