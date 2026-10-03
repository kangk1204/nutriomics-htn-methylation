import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/download_endocrine_ranges.py"
spec = importlib.util.spec_from_file_location("endocrine_ranges", SCRIPT)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_range_rejects_full_response_wrong_offset_encoding_or_new_version():
    good = {"Content-Range": "bytes 8-15/100", "ETag": "fixed"}
    module.check_range_headers(SimpleNamespace(status_code=206, headers=good), 8, 15, 100, "fixed")
    variants = [(200, good), (206, {**good,"Content-Range":"bytes 0-7/100"}),
                (206, {**good,"Content-Encoding":"gzip"}), (206, {**good,"ETag":"new"})]
    for status, headers in variants:
        with pytest.raises(ValueError): module.check_range_headers(SimpleNamespace(status_code=status, headers=headers), 8, 15, 100, "fixed")


def test_crc_and_hash_detect_same_length_local_corruption(tmp_path):
    path = tmp_path / "chunk"; path.write_bytes(b"source chunk")
    before = module.digest(path); path.write_bytes(b"source chunX"); after = module.digest(path)
    assert before["bytes"] == after["bytes"]
    assert before["sha256"] != after["sha256"] and before["crc32"] != after["crc32"]
