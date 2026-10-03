"""Resume the official unnormalized S-BSST831 beta matrix with bounded ranges.

Run on an authorized Linux server only after checking the previous host's
acquisition process is stopped or stationary. This process refuses another
known acquisition on the same host. At most two HTTP worker threads are used.
CRC32 and SHA256 verify local resume chunks, not an unavailable official hash.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import zlib

import requests

API = "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BSST831"
RAW = "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/831/S-BSST831/Files/BetaValue_Raw.txt"
NAME = "BetaValue_Raw.txt"


def utc(): return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def digest(path):
    hasher = hashlib.sha256(); crc = 0; count = 0
    with path.open("rb") as source:
        for data in iter(lambda: source.read(1024 * 1024), b""):
            hasher.update(data); crc = zlib.crc32(data, crc); count += len(data)
    return {"bytes": count, "sha256": hasher.hexdigest(), "crc32": f"{crc & 0xffffffff:08x}"}


def check_range_headers(response, start, stop, total, etag):
    if response.status_code != 206:
        raise ValueError(f"Server ignored Range or failed: HTTP {response.status_code}")
    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
    if match is None or tuple(map(int, match.groups())) != (start, stop, total):
        raise ValueError("Content-Range differs from requested official interval")
    if response.headers.get("Content-Encoding", "identity").lower() != "identity":
        raise ValueError("Encoded transport is not admitted for byte-range stitching")
    if etag and response.headers.get("ETag") != etag:
        raise ValueError("Official file ETag changed during acquisition")


def official_file(obj):
    if isinstance(obj, dict):
        if obj.get("path", "").split("/")[-1] == NAME and isinstance(obj.get("size"), int):
            yield obj
        for value in obj.values(): yield from official_file(value)
    elif isinstance(obj, list):
        for value in obj: yield from official_file(value)


def refuse_same_host_duplicate():
    proc = Path("/proc")
    if not proc.exists(): raise ValueError("Run acquisition on the authorized Linux host, not a third local compute host")
    for process in proc.iterdir():
        if not process.name.isdigit() or int(process.name) == os.getpid(): continue
        try: command = (process / "cmdline").read_bytes().decode(errors="replace").replace("\0", " ")
        except (OSError, PermissionError): continue
        if "python" in command and "inspect_extension_sources.py" in command and "--download-beta" in command:
            raise RuntimeError(f"Prior acquisition PID {process.name} still exists; inspect progress before starting another")


def acquire_chunk(folder, number, start, stop, total, etag, retries):
    target = folder / f"{number:05d}.chunk"; receipt = folder / f"{number:05d}.json"
    if target.exists() and receipt.exists():
        old = json.loads(receipt.read_text(encoding="utf-8"))
        if old.get("etag") == etag and old.get("start") == start and old.get("stop") == stop:
            current = digest(target)
            if all(current[k] == old[k] for k in current): return old
    error = None
    for attempt in range(retries):
        temporary = target.with_suffix(".partial")
        try:
            with requests.get(RAW, headers={"Range": f"bytes={start}-{stop}", "Accept-Encoding": "identity"},
                              timeout=(20, 90), stream=True) as response:
                check_range_headers(response, start, stop, total, etag)
                count = 0
                with temporary.open("wb") as destination:
                    for data in response.iter_content(1024 * 1024):
                        count += len(data)
                        if count > stop-start+1: raise ValueError("Range body exceeds declared interval")
                        destination.write(data)
                if count != stop-start+1: raise ValueError("Truncated range body")
            info = {"number": number, "start": start, "stop": stop, "etag": etag, **digest(temporary), "completed_utc": utc()}
            temporary.replace(target); write_json(receipt, info); return info
        except (requests.RequestException, ValueError, OSError) as exc:
            error = str(exc)
            print(json.dumps({"chunk": number, "attempt": attempt+1, "error": error}), flush=True)
            if attempt+1 < retries: time.sleep(min(2 ** attempt, 10))
    raise RuntimeError(f"Range {number} exhausted retries: {error}")


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, choices=[1, 2], default=2)
    parser.add_argument("--chunk-mib", type=int, default=8); parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args(); folder = args.output.resolve(); folder.mkdir(parents=True, exist_ok=True)
    if not 1 <= args.chunk_mib <= 64: raise ValueError("Bound ranges to 1-64 MiB")
    refuse_same_host_duplicate()
    import fcntl
    with (folder / ".range-acquisition.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = folder / "range_acquisition_status.json"
        response = requests.get(API, timeout=(20, 60)); response.raise_for_status(); study = response.json()
        candidates = list(official_file(study))
        if len(candidates) != 1: raise ValueError("Official API must declare exactly one native raw-beta file")
        total = candidates[0]["size"]
        if shutil.disk_usage(folder).free < 2 * total + 1024 ** 3: raise ValueError("Need space for verified chunks and atomic full-file assembly")
        with requests.get(RAW, headers={"Range": "bytes=0-1023", "Accept-Encoding": "identity"}, timeout=(20, 90)) as probe:
            etag = probe.headers.get("ETag"); check_range_headers(probe, 0, 1023, total, etag)
            if not probe.content.startswith(b'"S_1"\t"S_2"'): raise ValueError("Native S_n header does not match the inspected raw matrix")
        if not etag: raise ValueError("A stable ETag is required to pin all resumed chunks")
        write_json(folder / "S-BSST831_study.json", study)
        chunks = folder / ".beta_ranges"; chunks.mkdir(exist_ok=True)
        width = args.chunk_mib * 1024 ** 2
        ranges = [(i, start, min(start+width-1, total-1)) for i, start in enumerate(range(0, total, width))]
        status = {"state": "running", "pid": os.getpid(), "started_utc": utc(), "source_url": RAW,
                  "official_bytes": total, "etag": etag, "threads": args.threads, "chunk_mib": args.chunk_mib,
                  "completed_chunks": 0, "verified_bytes": 0, "official_hash_available": False}
        write_json(status_path, status); records = []
        executor = ThreadPoolExecutor(max_workers=args.threads)
        try:
            futures = [executor.submit(acquire_chunk, chunks, i, start, stop, total, etag, args.retries) for i, start, stop in ranges]
            for future in as_completed(futures):
                records.append(future.result()); status.update(completed_chunks=len(records), verified_bytes=sum(r["bytes"] for r in records), updated_utc=utc())
                write_json(status_path, status); print(json.dumps(status), flush=True)
            temporary = folder / (NAME + ".assembled.partial")
            with temporary.open("wb") as destination:
                for number, _, _ in ranges:
                    with (chunks / f"{number:05d}.chunk").open("rb") as source: shutil.copyfileobj(source, destination, 1024 * 1024)
            full = digest(temporary)
            if full["bytes"] != total: raise ValueError("Assembled file differs from API size")
            temporary.replace(folder / NAME)
            status.update(state="complete", completed_utc=utc(), final_file=NAME, **full)
            write_json(folder / "range_chunk_manifest.json", sorted(records, key=lambda r:r["number"]))
            write_json(status_path, status)
        except Exception as exc:
            status.update(state="failed_resumable", error=str(exc), updated_utc=utc()); write_json(status_path, status); raise
        finally:
            # Pending ranges must not continue after a failed integrity/version check.
            executor.shutdown(wait=True, cancel_futures=True)


if __name__ == "__main__": main()
