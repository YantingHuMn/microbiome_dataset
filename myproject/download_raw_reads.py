#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import signal
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable



COMPLETED_CSV = "download_completed.csv"
FAILED_CSV = "download_failed.csv"
PROGRESS_COLUMNS = [
    "doi", "paper_folder", "title", "status", "n_fastq_expected", "n_downloaded",
    "n_already_verified", "n_failed_files", "unresolved_ena_accessions",
    "unsupported_non_ena_accessions", "checksum_method", "error", "updated_utc",
]
FAILED_COLUMNS = PROGRESS_COLUMNS + ["attempts"]

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession,sample_accession,study_accession,library_layout,"
    "fastq_ftp,fastq_md5,fastq_bytes"
)

# Chunks explicitly describing a different/irrelevant data product are ignored.
# This prevents entries such as "PRJNA494608 (isolate genomes, previous study)"
# from being downloaded merely because an accession appears in the notes.
SKIP_CHUNK = re.compile(
    r"previous study|isolate genome|isolate 16s|not abundance|reference genome|"
    r"assembled genome|genbank",
    re.IGNORECASE,
)
SUPPORTED_ACCESSION = re.compile(
    r"\b(?:PRJNA\d+|PRJEB\d+|PRJDB\d+|SRP\d+|ERP\d+|DRP\d+|"
    r"SRR\d+|ERR\d+|DRR\d+|SAMN\d+|SAMEA\d+|SAMD\d+)\b",
    re.IGNORECASE,
)
UNSUPPORTED_ACCESSION = re.compile(
    r"\b(?:CRA\d+|PRJCA\d+|SUB\d+)\b", re.IGNORECASE
)
RANGE_ACCESSION = re.compile(
    r"\b(PRJNA|PRJEB|PRJDB|SRR|ERR|DRR)(\d+)-(?:\1)?(\d+)\b",
    re.IGNORECASE,
)


class RawReadDownloadError(RuntimeError):
    """Base error for this download step."""


class UnsupportedAccessionError(RawReadDownloadError):
    """Raised when the paper has only a repository type not implemented here."""


@dataclass(frozen=True)
class FastqFile:
    paper_folder: str
    doi: str
    query_accession: str
    run_accession: str
    sample_accession: str
    study_accession: str
    library_layout: str
    url: str
    md5: str
    size_bytes: int

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", 1)[-1]


def normalize_doi(doi: str) -> str:
    value = doi.strip()
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value, flags=re.I)
    return value.rstrip("/").lower()


def read_status_table(status_csv: str | Path) -> list[dict[str, str]]:
    with status_csv.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    seen: dict[str, str] = {}
    for row in rows:
        doi = normalize_doi(row.get("doi", ""))
        folder = (row.get("paper_folder") or "").strip()
        if not doi or not folder:
            raise RawReadDownloadError(
                f"{status_csv}: every row must have a DOI and paper_folder"
            )
        if doi in seen:
            raise RawReadDownloadError(
                f"{status_csv}: duplicate DOI {doi!r} ({seen[doi]!r}, {folder!r})"
            )
        seen[doi] = folder
    return rows


def paper_for_doi(doi: str, status_csv: str | Path) -> dict[str, str]:
    wanted = normalize_doi(doi)
    matches = [r for r in read_status_table(status_csv)
               if normalize_doi(r.get("doi", "")) == wanted]
    if not matches:
        raise RawReadDownloadError(
            f"DOI {doi!r} is not present in {status_csv.name}"
        )
    return matches[0]


def _expand_ranges(text: str) -> str:
    def expand(match: re.Match[str]) -> str:
        prefix, start_s, end_s = match.groups()
        start, end = int(start_s), int(end_s)
        if end < start or end - start > 10_000:
            raise RawReadDownloadError(f"Invalid or excessive accession range: {match.group(0)}")
        return "; ".join(f"{prefix.upper()}{n}" for n in range(start, end + 1))
    return RANGE_ACCESSION.sub(expand, text)


def parse_accessions(raw: str) -> tuple[list[str], list[str]]:
    """Return (ENA-compatible, unsupported) accessions from a status cell."""
    supported: list[str] = []
    unsupported: list[str] = []
    for chunk in re.split(r"[;\n]", _expand_ranges(raw or "")):
        if SKIP_CHUNK.search(chunk):
            continue
        supported.extend(m.group(0).upper() for m in SUPPORTED_ACCESSION.finditer(chunk))
        unsupported.extend(m.group(0).upper() for m in UNSUPPORTED_ACCESSION.finditer(chunk))
    return list(dict.fromkeys(supported)), list(dict.fromkeys(unsupported))


def _open_url(request: urllib.request.Request, timeout: int, retries: int = 4):
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            return urllib.request.urlopen(request, timeout=timeout)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RawReadDownloadError(f"Network request failed: {request.full_url}: {last_error}")


def resolve_ena_accession(
    accession: str,
    paper_folder: str,
    doi: str,
    timeout: int = 120,
) -> list[FastqFile]:
    params = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": ENA_FIELDS,
        "format": "tsv",
        "download": "true",
    })
    request = urllib.request.Request(
        f"{ENA_API}?{params}", headers={"User-Agent": "microbiome-raw-read-downloader/1.0"}
    )
    with _open_url(request, timeout=timeout) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    files: list[FastqFile] = []
    for row in rows:
        urls = [x for x in (row.get("fastq_ftp") or "").split(";") if x]
        md5s = (row.get("fastq_md5") or "").split(";")
        sizes = (row.get("fastq_bytes") or "").split(";")
        if not urls:
            continue
        if not (len(urls) == len(md5s) == len(sizes)):
            raise RawReadDownloadError(
                f"ENA returned mismatched URL/MD5/size lists for "
                f"{row.get('run_accession') or accession}"
            )
        for url, md5, size in zip(urls, md5s, sizes):
            if url.startswith("ftp://"):
                url = "https://" + url[6:]
            elif not url.startswith(("http://", "https://")):
                url = "https://" + url
            files.append(FastqFile(
                paper_folder=paper_folder,
                doi=normalize_doi(doi),
                query_accession=accession,
                run_accession=(row.get("run_accession") or accession).strip(),
                sample_accession=(row.get("sample_accession") or "").strip(),
                study_accession=(row.get("study_accession") or "").strip(),
                library_layout=(row.get("library_layout") or "").strip(),
                url=url,
                md5=md5.strip().lower(),
                size_bytes=int(size),
            ))
    return files


def resolve_paper_fastqs(
    doi: str,
    status_csv: str | Path,
    timeout: int = 120,
) -> tuple[dict[str, str], list[FastqFile], list[str], list[str]]:
    row = paper_for_doi(doi, status_csv)
    accessions, unsupported = parse_accessions(row.get("raw_data_accessions", ""))
    if not accessions and unsupported:
        raise UnsupportedAccessionError(
            f"{row['paper_folder']} has unsupported repository accessions only: "
            f"{', '.join(unsupported)}"
        )
    if not accessions:
        raise RawReadDownloadError(
            f"No ENA-compatible raw-read accession found for {row['paper_folder']}"
        )

    resolved: list[FastqFile] = []
    unresolved: list[str] = []
    for accession in accessions:
        found = resolve_ena_accession(
            accession, row["paper_folder"], row["doi"], timeout=timeout
        )
        if found:
            resolved.extend(found)
        else:
            unresolved.append(accession)

    # The same run may be reached through both a BioProject and an explicitly
    # listed SRR/ERR.  Deduplicate by the physical FASTQ URL.
    unique = {f.url: f for f in resolved}
    files = sorted(unique.values(), key=lambda f: (f.run_accession, f.filename))
    if not files:
        details = ", ".join(unresolved + unsupported) or "none"
        raise RawReadDownloadError(
            f"No public FASTQ URLs resolved for {row['paper_folder']}; "
            f"unresolved/unsupported accessions: {details}"
        )
    return row, files, unresolved, unsupported


def _md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _verify(path: Path, item: FastqFile, checksum: bool = True) -> None:
    actual_size = path.stat().st_size
    if actual_size != item.size_bytes:
        raise RawReadDownloadError(
            f"Size mismatch for {path}: expected {item.size_bytes}, got {actual_size}"
        )
    if checksum and item.md5:
        actual_md5 = _md5(path)
        if actual_md5 != item.md5:
            raise RawReadDownloadError(
                f"MD5 mismatch for {path}: expected {item.md5}, got {actual_md5}"
            )


def _download_one(
    item: FastqFile,
    paper_dir: Path,
    timeout: int,
    checksum: bool,
) -> tuple[Path, str]:
    run_dir = paper_dir / "raw_reads" / item.run_accession
    run_dir.mkdir(parents=True, exist_ok=True)
    final = run_dir / item.filename
    partial = final.with_name(final.name + ".part")

    if final.exists():
        _verify(final, item, checksum=checksum)
        return final, "already_verified"

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > item.size_bytes:
        raise RawReadDownloadError(
            f"Partial file is larger than expected; inspect manually: {partial}"
        )
    headers = {"User-Agent": "microbiome-raw-read-downloader/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = urllib.request.Request(item.url, headers=headers)
    with _open_url(request, timeout=timeout) as response:
        status = getattr(response, "status", None)
        if offset and status != 206:
            # Server ignored Range. Restart safely rather than append duplicate data.
            offset = 0
            mode = "wb"
        else:
            mode = "ab" if offset else "wb"
        with partial.open(mode) as handle:
            shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)

    _verify(partial, item, checksum=checksum)
    os.replace(partial, final)
    return final, "downloaded"


class Progress:
    """Two CSVs under output_root, rewritten atomically after every paper.

    Keyed by normalised DOI. A paper is written as `in_progress` before work
    starts, so a hard kill still leaves a marker that resume can act on.
    """

    def __init__(self, output_root: Path, status_rows: list[dict[str, str]] | None = None):
        self.root = Path(output_root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.titles = {normalize_doi(r.get("doi", "")): r for r in (status_rows or [])}
        self.done = self._load(self.root / COMPLETED_CSV)
        self.bad = self._load(self.root / FAILED_CSV)
        # A DOI in both tables is ambiguous; keep the unfinished record and let
        # it be reprocessed rather than silently trusting the completed one.
        for doi in set(self.done) & set(self.bad):
            print(f"  ! {doi} is in both progress tables; will reprocess", flush=True)
            self.done.pop(doi, None)

    @staticmethod
    def _load(path: Path) -> dict[str, dict]:
        if not path.exists():
            return {}
        with path.open(encoding="utf-8-sig", newline="") as fh:
            return {normalize_doi(r.get("doi", "")): r for r in csv.DictReader(fh)
                    if (r.get("doi") or "").strip()}

    @staticmethod
    def _write(path: Path, rows: list[dict], columns: list[str]) -> None:
        temp = path.with_suffix(path.suffix + ".tmp")
        with temp.open("w", encoding="utf-8-sig", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp, path)

    def save(self) -> None:
        key = lambda r: (r.get("paper_folder") or "").lower()
        self._write(self.root / COMPLETED_CSV,
                    sorted(self.done.values(), key=key), PROGRESS_COLUMNS)
        self._write(self.root / FAILED_CSV,
                    sorted(self.bad.values(), key=key), FAILED_COLUMNS)

    def _base(self, doi: str, row: dict[str, str] | None) -> dict:
        doi = normalize_doi(doi)
        src = row or self.titles.get(doi) or {}
        return {
            "doi": doi,
            "paper_folder": (src.get("paper_folder") or "").strip(),
            "title": (src.get("title") or "").strip(),
            "n_fastq_expected": 0, "n_downloaded": 0, "n_already_verified": 0,
            "n_failed_files": 0, "unresolved_ena_accessions": "",
            "unsupported_non_ena_accessions": "", "checksum_method": "",
            "error": "", "updated_utc": _now(),
        }

    def start(self, doi: str, row: dict[str, str] | None = None) -> None:
        """Mark in_progress before any work. Not a failure - just 'not confirmed done'."""
        doi = normalize_doi(doi)
        prev = self.bad.get(doi, {})
        record = {**self._base(doi, row), **{k: v for k, v in prev.items() if v},
                  "status": "in_progress", "error": "", "updated_utc": _now()}
        record["attempts"] = str(int(prev.get("attempts") or 0) + 1)
        self.bad[doi] = record
        self.done.pop(doi, None)
        self.save()

    def finish(self, doi: str, status: str, stats: dict,
               row: dict[str, str] | None = None) -> None:
        doi = normalize_doi(doi)
        prev = self.bad.get(doi, {})
        record = {**self._base(doi, row), **stats,
                  "status": status, "updated_utc": _now()}
        if status == "complete":
            self.bad.pop(doi, None)
            record.pop("attempts", None)
            self.done[doi] = record
        else:
            record["attempts"] = prev.get("attempts") or "1"
            self.done.pop(doi, None)
            self.bad[doi] = record
        self.save()

    def state_of(self, doi: str) -> str:
        doi = normalize_doi(doi)
        if doi in self.done:
            return "complete"
        return (self.bad.get(doi, {}).get("status") or "") if doi in self.bad else ""


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def manifest_intact(paper_dir: Path) -> bool:
    """Re-check a `complete` paper cheaply: manifest present, every file the right size."""
    manifest = paper_dir / "raw_reads" / "download_manifest.json"
    if not manifest.exists():
        return False
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    files = payload.get("files") or []
    if not files:
        return False
    for item in files:
        path = paper_dir / "raw_reads" / item.get("run_accession", "") / \
            str(item.get("url", "")).rsplit("/", 1)[-1]
        if not path.exists() or path.stat().st_size != int(item.get("size_bytes") or 0):
            return False
    return True


def _write_manifest(
    paper_dir: Path,
    status_csv: str | Path,
    row: dict[str, str],
    files: list[FastqFile],
    unresolved: list[str],
    unsupported: list[str],
) -> Path:
    raw_dir = paper_dir / "raw_reads"
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = raw_dir / "download_manifest.json"
    payload = {
        "paper_folder": row["paper_folder"],
        "doi": normalize_doi(row["doi"]),
        "source_status_file": str(status_csv),
        "raw_data_accessions_original": row.get("raw_data_accessions", ""),
        "unresolved_ena_accessions": unresolved,
        "unsupported_non_ena_accessions": unsupported,
        "n_fastq_files": len(files),
        "total_bytes": sum(f.size_bytes for f in files),
        "files": [asdict(f) for f in files],
    }
    temp = manifest.with_suffix(".json.tmp")
    temp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temp, manifest)
    return manifest


def download_raw_reads(
    doi: str,
    output_root: str | Path,
    status_csv: str | Path,
    *,
    workers: int = 3,
    timeout: int = 120,
    checksum: bool = True,
    manifest_only: bool = False,
) -> dict[str, object]:
    """Download public ENA FASTQs for one DOI into its paper folder.

    Parameters
    ----------
    doi:
        DOI present in ``1_2_status_raw_reads_only.csv``.
    output_root:
        Parent directory. Files go to
        ``output_root/<paper_folder>/raw_reads/<run_accession>/``.
    workers:
        Number of concurrent FASTQ downloads for this paper.
    checksum:
        Verify ENA MD5 after download (strongly recommended).
    manifest_only:
        Resolve and record the FASTQ manifest without downloading file content.
    """
    status_path = Path(status_csv)
    row, files, unresolved, unsupported = resolve_paper_fastqs(
        doi, status_csv=status_path, timeout=timeout
    )
    paper_dir = Path(output_root) / row["paper_folder"]
    manifest = _write_manifest(paper_dir, status_csv, row, files, unresolved, unsupported)
    required_bytes = sum(f.size_bytes for f in files)

    result: dict[str, object] = {
        "paper_folder": row["paper_folder"],
        "doi": normalize_doi(row["doi"]),
        "paper_dir": str(paper_dir),
        "manifest": str(manifest),
        "n_fastq_files": len(files),
        "total_bytes": required_bytes,
        "unresolved_ena_accessions": unresolved,
        "unsupported_non_ena_accessions": unsupported,
        "downloaded": 0,
        "already_verified": 0,
        "failed_files": 0,
        "status": "",
        "error": "",
        "checksum_method": "md5" if checksum else "size_only (no MD5)",
    }
    if manifest_only:
        # Resolving a manifest proves nothing about the bytes on disk.
        result["status"] = "manifest_only"
        return result

    usage = shutil.disk_usage(paper_dir)
    existing_bytes = sum(
        min((paper_dir / "raw_reads" / f.run_accession / f.filename).stat().st_size,
            f.size_bytes)
        for f in files
        if (paper_dir / "raw_reads" / f.run_accession / f.filename).exists()
    )
    need = required_bytes - existing_bytes
    if need > usage.free:
        raise RawReadDownloadError(
            f"Insufficient free space for {row['paper_folder']}: need approximately "
            f"{need:,} bytes, free {usage.free:,} bytes"
        )

    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {
            pool.submit(_download_one, item, paper_dir, timeout, checksum): item
            for item in files
        }
        for future in as_completed(futures):
            item = futures[future]
            try:
                path, state = future.result()
                result[state] = int(result[state]) + 1
                print(f"[{row['paper_folder']}] {state}: {path}", flush=True)
            except Exception as exc:  # collect all failures from concurrent jobs
                errors.append(f"{item.run_accession}/{item.filename}: {exc}")
    # Per-file failures are reported, not raised: the caller must be able to tell
    # "some files missing" (partial/failed, retryable) from "nothing resolved".
    result["failed_files"] = len(errors)
    if errors:
        result["status"] = "failed"
        result["error"] = (f"{len(errors)} FASTQ download(s) failed: "
                           + "; ".join(errors))[:1000]
    elif unresolved or unsupported:
        # Every byte we knew about arrived, but the paper lists more than we reached.
        result["status"] = "partial"
        result["error"] = ("all resolved FASTQs OK; unreachable accessions: "
                           + ", ".join(unresolved + unsupported))[:1000]
    else:
        result["status"] = "complete"
    return result


def process_paper(
    row: dict[str, str],
    output_root: str | Path,
    progress: Progress,
    status_csv: str | Path,
    *,
    workers: int = 3,
    timeout: int = 120,
    checksum: bool = True,
    retry_failed: bool = False,
) -> tuple[str, str]:
    """Download one paper under progress tracking. -> (status, one-line summary).

    Skip rules: `complete` papers are skipped only if their manifest and files
    are still intact; `failed`/`partial` are skipped unless retry_failed;
    `in_progress` (a previous interruption) is always resumed.
    """
    doi = normalize_doi(row.get("doi", ""))
    folder = (row.get("paper_folder") or "").strip()
    state = progress.state_of(doi)
    paper_dir = Path(output_root) / folder

    if state == "complete":
        if manifest_intact(paper_dir):
            return "skipped_complete", "already complete (manifest and files verified)"
        print(f"  ! {folder} marked complete but files are missing or wrong size; "
              f"reprocessing", flush=True)
    elif state in ("failed", "partial") and not retry_failed:
        return "skipped_" + state, f"previously {state}; use --retry-failed to retry"

    progress.start(doi, row)          # persisted before any work begins
    try:
        result = download_raw_reads(
            doi, output_root, status_csv, workers=workers,
            timeout=timeout, checksum=checksum, manifest_only=False,
        )
    except Exception as exc:
        stats = {"error": f"{type(exc).__name__}: {exc}"[:1000],
                 "checksum_method": "md5" if checksum else "size_only (no MD5)"}
        progress.finish(doi, "failed", stats, row)
        return "failed", f"{type(exc).__name__}: {exc}"

    stats = {
        "n_fastq_expected": result["n_fastq_files"],
        "n_downloaded": result["downloaded"],
        "n_already_verified": result["already_verified"],
        "n_failed_files": result["failed_files"],
        "unresolved_ena_accessions": "; ".join(result["unresolved_ena_accessions"]),
        "unsupported_non_ena_accessions": "; ".join(
            result["unsupported_non_ena_accessions"]),
        "checksum_method": result["checksum_method"],
        "error": result["error"],
    }
    progress.finish(doi, result["status"], stats, row)
    return result["status"], (
        f"{result['n_fastq_files']} files, downloaded={result['downloaded']} "
        f"verified={result['already_verified']} failed={result['failed_files']}")


def _human_bytes(value: int) -> str:
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if size < 1024 or unit == "PiB":
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{value} B"


def install_signal_handlers(progress: Progress) -> None:
    """Flush progress on SIGINT/SIGTERM. A SIGKILL cannot be caught, but the
    `in_progress` row written before each paper still makes resume correct."""
    def handler(signum, _frame):
        print(f"\nsignal {signum}: saving progress and stopping", flush=True)
        progress.save()
        raise SystemExit(130)
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, handler)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--doi", help="download one DOI from the raw-reads status table")
    target.add_argument("--all", action="store_true", help="process every DOI in the table")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--status-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--no-checksum", action="store_true")
    parser.add_argument("--retry-failed", action="store_true",
                        help="also retry papers recorded as failed or partial")
    args = parser.parse_args()

    rows = read_status_table(args.status_csv)
    wanted = rows if args.all else [paper_for_doi(args.doi, args.status_csv)]

    if args.manifest_only:
        # Never touches the progress tables: resolving is not downloading.
        for row in wanted:
            result = download_raw_reads(
                row["doi"], args.output_root, args.status_csv,
                timeout=args.timeout, manifest_only=True)
            print(f"{result['paper_folder']}: {result['n_fastq_files']} FASTQ files, "
                  f"{_human_bytes(int(result['total_bytes']))}", flush=True)
        return 0

    progress = Progress(args.output_root, rows)
    install_signal_handlers(progress)
    counts: dict[str, int] = {}
    for i, row in enumerate(wanted, 1):
        print(f"[{i}/{len(wanted)}] {row['paper_folder']}", flush=True)
        status, summary = process_paper(
            row, args.output_root, progress, args.status_csv,
            workers=args.workers, timeout=args.timeout,
            checksum=not args.no_checksum, retry_failed=args.retry_failed)
        counts[status] = counts.get(status, 0) + 1
        print(f"    {status}: {summary}", flush=True)

    print("\n" + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"completed table -> {args.output_root / COMPLETED_CSV}")
    print(f"failed table    -> {args.output_root / FAILED_CSV}")
    return 1 if progress.bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
