import csv
import json
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "myproject"))
from find_data_sources import find_data_sources
from find_candidate_files import find_candidate_files

parser = argparse.ArgumentParser()
parser.add_argument("--input_csv", type=Path, required=True)
parser.add_argument("--out_dir", type=Path, required=True)
parser.add_argument("--workers", type=int, default=8)
args = parser.parse_args()

if args.workers < 1:
    parser.error("--workers must >1")

input_csv = args.input_csv
out_dir = args.out_dir
out_dir.mkdir(parents=True, exist_ok=True)

WORKERS = args.workers

log_lock = threading.Lock()


def record_paper(filename, paper_id, detail=""):
    detail = str(detail).replace("\n", " ").replace("\r", " ").replace("\t", " ")

    with log_lock:
        path = out_dir / filename
        records = {}

        if path.exists():
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    pid = line.split("\t", 1)[0].strip().lower()
                    records[pid] = line

        records[paper_id.strip().lower()] = (
            f"{paper_id}\t{detail}" if detail else paper_id
        )
        path.write_text(
            "\n".join(records.values()) + "\n",
            encoding="utf-8",
        )


def process_paper(paper_id):
    stage = "find_data_sources"

    try:
        sources = find_data_sources(paper_id)

        if not sources:
            record_paper("no_data_sources.txt", paper_id)
            return paper_id, "no_data_sources"

        stage = "find_candidate_files"
        result = find_candidate_files(
            paper_id=paper_id,
            output_dir=out_dir,
            data_sources=sources,
        )

        issues = result.get("issues", [])
        files = result.get("files", [])
        downloaded = any(item.get("local_path") for item in files)
        raw_reads_only = bool(files) and all(
            item.get("candidate_type") == "raw_reads" for item in files
        )

        if issues:
            details = []
            for issue in issues:
                location = issue.get("link") or issue.get("local_path") or ""
                reason = issue.get("reason", "")
                details.append(f"{location}: {reason}")

            record_paper(
                "candidate_download_failed.txt",
                paper_id,
                " | ".join(dict.fromkeys(details)),
            )

        if not downloaded and not issues:
            record_paper(
                "raw_reads_only.txt" if raw_reads_only
                else "no_candidate_files.txt",
                paper_id,
            )

        with log_lock:
            print(json.dumps({
                "paper_id": paper_id,
                "sources": sources,
                "result": result,
            }, ensure_ascii=False, indent=2), flush=True)

        status = (
            "downloaded_with_issues" if downloaded and issues
            else "downloaded" if downloaded
            else "candidate_download_failed" if issues
            else "raw_reads_only" if raw_reads_only
            else "no_candidate_files"
        )
        return paper_id, status

    except Exception as exc:
        record_paper("paper_errors.txt", paper_id, f"{stage}: {exc}")
        return paper_id, f"error: {exc}"


def main():
    with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if "paper_id" not in (reader.fieldnames or []):
            raise ValueError("CSV 缺少 paper_id 列")

        unique_ids = {}
        for row in reader:
            paper_id = (row.get("paper_id") or "").strip()
            if paper_id:
                unique_ids.setdefault(paper_id.lower(), paper_id)

    paper_ids = list(unique_ids.values())
    print(f"Total papers: {len(paper_ids)}; workers: {WORKERS}", flush=True)

    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = [
            executor.submit(process_paper, paper_id)
            for paper_id in paper_ids
        ]
        for index, future in enumerate(as_completed(futures), 1):
            paper_id, status = future.result()
            with log_lock:
                print(
                    f"[{index}/{len(paper_ids)}] {paper_id}: {status}",
                    flush=True,
                )


if __name__ == "__main__":
    main()