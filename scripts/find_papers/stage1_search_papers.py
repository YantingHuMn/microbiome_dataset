#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _netutil import GLOBAL_THROTTLE  # noqa: E402

csv.field_size_limit(sys.maxsize)  # abstract/matched_queries can exceed the 131072-byte default

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "microbiome-paper-discovery/1.0 (academic research)"}
FIELDS = ["paper_id", "pmid", "pmcid", "doi", "title", "abstract", "journal",
          "publication_year", "publication_type", "is_open_access", "cited_by_count",
          "discovery_source", "matched_queries", "microbe_scope_hint", "assay_group_hint"]


def valid_search_response(data: dict) -> bool:
    if not isinstance(data, dict) or "hitCount" not in data:
        return False
    if int(data["hitCount"]) == 0:
        return True
    return isinstance(data.get("resultList", {}).get("result"), list)


def normalize_identifier(value: str) -> str:
    value = (value or "").strip().lower()
    if value.startswith("https://doi.org/"):
        return value.replace("https://doi.org/", "")
    if value.startswith("http://doi.org/"):
        return value.replace("http://doi.org/", "")
    return value


def merge_value_set(current, incoming):
    if current is None:
        return incoming
    if not incoming:
        return current
    if isinstance(current, set):
        current |= set(str(v).strip() for v in (incoming if isinstance(incoming, (list, tuple, set)) else [incoming]) if str(v).strip())
        return current
    if isinstance(current, str):
        s = {p for p in current.split(";") if p.strip()}
        s |= {str(v).strip() for v in (incoming if isinstance(incoming, (list, tuple, set)) else [incoming]) if str(v).strip()}
        return ";".join(sorted(s))
    return current


def get_json(url: str, attempts: int = 4) -> dict:
    for n in range(attempts):
        try:
            GLOBAL_THROTTLE.wait(url)
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                data = json.loads(r.read())
            if not valid_search_response(data):
                raise ValueError(f"invalid Europe PMC search response: {data}")
            return data
        except Exception as e:
            if n + 1 == attempts:
                raise
            print(f"Europe PMC request failed ({e}); retrying {n + 2}/{attempts}", flush=True)
            time.sleep(2 ** n)
    return {}


def paper_key(r: dict) -> str:
    return normalize_identifier(r.get("doi") or r.get("pmid") or r.get("pmcid") or r.get("id") or "")


def merge_manual_allowlist(papers: dict[str, dict], path: Path) -> None:
    if not path.exists():
        return
    identifier_index = {}
    for existing_key, existing in papers.items():
        for value in (existing_key, existing.get("paper_id"), existing.get("doi"), existing.get("pmid"), existing.get("pmcid")):
            identifier = normalize_identifier(value)
            if identifier:
                identifier_index[identifier] = existing_key
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if not row or not any((v or "").strip() for v in row.values()):
                continue
            identifiers = [normalize_identifier(row.get(field, "")) for field in ("paper_id", "doi", "pmid", "pmcid")]
            identifiers = [identifier for identifier in identifiers if identifier]
            if any(identifier in identifier_index for identifier in identifiers):
                continue
            key = identifiers[0] if identifiers else ""
            if not key:
                continue
            rec = papers.setdefault(key, {
                "paper_id": key,
                "pmid": row.get("pmid", ""),
                "pmcid": row.get("pmcid", ""),
                "doi": row.get("doi", ""),
                "title": row.get("title", ""),
                "abstract": row.get("abstract", ""),
                "journal": row.get("journal", ""),
                "publication_year": row.get("publication_year", ""),
                "publication_type": row.get("publication_type", ""),
                "is_open_access": row.get("is_open_access", ""),
                "cited_by_count": row.get("cited_by_count", ""),
                "discovery_source": row.get("discovery_source", "manual_whitelist"),
                "matched_queries": set(),
                "microbe_scope_hint": set(),
                "assay_group_hint": set(),
            })
            for field, value in row.items():
                if field in {"matched_queries", "microbe_scope_hint", "assay_group_hint"}:
                    if value:
                        for part in str(value).split(";"):
                            if part.strip():
                                rec[field].add(part.strip())
                    continue
                if field not in rec:
                    continue
                if field in {"paper_id", "pmid", "pmcid", "doi"} and rec[field] and value and rec[field] != value:
                    continue
                if not rec[field] and value:
                    rec[field] = value
            if not rec["paper_id"]:
                rec["paper_id"] = key
            rec["discovery_source"] = rec.get("discovery_source") or "manual_whitelist"
            rec["matched_queries"].add("manual_whitelist")
            rec["microbe_scope_hint"].add((row.get("microbe_scope_hint") or "").strip() or "manual")
            rec["assay_group_hint"].add((row.get("assay_group_hint") or "").strip() or "manual")
            for value in (key, rec.get("paper_id"), rec.get("doi"), rec.get("pmid"), rec.get("pmcid")):
                identifier = normalize_identifier(value)
                if identifier:
                    identifier_index[identifier] = key


def run_query(qrow: dict, cache: Path, page_size: int, max_pages: int, sleep: float,
             qi: int, total: int) -> dict[str, dict]:
    """One query's full cursor-paginated fetch. Pagination within a query is
    inherently sequential (page N+1 needs page N's cursor); DIFFERENT
    queries are independent, which is what --workers parallelises across."""
    local: dict[str, dict] = {}
    cursor = "*"; page = 0; query_hits = 0; hit_count = None
    query_hash = hashlib.sha256(qrow["query"].encode("utf-8")).hexdigest()[:8]
    while cursor and page < max_pages:
        page += 1
        cp = cache / f"{qrow['query_id']}_{query_hash}_{page:04d}.json"
        params = {"query": qrow["query"], "format": "json", "pageSize": page_size,
                  "cursorMark": cursor, "resultType": "core"}
        url = API + "?" + urllib.parse.urlencode(params)
        if cp.exists():
            try:
                data = json.loads(cp.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                data = {}
            if not valid_search_response(data):
                print(f"invalid cache ignored: {cp}", flush=True)
                data = get_json(url)
                cp.write_text(json.dumps(data), encoding="utf-8")
                if sleep: time.sleep(sleep)
        else:
            data = get_json(url)
            cp.write_text(json.dumps(data), encoding="utf-8")
            if sleep: time.sleep(sleep)
        if page == 1:
            hit_count = int(data.get("hitCount", 0))
            print(f"[{qi}/{total}] {qrow['query_id']}: hitCount={hit_count}", flush=True)
        result = data.get("resultList", {}).get("result", [])
        if not result:
            break
        query_hits += len(result)
        for r in result:
            key = paper_key(r)
            if not key:
                continue
            row = local.setdefault(key, {
                "paper_id": key, "pmid": r.get("pmid", ""), "pmcid": r.get("pmcid", ""),
                "doi": r.get("doi", ""), "title": r.get("title", ""),
                "abstract": r.get("abstractText", ""), "journal": r.get("journalTitle", ""),
                "publication_year": r.get("pubYear", ""),
                "publication_type": r.get("pubType", ""),
                "is_open_access": r.get("isOpenAccess", ""),
                "cited_by_count": r.get("citedByCount", ""), "discovery_source": "europe_pmc",
                "matched_queries": set(), "microbe_scope_hint": set(), "assay_group_hint": set()})
            row["matched_queries"].add(qrow["query_id"])
            row["microbe_scope_hint"].add(qrow["microbe_scope"])
            row["assay_group_hint"].add(qrow["assay_group"])
        nxt = data.get("nextCursorMark", "")
        if not nxt or nxt == cursor:
            break
        cursor = nxt
    expected_hits = min(hit_count, page_size * max_pages) if hit_count is not None else 0
    if query_hits < expected_hits:
        raise RuntimeError(f"incomplete Europe PMC results for {qrow['query_id']}: "
                           f"retrieved {query_hits} of {hit_count} hits")
    suffix = " TRUNCATED warning: hitCount exceeds page-size x max-pages" if hit_count is not None and hit_count > page_size * max_pages else ""
    print(f"[{qi}/{total}] {qrow['query_id']}: {query_hits} hits{suffix}", flush=True)
    return local


def merge_paper_dict(dest: dict[str, dict], src: dict[str, dict]) -> None:
    """Merge one query's local paper dict into the shared accumulator,
    union-ing the hint sets rather than overwriting when a paper was
    matched by more than one query."""
    for key, row in src.items():
        if key not in dest:
            dest[key] = row
            continue
        d = dest[key]
        for field in ("matched_queries", "microbe_scope_hint", "assay_group_hint"):
            d[field] |= row[field]
        for field in ("pmid", "pmcid", "doi", "title", "abstract", "journal",
                      "publication_year", "publication_type", "is_open_access", "cited_by_count"):
            if not d.get(field) and row.get(field):
                d[field] = row[field]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--queries", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--allowlist", default=None)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--page-size", type=int, default=1000)
    ap.add_argument("--max-pages", type=int, default=100, help="Safety cap per query; truncation is reported")
    ap.add_argument("--sleep", type=float, default=0.2,
                    help="ignored when --workers > 1 -- the shared Throttle paces requests instead")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel queries in flight. Pagination WITHIN one query stays sequential "
                         "(each page needs the previous page's cursor); --workers parallelises across "
                         "DIFFERENT queries instead. 4-8 is reasonable.")
    args = ap.parse_args()

    base = Path(args.out_dir)
    args.queries = args.queries or str(base / "stage0_queries/search_queries.tsv")
    args.out = args.out or str(base / "stage1_papers/papers_master.csv")
    args.allowlist = args.allowlist or str(Path(__file__).resolve().parent / "manual_whitelist.tsv")
    args.cache_dir = args.cache_dir or str(base / "cache/europe_pmc_search")
    
    queries = list(csv.DictReader(open(args.queries, encoding="utf-8"), delimiter="\t"))
    cache = Path(args.cache_dir); cache.mkdir(parents=True, exist_ok=True)
    papers: dict[str, dict] = {}
    sleep = 0.0 if args.workers > 1 else args.sleep
    if args.workers <= 1:
        for qi, qrow in enumerate(queries, 1):
            local = run_query(qrow, cache, args.page_size, args.max_pages, sleep, qi, len(queries))
            merge_paper_dict(papers, local)
    else:
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(run_query, qrow, cache, args.page_size, args.max_pages, sleep, qi, len(queries)): qrow
                   for qi, qrow in enumerate(queries, 1)}
            for fu in as_completed(futs):
                merge_paper_dict(papers, fu.result())  # raises immediately if any query was incomplete
    merge_manual_allowlist(papers, Path(args.allowlist))
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader()
        for row in papers.values():
            for k in ("matched_queries", "microbe_scope_hint", "assay_group_hint"):
                row[k] = ";".join(sorted(row[k]))
            w.writerow(row)
    print(f"wrote {len(papers)} unique papers -> {out}")


if __name__ == "__main__":
    main()
