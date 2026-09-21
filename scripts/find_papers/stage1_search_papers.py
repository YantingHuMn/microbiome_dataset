#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "microbiome-paper-discovery/1.0 (academic research)"}
FIELDS = ["paper_id", "pmid", "pmcid", "doi", "title", "abstract", "journal",
          "publication_year", "publication_type", "is_open_access", "cited_by_count",
          "discovery_source", "matched_queries", "microbe_scope_hint", "assay_group_hint"]


def get_json(url: str, attempts: int = 4) -> dict:
    for n in range(attempts):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return json.loads(r.read())
        except Exception:
            if n + 1 == attempts:
                raise
            time.sleep(2 ** n)
    return {}


def paper_key(r: dict) -> str:
    return (r.get("doi") or r.get("pmid") or r.get("pmcid") or r.get("id") or "").lower()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--queries", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--page-size", type=int, default=1000)
    ap.add_argument("--max-pages", type=int, default=100, help="Safety cap per query; truncation is reported")
    ap.add_argument("--sleep", type=float, default=0.2)
    args = ap.parse_args()

    base = Path(args.out_dir)
    args.queries = args.queries or str(base / "stage0_queries/search_queries.tsv")
    args.out = args.out or str(base / "stage1_papers/papers_master.csv")
    args.cache_dir = args.cache_dir or str(base / "cache/europe_pmc_search")
    
    queries = list(csv.DictReader(open(args.queries, encoding="utf-8"), delimiter="\t"))
    cache = Path(args.cache_dir); cache.mkdir(parents=True, exist_ok=True)
    papers: dict[str, dict] = {}
    for qi, qrow in enumerate(queries, 1):
        cursor = "*"; page = 0; query_hits = 0; hit_page_bound = False
        while cursor and page < args.max_pages:
            page += 1
            cp = cache / f"{qrow['query_id']}_{page:04d}.json"
            if cp.exists():
                data = json.loads(cp.read_text(encoding="utf-8"))
            else:
                params = {"query": qrow["query"], "format": "json", "pageSize": args.page_size,
                          "cursorMark": cursor, "resultType": "core"}
                data = get_json(API + "?" + urllib.parse.urlencode(params))
                cp.write_text(json.dumps(data), encoding="utf-8")
                time.sleep(args.sleep)
            result = data.get("resultList", {}).get("result", [])
            if not result:
                break
            query_hits += len(result)
            for r in result:
                key = paper_key(r)
                if not key:
                    continue
                row = papers.setdefault(key, {
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
            if page == args.max_pages:
                hit_page_bound = len(result) >= args.page_size
            cursor = nxt
        suffix = " TRUNCATED" if hit_page_bound else ""
        print(f"[{qi}/{len(queries)}] {qrow['query_id']}: {query_hits} hits{suffix}", flush=True)
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
