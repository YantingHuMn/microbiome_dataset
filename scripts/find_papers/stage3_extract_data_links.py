#!/usr/bin/env python3
"""extract accessions, data URLs and availability text from papers.

--workers > 1 fetches multiple papers' PMC full text concurrently (this
stage is network-latency bound, not CPU bound: each paper is an
independent HTTP fetch + local regex pass). A shared per-host Throttle
(_netutil.py) keeps every worker thread paced to the SAME per-host rate
regardless of --workers, so raising --workers shortens wall-clock time by
overlapping wait time across papers, without hitting Europe PMC harder
per second than a single worker would.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _netutil import GLOBAL_THROTTLE  # noqa: E402

csv.field_size_limit(sys.maxsize)  # data_availability_text can exceed the 131072-byte default

UA = {"User-Agent": "microbiome-paper-link-extractor/1.0 (academic research)"}
BIOPROJECT_NUM = re.compile(r"(?:\bbioproject\s*#?\s*|https?://(?:www\.)?ncbi\.nlm\.nih\.gov/bioproject/)(\d{4,10})\b", re.I)
PATTERNS = [
    ("bioproject", re.compile(r"\b(PRJ(?:NA|EB|DB|DA|CA)\d{4,10})\b", re.I)),
    ("bioproject", BIOPROJECT_NUM),
    ("sra_study", re.compile(r"\b([SED]RP\d{5,12})\b", re.I)),
    ("sra_run", re.compile(r"\b([SED]RR\d{5,12})\b", re.I)),
    ("sra_experiment", re.compile(r"\b([SED]RX\d{5,12})\b", re.I)),
    ("biosample", re.compile(r"\b(SAM[NED][A-Z]?\d{4,14})\b", re.I)),
    ("mgnify_study", re.compile(r"\b(MGYS\d{8,})\b", re.I)),
    ("mgnify_analysis", re.compile(r"\b(MGYA\d{8,})\b", re.I)),
    ("arrayexpress", re.compile(r"\b(E-MTAB-\d{2,8})\b", re.I)),
    ("biostudies", re.compile(r"\b(S-(?:BSST|BIAD)\d{2,8})\b", re.I)),
    ("geo", re.compile(r"\b(GSE\d{3,9})\b", re.I)),
    ("gsa", re.compile(r"\b((?:CRA|HRA|PRJCA)\d{4,10})\b", re.I)),
    ("omix", re.compile(r"\b(OMIX\d{3,8})\b", re.I)),
    ("zenodo", re.compile(r"(?:10\.5281/zenodo\.|zenodo\.org/(?:record|records)/)(\d{3,})", re.I)),
    ("figshare", re.compile(r"(?:10\.6084/m9\.figshare\.|figshare\.com/articles/\S*?/)(\d{5,})", re.I)),
    ("dryad", re.compile(r"10\.5061/dryad\.([a-z0-9]+)", re.I)),
    ("mendeley", re.compile(r"data\.mendeley\.com/datasets/([a-z0-9]+)", re.I)),
    ("osf", re.compile(r"osf\.io/([a-z0-9]{5})\b", re.I)),
    ("github", re.compile(r"github\.com/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)", re.I)),
]

# GitHub repos that are analysis tools, not data deposits. In a manual audit
# of an earlier paper batch, 26/61 papers cited one of these; unfiltered, the
# github channel is almost pure noise (pipeline utilities, not abundance
# tables). Extend this set as new tool repos turn up.
GITHUB_TOOL_DENYLIST = {
    "raivokolde/pheatmap", "pmartinezarbizu/pairwiseadonis", "mikemc/speedyseq",
    "joey711/phyloseq", "benjjneb/dada2", "microbiome/microbiome",
    "jfq3/qsrutils", "vegandevs/vegan", "twbattaglia/btools",
    "zdk123/spieceasi", "hallucigenia-sparsa/seqtime", "biobakery/humann",
    "biobakery/metaphlan", "qiime2/qiime2", "ropensci/taxize",
}
URL_RE = re.compile(r"https?://[^\s<>\"']+", re.I)
DATA_HOST = re.compile(r"ncbi\.nlm\.nih\.gov|ebi\.ac\.uk|ena|mgnify|zenodo|figshare|dryad|mendeley|osf\.io|github\.com|dataverse", re.I)


def fetch_xml(pmcid: str, cache: Path) -> str:
    p = cache / f"{pmcid}.xml"
    if p.exists(): return p.read_text(encoding="utf-8", errors="replace")
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
        GLOBAL_THROTTLE.wait(url)
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
            txt = r.read().decode("utf-8", "replace")
        p.write_text(txt, encoding="utf-8"); return txt
    except Exception as e:
        (cache / f"{pmcid}.error.json").write_text(json.dumps({"error": str(e)}))
        return ""


def xml_text_and_availability(raw: str) -> tuple[str, str]:
    if not raw: return "", ""
    try: root = ET.fromstring(raw)
    except ET.ParseError: return re.sub(r"<[^>]+>", " ", raw), ""
    full = " ".join("".join(root.itertext()).split())
    sections = []
    for sec in root.findall(".//sec"):
        title = " ".join("".join(sec.find("title").itertext()).split()) if sec.find("title") is not None else ""
        if re.search(r"data|code|availability|accession", title, re.I):
            sections.append(" ".join("".join(sec.itertext()).split()))
    return full, " | ".join(sections)[:20000]


def process_paper(r: dict, cache: Path) -> tuple[dict, list[dict]]:
    """Everything for one paper: fetch full text, extract links. Independent
    across papers -- this is the unit of work handed to worker threads."""
    links = []
    raw = fetch_xml(r.get("pmcid", ""), cache) if r.get("pmcid") else ""
    full, availability = xml_text_and_availability(raw)
    search_text = " ".join([r.get("title", ""), r.get("abstract", ""), full])
    seen = set()
    for repo, rx in PATTERNS:
        for m in rx.finditer(search_text):
            acc = m.group(1).upper() if repo not in {"zenodo", "figshare", "dryad", "mendeley", "osf", "github"} else m.group(1)
            if repo == "github" and acc.lower() in GITHUB_TOOL_DENYLIST:
                continue
            key = (repo, acc)
            if key not in seen:
                seen.add(key); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                    "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": repo,
                    "accession": acc, "raw_url": "", "link_type": "accession", "source": "fulltext" if raw else "abstract"})
    if r.get("pmcid"):
        key = ("europepmc_supp", r["pmcid"])
        if key not in seen:
            seen.add(key); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": "europepmc_supp",
                "accession": r["pmcid"], "raw_url": "", "link_type": "accession", "source": "pmcid"})
    for m in BIOPROJECT_NUM.finditer(html.unescape(raw)):
        acc = m.group(1)
        key = ("bioproject", acc)
        if key not in seen:
            seen.add(key); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": "bioproject",
                "accession": acc, "raw_url": "", "link_type": "accession", "source": "fulltext"})
    for u in URL_RE.findall(html.unescape(raw)):
        u = u.rstrip('.,;)\\"]')
        if DATA_HOST.search(u) and ("url", u) not in seen:
            seen.add(("url", u)); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": "url",
                "accession": "", "raw_url": u, "link_type": "data_url", "source": "fulltext"})
    rr = dict(r); rr.update({"fulltext_available": bool(raw), "data_availability_text": availability,
                            "n_data_links": len(seen)})
    return rr, links


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--papers", default=str(base / "stage2_screen/papers_screened.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 3 outputs and PMC XML cache")
    ap.add_argument("--statuses", default="likely_relevant,manual_review")
    ap.add_argument("--sleep", type=float, default=0.15,
                    help="ignored when --workers > 1 -- the shared Throttle paces requests instead")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel papers in flight. Each worker still hits any one host no "
                         "faster than a single worker would (see _netutil.Throttle); raising this "
                         "mainly overlaps wait time across DIFFERENT papers/hosts. 8-16 is reasonable.")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    cache = outdir / "cache" / "pmc_xml"
    wanted = set(args.statuses.split(",")); cache.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader(open(args.papers, newline="", encoding="utf-8"))
           if r.get("screen_status") in wanted]

    links, paper_rows = [], []
    if args.workers <= 1:
        for i, r in enumerate(rows, 1):
            rr, ls = process_paper(r, cache)
            paper_rows.append(rr); links.extend(ls)
            if i % 100 == 0: print(f"processed {i}/{len(rows)}", flush=True)
            if rr["fulltext_available"]: time.sleep(args.sleep)
    else:
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(process_paper, r, cache): r for r in rows}
            for i, fu in enumerate(as_completed(futs), 1):
                rr, ls = fu.result()
                paper_rows.append(rr); links.extend(ls)
                if i % 100 == 0: print(f"processed {i}/{len(rows)}", flush=True)
    out = outdir / "paper_data_links.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["paper_id", "pmid", "pmcid", "doi", "repository", "accession", "raw_url", "link_type", "source"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(links)
    pout = outdir / "papers_with_data_text.csv"
    with pout.open("w", newline="", encoding="utf-8") as fh:
        fields2 = list(paper_rows[0]) if paper_rows else []
        w = csv.DictWriter(fh, fieldnames=fields2); w.writeheader(); w.writerows(paper_rows)
    print(f"wrote {len(links)} links -> {out}")


if __name__ == "__main__":
    main()
