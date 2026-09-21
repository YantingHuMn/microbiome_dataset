#!/usr/bin/env python3
"""extract accessions, data URLs and availability text from papers."""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

UA = {"User-Agent": "microbiome-paper-link-extractor/1.0 (academic research)"}
PATTERNS = [
    ("bioproject", re.compile(r"\b(PRJ(?:NA|EB|DB|DA|CA)\d{4,10})\b", re.I)),
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
]
URL_RE = re.compile(r"https?://[^\s<>\"']+", re.I)
DATA_HOST = re.compile(r"ncbi\.nlm\.nih\.gov|ebi\.ac\.uk|ena|mgnify|zenodo|figshare|dryad|mendeley|osf\.io|github\.com|dataverse", re.I)


def fetch_xml(pmcid: str, cache: Path) -> str:
    p = cache / f"{pmcid}.xml"
    if p.exists(): return p.read_text(encoding="utf-8", errors="replace")
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
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


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--papers", default=str(base / "stage2_screen/papers_screened.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 3 outputs and PMC XML cache")
    ap.add_argument("--statuses", default="likely_relevant,manual_review")
    ap.add_argument("--sleep", type=float, default=0.15)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    cache = outdir / "cache" / "pmc_xml"
    wanted = set(args.statuses.split(",")); cache.mkdir(parents=True, exist_ok=True)
    links, paper_rows = [], []
    for i, r in enumerate(csv.DictReader(open(args.papers, newline="", encoding="utf-8")), 1):
        if r.get("screen_status") not in wanted: continue
        raw = fetch_xml(r.get("pmcid", ""), cache) if r.get("pmcid") else ""
        full, availability = xml_text_and_availability(raw)
        search_text = " ".join([r.get("title", ""), r.get("abstract", ""), full])
        seen = set()
        for repo, rx in PATTERNS:
            for m in rx.finditer(search_text):
                acc = m.group(1).upper() if repo not in {"zenodo", "figshare", "dryad", "mendeley"} else m.group(1)
                key = (repo, acc)
                if key not in seen:
                    seen.add(key); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                        "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": repo,
                        "accession": acc, "raw_url": "", "link_type": "accession", "source": "fulltext" if raw else "abstract"})
        for u in URL_RE.findall(html.unescape(raw)):
            u = u.rstrip('.,;)\\"]')
            if DATA_HOST.search(u) and ("url", u) not in seen:
                seen.add(("url", u)); links.append({"paper_id": r["paper_id"], "pmid": r.get("pmid", ""),
                    "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""), "repository": "url",
                    "accession": "", "raw_url": u, "link_type": "data_url", "source": "fulltext"})
        rr = dict(r); rr.update({"fulltext_available": bool(raw), "data_availability_text": availability,
                                "n_data_links": len(seen)})
        paper_rows.append(rr)
        if i % 100 == 0: print(f"processed {i}", flush=True)
        if raw: time.sleep(args.sleep)
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