#!/usr/bin/env python3
"""high-recall title/abstract screening with auditable rules"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

COMMUNITY = re.compile(r"microbiom|microbiota|metagenom|virom|phageom|microbial community|bacterial community|viral community", re.I)
ASSAYS = {
    "16S_amplicon": re.compile(r"16s(?:\s+rrna)?|amplicon sequencing|v[1-9][\-–]?v[1-9]", re.I),
    "shotgun_metagenomics": re.compile(r"shotgun metagenom|whole.metagenom|metagenomic sequencing", re.I),
    "viral_metagenomics": re.compile(r"viral metagenom|virome sequencing|virus.like particle|\bvlp\b", re.I),
    "metatranscriptomics": re.compile(r"metatranscriptom", re.I),
}
PRODUCT = re.compile(r"abundance|taxonomic profil|taxonomic composition|otu|asv|feature table|species profil|genus profil", re.I)
DATA = re.compile(r"bioproject|sequence read archive|\bsra\b|\bena\b|mgnify|zenodo|figshare|dryad|supplement", re.I)
EXCLUDE_TYPE = re.compile(r"\breview\b|systematic review|meta.analysis|editorial|commentary|protocol", re.I)
ISOLATE = re.compile(r"complete genome|whole genome sequenc(?:e|ing).{0,30}(?:isolate|strain)|single isolate|genome announcement|outbreak surveillance", re.I)
HOST_ONLY = re.compile(r"single.cell rna|scrna.seq|whole.exome|host transcriptom|poly\(a\).{0,20}rna.seq", re.I)

# Descriptive metadata only: disease matches never affect screen_status.
DISEASES = [
    ("Crohn's disease", re.compile(r"\bcrohn(?:'s|s)? disease\b|\bCD patients?\b", re.I)),
    ("ulcerative colitis", re.compile(r"\bulcerative colitis\b|\bUC patients?\b", re.I)),
    ("inflammatory bowel disease", re.compile(r"\binflammatory bowel disease\b|\bIBD\b", re.I)),
    ("irritable bowel syndrome", re.compile(r"\birritable bowel syndrome\b|\bIBS\b", re.I)),
    ("colorectal cancer", re.compile(r"\bcolorectal (?:cancer|carcinoma)\b|\bcolon cancer\b|\bCRC\b", re.I)),
    ("gastric cancer", re.compile(r"\bgastric (?:cancer|carcinoma)\b|\bstomach cancer\b", re.I)),
    ("liver cirrhosis", re.compile(r"\b(?:liver|hepatic) cirrhosis\b", re.I)),
    ("non-alcoholic fatty liver disease", re.compile(r"\bnon[- ]alcoholic fatty liver disease\b|\bNAFLD\b|\bNASH\b", re.I)),
    ("alcoholic liver disease", re.compile(r"\balcohol(?:ic|-associated) liver disease\b", re.I)),
    ("type 1 diabetes", re.compile(r"\btype 1 diabetes(?: mellitus)?\b|\bT1D\b", re.I)),
    ("type 2 diabetes", re.compile(r"\btype 2 diabetes(?: mellitus)?\b|\bT2D\b", re.I)),
    ("obesity", re.compile(r"\bobes(?:e|ity)\b", re.I)),
    ("HIV/AIDS", re.compile(r"\bHIV(?:-1)?\b|\bAIDS\b|human immunodeficiency virus", re.I)),
    ("COVID-19", re.compile(r"\bCOVID-?19\b|\bSARS-CoV-2\b", re.I)),
    ("autism spectrum disorder", re.compile(r"\bautism spectrum disorder\b|\bASD\b", re.I)),
    ("Parkinson's disease", re.compile(r"\bParkinson(?:'s|s)? disease\b", re.I)),
    ("Alzheimer's disease", re.compile(r"\bAlzheimer(?:'s|s)? disease\b", re.I)),
    ("multiple sclerosis", re.compile(r"\bmultiple sclerosis\b", re.I)),
    ("rheumatoid arthritis", re.compile(r"\brheumatoid arthritis\b", re.I)),
    ("celiac disease", re.compile(r"\b(?:celiac|coeliac) disease\b", re.I)),
    ("asthma", re.compile(r"\basthma(?:tic)?\b", re.I)),
    ("atopic dermatitis", re.compile(r"\batopic dermatitis\b|\beczema\b", re.I)),
    ("chronic kidney disease", re.compile(r"\bchronic kidney disease\b|\bCKD\b", re.I)),
    ("cardiovascular disease", re.compile(r"\bcardiovascular disease\b", re.I)),
    ("hypertension", re.compile(r"\bhypertension\b", re.I)),
    ("periodontitis", re.compile(r"\bperiodontitis\b|\bperiodontal disease\b", re.I)),
    ("necrotizing enterocolitis", re.compile(r"\bnecrotizing enterocolitis\b|\bNEC\b", re.I)),
    ("cystic fibrosis", re.compile(r"\bcystic fibrosis\b", re.I)),
    ("depression", re.compile(r"\bmajor depressive disorder\b|\bdepression\b", re.I)),
    ("cancer", re.compile(r"\b(?:cancer|carcinoma|tumou?r)\b", re.I)),
]


def extract_diseases(text: str) -> str:
    matches = [name for name, pattern in DISEASES if pattern.search(text)]
    return ";".join(matches) if matches else "unknown"


def main() -> None:
    base = Path(__file__).resolve().parent / "results"
    ap = argparse.ArgumentParser()
    ap.add_argument("--papers", default=str(base / "stage1_papers/papers_master.csv"))
    ap.add_argument("--outdir", required=True, help="Directory for stage 2 outputs")
    args = ap.parse_args()
    rows = []
    for r in csv.DictReader(open(args.papers, newline="", encoding="utf-8")):
        text = f"{r.get('title','')} {r.get('abstract','')}"
        assays = [k for k, rx in ASSAYS.items() if rx.search(text)]
        evidence = []
        if COMMUNITY.search(text): evidence.append("community_terms")
        if assays: evidence.append("assay:" + ";".join(assays))
        if PRODUCT.search(text): evidence.append("abundance_terms")
        if DATA.search(text): evidence.append("data_terms")
        reasons = []
        if EXCLUDE_TYPE.search(text) or EXCLUDE_TYPE.search(r.get("publication_type", "")):
            reasons.append("review_or_non_primary")
        if ISOLATE.search(text): reasons.append("isolate_or_single_genome")
        if HOST_ONLY.search(text) and not COMMUNITY.search(text): reasons.append("host_only_assay")
        if reasons:
            status = "exclude_likely"
        elif COMMUNITY.search(text) and assays:
            status = "likely_relevant"
        elif COMMUNITY.search(text):
            status = "manual_review"
        else:
            status = "low_priority"
        r.update({"screen_status": status, "disease": extract_diseases(text),
                  "assay_prediction": ";".join(assays) or "unknown",
                  "processed_data_hint": bool(PRODUCT.search(text)),
                  "data_link_hint": bool(DATA.search(text)), "screen_evidence": ";".join(evidence),
                  "screen_reason": ";".join(reasons)})
        rows.append(r)
    fields = list(rows[0]) if rows else []
    out = Path(args.outdir) / "papers_screened.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} screened papers -> {out}")


if __name__ == "__main__":
    main()
