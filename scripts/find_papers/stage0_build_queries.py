#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

OBJECT = {
    "bacteria": [
        "microbiome",
        "microbiota",
        "bacteriome",
        '"bacterial community"',
        '"microbial community"'
    ],
    "virus": [
        "virome",
        '"viral metagenome"',
        '"viral community"',
        "phageome",
        "bacteriophage"
    ],
}

ASSAY = {
    "amplicon": [
        '"16S rRNA"',
        '"16S sequencing"',
        '"amplicon sequencing"'
    ],
    "shotgun": [
        '"shotgun metagenomics"',
        '"whole metagenome sequencing"',
        "metagenomic"
    ],
    "virome": [
        '"virome sequencing"',
        '"viral metagenomics"',
        '"virus-like particle sequencing"',
        '"VLP sequencing"'
    ],
    "metatranscriptome": [
        "metatranscriptomic",
        "metatranscriptome"
    ],
}

PRODUCT = [
    "abundance",
    '"relative abundance"',
    '"read count"',
    '"count table"',
    '"taxonomic profile"',
    '"taxonomic composition"',
    '"feature table"',
    '"OTU table"',
    '"ASV table"',
    '"species profile"',
    '"genus profile"',
    "BIOM"
]

RAW_DATA = [
    "FASTQ",
    '"raw reads"',
    '"sequencing reads"',
    '"Sequence Read Archive"',
    "SRA",
    "ENA",
    "BioProject",
    "BioSample"
]


def joined(xs: list[str]) -> str:
    return "(" + " OR ".join(xs) + ")"


def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--out",
        default=str(
            Path(__file__).resolve().parent
            / "results/stage0_queries/search_queries.tsv"
        )
    )

    ap.add_argument(
        "--extra-filter",
        default="",
        help="Optional clause, e.g. AND (human OR patient)"
    )

    args = ap.parse_args()
    rows = []

    for microbe, object_terms in OBJECT.items():

        for assay, assay_terms in ASSAY.items():

            if microbe == "bacteria" and assay == "virome":
                continue

            if microbe == "virus" and assay == "amplicon":
                continue

            rows.append((
                f"{microbe}_{assay}",
                microbe,
                assay,
                "assay",
                f"{joined(object_terms)} AND "
                f"{joined(assay_terms)} "
                f"{args.extra_filter}".strip()
            ))

        rows.append((
            f"{microbe}_processed",
            microbe,
            "unspecified",
            "processed",
            f"{joined(object_terms)} AND "
            f"{joined(PRODUCT)} "
            f"{args.extra_filter}".strip()
        ))

        rows.append((
            f"{microbe}_raw_data",
            microbe,
            "unspecified",
            "raw_data",
            f"{joined(object_terms)} AND "
            f"{joined(RAW_DATA)} "
            f"{args.extra_filter}".strip()
        ))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "query_id",
            "microbe_scope",
            "assay_group",
            "query_type",
            "query"
        ])
        writer.writerows(rows)

    print(f"wrote {len(rows)} queries -> {out}")


if __name__ == "__main__":
    main()