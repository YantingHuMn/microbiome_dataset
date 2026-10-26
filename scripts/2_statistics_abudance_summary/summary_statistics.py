#!/usr/bin/env python3
"""Database-wide summary statistics for the integrated abundance matrices.

Cohort
------
Start from 1_1_status_integrated.csv (every paper that produced a matrix) and drop
the papers parked in the three triage files:

    1_1_1_manual_review_queue.csv          grouping/labels need a human decision
    1_1_2_status_missing_sample_data.csv   no usable per-sample metadata
    1_1_3_status_missing_feature_data.csv  features are opaque IDs (no taxonomy, no sequence)

What it reports
---------------
Three levels, because they answer different downstream questions:

    paper level    what kind of study / assay / unit / rank is in the DB
    sample level   how many samples, from what body site / host / disease state
    feature level  how many features, and how they split across
                   Bacteria / Archaea / Viruses / Fungi / Eukaryota

Plus a cross-paper taxon catalogue (how many distinct phyla / genera / species,
and in how many papers each one shows up) -- that is the number that tells you
whether cross-study merging at a given rank is realistic.

Outputs: tables as CSV, a readable report as Markdown, figures as SVG
(and PNG as well when matplotlib happens to be importable).

Dependencies: standard library only. matplotlib is optional.

Usage
-----
    python3 summary_statistics.py                  # everything
    python3 summary_statistics.py --no-features    # skip the 1.2 GB feature scan (fast)
    python3 summary_statistics.py --out-dir DIR    # write somewhere else
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------

ROOT = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration")

status_integrated = ROOT / "1_1_status_integrated.csv"
manual_review_queue = ROOT / "1_1_1_manual_review_queue.csv"
missing_feature = ROOT / "1_1_3_status_missing_feature_data.csv"
missing_sample = ROOT / "1_1_2_status_missing_sample_data.csv"

out_dir = Path("/hickory/proj/didonglab/dataset/virus/Database/results/2_statistics_abudance_summary")

# Some feature_metadata.csv rows carry full-length sequences.
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:                                    # 32-bit guard
    csv.field_size_limit(2 ** 31 - 1)


# --------------------------------------------------------------------------
# bucketing rules
#
# Every rule list is ordered: the FIRST pattern that matches wins, so put the
# specific pattern above the general one. Values that match nothing land in the
# trailing "Other" bucket and are dumped verbatim to 99_unmapped_values.csv so
# the rules can be extended against real data instead of guesses.
# --------------------------------------------------------------------------

DOMAIN_RULES = [
    # viruses first: "bacteriophage" contains "bacteria"
    ("Viruses", r"virus|viral|phage|virion|caudoviric|duplodnaviria|monodnaviria|"
                r"varidnaviria|riboviria|uroviricota|anellovir|votu|microviridae"),
    ("Archaea", r"\barchae"),
    ("Bacteria", r"bacteri|cyanobacteri|proteobacteri|firmicut|bacteroidot|actinobacteri"),
    ("Fungi", r"\bfungi|fungal|ascomycot|basidiomycot|mucoromycot|chytridiomycot|yeast"),
    ("Eukaryota", r"eukaryot|eukarya|metazoa|viridiplantae|protist|chromista|"
                  r"alveolata|stramenopil|rhizaria|apicomplexa|opisthokonta"),
]

# paper-level fallback, read off feature_level + reason when per-feature
# taxonomy is unavailable
PAPER_DOMAIN_RULES = [
    # The viral rule must name a viral UNIT, not merely mention viruses: a
    # bacterial paper whose feature_level lists "bacteria, archaea, fungi and
    # viruses" is not a virome, and neither is one whose provenance note says a
    # viral side-table was kept as an alternative.
    ("Viruses", r"\bvotu\b|viral genome|viral cluster|viral family|viral taxon|"
                r"viral operational|viral scaffold|virome|phage genome|"
                r"phage .{0,14}cluster|\bvmag\b|anellovir"),
    ("Fungi", r"\\bits[12]\\b|\\bits[-_ ]?(region|amplicon|gene|rrna|rdna|sequencing|primer)\\b|internal transcribed spacer|fungal|\\bfungi\\b|yeast"),
    ("Eukaryota", r"\b18s\b|protist|eukaryot"),
    ("Bacteria", r"16s|\botu|\basv|rrna|metaphlan|kraken|bracken|motu|\bsgb\b|\bmag\b|"
                 r"genus|species|phylum|family|bacteri|rpob|cpn60|nifh"),
]

ASSAY_RULES = [
    ("Viral metagenomics (virome)", r"\bvotu|virome|viral metagenom|vlp\b|viral genome|vmag"),
    ("Metatranscriptomics", r"metatranscript|rna-seq|metatranscriptom|total rna"),
    ("Metaproteomics", r"metaproteom|spectral|metaprotein"),
    ("qPCR / targeted panel", r"qpcr|upip|taqman|digital pcr|ddpcr"),
    ("Shotgun metagenomics", r"shotgun|metagenom|metaphlan|kraken|bracken|humann|motus|"
                             r"\bmag\b|\bsgb\b|uhgg|gtdb|kaiju"),
    ("ITS amplicon (fungi)", r"\\bits[12]\\b|\\bits[-_ ]?(region|amplicon|gene|rrna|rdna|sequencing|primer)\\b|internal transcribed spacer"),
    ("18S rRNA amplicon", r"\b18s\b"),
    ("Other marker gene", r"cpn60|nifh|rpob|\bcoi\b|amoa"),
    ("16S rRNA amplicon", r"16s|\bv3-v4\b|\bv4\b|\bv1-v2\b|\bv3v4\b|amplicon|\basv\b|\botu\b"),
]

UNIT_RULES = [
    ("Relative abundance (proportion / percent)",
     r"relative abundance|relative frequency|proportion|percent|%|fraction|tss|"
     r"share of reads|compositional"),
    ("Read count (rarefied / equal depth)",
     r"rarefi|subsampl|equal depth|resampl|even depth"),
    ("Coverage / RPKM / normalised",
     r"rpkm|fpkm|\bcpm\b|\btpm|coverage|normalis|normaliz|mean depth|log10|"
     r"spectral abundance|magnitude"),
    ("Read count (raw / filtered)",
     r"count|reads|read-pair|read pair"),
]

RANK_RULES = [
    ("vOTU / viral genome", r"\bvotu|viral genome|viral cluster|viral scaffold|phage genome|vmag"),
    ("MAG / genome bin", r"\bmag\b|genome bin|genome cluster|metagenome-assembled|"
                         r"representative genome|\bsgb\b|\bmsp\b|\bosu\b|uhgg"),
    ("ASV / exact sequence variant", r"\basv\b|exact .*sequence|sequence variant|zotu|\bsv\b"),
    ("OTU", r"\botu\b|operational taxonomic"),
    ("Species", r"species"),
    ("Strain / isolate", r"strain|isolate|reference cluster"),
    ("Genus", r"genus|genera"),
    ("Family", r"family|families"),
    ("Order", r"\border\b"),
    ("Class", r"\bclass\b"),
    ("Phylum", r"phylum|phyla|division"),
    ("Kingdom / Domain", r"kingdom|domain|superkingdom"),
]

# The seven Linnaean ranks, coarse -> fine. Used both to normalise the per-feature
# `taxonomic_rank` cell and to ask how far up a paper's matrix can be aggregated.
LINEAGE_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# A matrix row is not always a Linnaean rank: it can be an ASV, an OTU, a vOTU,
# a MAG... Those are listed first because "species-level vOTU" is a vOTU, and
# "species-level genome bin (SGB)" is a genome bin.
MATRIX_ROW_RANKS = ["Kingdom / Domain", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", "Strain / isolate", "OTU", "ASV / exact sequence variant",
                    "vOTU / viral genome", "MAG / genome bin", "Other / mixed"]

SAMPLE_TYPE_RULES = [
    ("Gut / stool", r"stool|feces|faec|fecal|\bgut\b|caecal|cecal|colon|intestin|rumen|"
                    r"rectal|digesta|midgut|hindgut|ileal|jejun|duoden|manure"),
    ("Oral / saliva", r"saliva|\boral\b|plaque|tongue|dental|gingiv|buccal|supragingival|"
                      r"subgingival|periodont|tooth|teeth|periapical"),
    ("Respiratory", r"sputum|lung|nasal|nasopharyn|oropharyn|airway|bronch|throat|"
                     r"respiratory|pharyn"),
    ("Vaginal / urogenital", r"vagin|cervic|urine|urinary|urethra|semen|penile|bladder|"
                             r"genital|reproductive tract"),
    ("Skin", r"\bskin\b|epiderm|forearm|elbow|axilla|sebum|scalp"),
    ("Milk", r"milk|colostrum"),
    ("Blood / plasma", r"blood|plasma|serum|cfrna"),
    ("Other host tissue / fluid", r"tissue|biopsy|organ|tumou?r|bile|liver|placenta|lesion|"
                                  r"aqueous humor|synovial|csf|mucosa|swab of"),
    ("Invertebrate / host organism", r"\bbee\b|bumble|\bant\b|mosquito|\btick\b|\bfly\b|larva|"
                                     r"insect|snail|isopod|sponge(?! biofilter)|coral|beetle|moth|"
                                     r"nematode|daphnia|arthropod|midgut|shrimp|oyster"),
    ("In vitro culture / bioreactor", r"in vitro|culture|microcosm|mesocosm|bioreactor|chemostat|"
                                      r"shime|synthetic (gut )?community|culture well|fermentation|"
                                      r"activated sludge|reactor"),
    ("Food / fermented product", r"\bfood\b|ferment|cheese|brine|wine|beer|yogurt|meal|mash|"
                                 r"silage|forage"),
    ("Plant / rhizosphere", r"rhizosphere|phyllosphere|\bleaf\b|leaves|\broot\b|\bseed\b|grain|"
                            r"flower|tomato|maize|wheat|\brice\b|plant|pepper|chili|poplar|"
                            r"endophyt|nodule"),
    ("Soil", r"\bsoil\b|permafrost|peat|compost"),
    ("Sediment", r"sediment"),
    ("Water / aquatic", r"water|seawater|marine|ocean|lake|river|effluent|wastewater|sludge|"
                        r"aquarium|groundwater|aquatic|biofilter|estuar"),
    ("Built environment / air", r"\bair\b|\bdust\b|surface|coupon|biofilm|shower|hospital|"
                                r"building|indoor"),
]

HOST_RULES = [
    ("Human", r"homo sapiens|\bhuman\b|\bman\b|patient|participant|infant|child|adult|women|men\b"),
    ("Mouse", r"\bmouse\b|\bmice\b|mus musculus|c57bl|balb/c|cba/j|5xfad|\bko\b mouse"),
    ("Rat", r"\brat\b|rattus|f344|fischer|sprague"),
    ("Pig", r"\bpig\b|swine|sus scrofa|piglet|porcine|boar"),
    ("Cattle / small ruminant", r"cattle|bovine|\bcow\b|bos taurus|sheep|ovis|goat|capra|calf|calves"),
    ("Poultry / bird", r"chicken|gallus|broiler|\bbird\b|poultry|duck|turkey|cockatiel|passerin|"
                       r"avian|penguin"),
    ("Fish", r"\bfish\b|salmon|trout|zebrafish|tilapia|carp\b"),
    ("Dog / cat", r"\bdog\b|canis|beagle|\bcat\b|feline"),
    ("Horse", r"horse|equus|equine|mare|foal"),
    ("Insect / arthropod", r"\bbee\b|bumble|\bant\b|mosquito|anopheles|aedes|\btick\b|drosophila|"
                           r"\bfly\b|insect|beetle|moth|larva|arthropod|sand ?fly|daphnia"),
    ("Other invertebrate", r"sponge|coral|snail|isopod|nematode|oyster|shrimp|mollusc|acropora"),
    ("Plant", r"plant|tomato|maize|wheat|\brice\b|poplar|pepper|chili|hevea|rubber tree|"
              r"arabidopsis|barley|soybean|bean\b"),
    ("Other mammal", r"bat\b|monkey|macaque|primate|deer|rabbit|mammal|seal\b|whale"),
]


# Positive evidence that a sample has no host at all. Only these make a sample
# "Environmental"; a blank host column means unknown, which is a different claim.
ENVIRONMENTAL_RULES = [
    ("Environmental / none", r"\bsoil\b|sediment|\bwater\b|seawater|marine|ocean|lake|river|"
                             r"groundwater|effluent|wastewater|sludge|\bair\b|\bdust\b|"
                             r"permafrost|peat|compost|bioreactor|microcosm|mesocosm|"
                             r"chemostat|activated sludge|biofilter|aquarium|surface|coupon"),
]


def host_of(host_value: str, sample_type: str, paper_hint: str) -> str:
    """Resolve a sample's host, keeping 'unknown' distinct from 'no host'.

    Order: an explicit host column -> a host named inside sample_type ->
    positive environmental evidence -> a paper-level hint -> unknown.
    """
    if (host_value or "").strip():
        return bucket(host_value, HOST_RULES, "Other / unclear")
    named = bucket(sample_type, HOST_RULES, "")
    if named and named != "Not recorded":
        return named
    env = bucket(sample_type, ENVIRONMENTAL_RULES, "")
    if env and env != "Not recorded":
        return env
    return paper_hint or "Not recorded"


# Disease names are free text and the same disease arrives under several
# spellings ("CRC" / "colorectal cancer" / "colon cancer"). Only unambiguous
# synonyms are merged. Abbreviations the source never defined (MS, ME, AS,
# BoneDisease -- all from Yin2025, whose README states they were deliberately
# left unexpanded) stay verbatim and are flagged, because guessing "MS" =
# multiple sclerosis would invent a fact the data does not carry.
DISEASE_SYNONYMS = {
    "crc": "Colorectal cancer",
    "colorectal cancer": "Colorectal cancer",
    "colon cancer (colorectal cancer)": "Colorectal cancer",
    "adenoma": "Colorectal adenoma",
    "colorectal adenoma": "Colorectal adenoma",
    "cd": "Crohn's disease",
    "crohn's disease": "Crohn's disease",
    "uc": "Ulcerative colitis",
    "ulcerative colitis": "Ulcerative colitis",
    "inflammatory bowel disease": "Inflammatory bowel disease (unspecified)",
    "ibd": "Inflammatory bowel disease (unspecified)",
    "t2d": "Type 2 diabetes",
    "type 2 diabetes": "Type 2 diabetes",
    "t1d": "Type 1 diabetes",
    "type 1 diabetes": "Type 1 diabetes",
    "obesity": "Obesity",
    "diet-induced obesity": "Obesity (diet-induced)",
    "impaired glucose tolerance": "Impaired glucose tolerance",
    "parkinson": "Parkinson's disease",
    "parkinson's disease": "Parkinson's disease",
    "autism spectrum disorder": "Autism spectrum disorder",
    "asd": "Autism spectrum disorder",
    "ra": "Rheumatoid arthritis",
    "rheumatoid arthritis": "Rheumatoid arthritis",
    "cystic fibrosis": "Cystic fibrosis",
    "liver cirrhosis": "Liver cirrhosis",
    "psoriasis": "Psoriasis",
}
DISEASE_AMBIGUOUS = {"ms", "me", "as", "bonedisease", "pd", "ad", "cf"}

DISEASE_CATEGORY_RULES = [
    ("Cancer / neoplasm", r"cancer|carcinoma|tumou?r|adenoma|neoplas|myeloma|lymphoma"),
    ("Gastrointestinal", r"crohn|colitis|\bibd\b|inflammatory bowel|\bibs\b|constipation|"
                         r"celiac|coeliac|gastr|enteriti|diverticul"),
    ("Metabolic / hepatic", r"diabet|obesity|glucose|nafld|metabolic|cirrhosis|liver|"
                            r"steato|cholesterol"),
    ("Neurological / psychiatric", r"parkinson|alzheimer|autism|depress|anxiety|cognitive|"
                                   r"dementia|sclerosis|encephalo|stroke|epilep|schizo"),
    ("Autoimmune / inflammatory", r"arthritis|psoriasis|spondyl|lupus|autoimmune|eczema|"
                                  r"dermatitis|allerg"),
    ("Infection", r"infection|malaria|\bhiv\b|covid|salmonella|sars|tubercul|sepsis|"
                  r"challenge|helicobacter|clostridio|\bcdi\b|myiasis"),
    ("Respiratory", r"cystic fibrosis|pulmonary|asthma|copd|apnea|apnoea|bronch|respiratory"),
    ("Cardiovascular", r"cardiovascular|atheroscler|hypertens|cardiac|coronary"),
    ("Oral / dental", r"periodont|caries|gingiv|endodontic|periapical"),
    ("Urogenital", r"\buti\b|urinary|vagin|cervic|prostat|kidney|renal|stone"),
]


def normalise_disease(value: str) -> tuple[str, str]:
    """-> (display name, how it was resolved). Never expands an undefined code."""
    raw = (value or "").strip()
    if not raw:
        return "", ""
    low = raw.lower()
    if low in DISEASE_SYNONYMS:
        return DISEASE_SYNONYMS[low], "merged synonym"
    if low in DISEASE_AMBIGUOUS:
        return raw, "verbatim (abbreviation undefined in source)"
    return raw[:1].upper() + raw[1:], "verbatim"


DOMAIN_NOUNS = r"bacteri|archae|fung|virus|viral|eukaryot|protist"


ALTERNATIVE_CLAUSE = re.compile(
    r"kept (as |separately|out)|retained separately|kept in alternative|"
    r"alternative_versions|as an? alternative|not merged|kept apart", re.I)


def drop_alternative_clauses(text: str) -> str:
    """Remove the clauses that describe tables deliberately kept OUT of the
    main matrix ("ITS fungal OTUs kept separately"), so they cannot decide what
    the main matrix contains."""
    parts = re.split(r"[;|]", text or "")
    keep = [p for p in parts if not ALTERNATIVE_CLAUSE.search(p)]
    return "; ".join(keep) if keep else (text or "")


def paper_domain_of(text: str) -> str:
    """Domain asserted for a whole matrix by its own unit description.

    A description that lists three or more domains ("bacteria, archaea, fungi,
    and viruses") resolves nothing, so it yields Unclassified rather than
    letting whichever rule happens to sit first take every feature.
    """
    low = re.sub(r"\bno[nt]?-?\s*(viral|virus|viruses|fungal|fungi|bacteria|archaea)\b",
                 " ", drop_alternative_clauses(text).lower())
    named = {m.group(0) for m in re.finditer(DOMAIN_NOUNS, low)}
    # bacteria + archaea is the normal prokaryote pairing, not a mixed table
    kinds = {"bacteri": "p", "archae": "p", "fung": "f",
             "virus": "v", "viral": "v", "eukaryot": "e", "protist": "e"}
    if len({kinds.get(n, n) for n in named}) >= 3:
        return "Unclassified"
    return bucket(low, PAPER_DOMAIN_RULES, "Unclassified")


def bucket(text: str, rules, default: str = "Other / unclear") -> str:
    """First matching rule wins; empty text is never forced into a bucket."""
    low = (text or "").lower()
    if not low.strip():
        return "Not recorded"
    for label, pattern in rules:
        if re.search(pattern, low):
            return label
    return default


def unit_bucket(value: str) -> str:
    """Bucket main_matrix_unit, honouring negations.

    Most unit strings read "read count (not rarefied, not normalised)". A plain
    keyword match sees rarefied/normalised and gets the answer exactly backwards.
    The negation can also cover a LIST - "not normalised or rarefied",
    "not rarefied/normalised" - so one `not` has to consume every term it governs,
    and it may be hyphenated ("un-normalised"). A unit that merely says the
    normalisation status is undocumented is not a normalised unit either.
    """
    TERM = r"rarefi\w*|rarif\w*|normali[sz]\w*|subsampl\w*|filtered|thresholded|scaled"
    low = (value or "").lower()
    # "normalization status not stated" -> says nothing about the unit
    low = re.sub(rf"({TERM})\s*(status|state)\b[^;.)]*?\b(not|un)\w*", " ", low)
    # "not/non/un X [or/and// , Y ...]" -> drop the whole negated list
    low = re.sub(rf"\b(not|non|un)[-\s]*({TERM})"
                 rf"((\s*[/,]\s*|\s+(or|and|nor)\s+)({TERM}))*", " ", low)
    return bucket(low, UNIT_RULES, "Other / unclear")


def disease_bucket(value: str) -> str:
    """Collapse the 200 free-text disease_status strings into 5 buckets."""
    low = (value or "").strip().lower()
    if not low:
        return "Not recorded"
    if re.match(r"^(not applicable|n/?a\b|no disease design)", low):
        return "Not applicable (environmental / no disease design)"
    # "not labelled ASD (0CTRL)" negates the DISEASE, and the 0CTRL/control marker
    # says which arm it is. Without this the control arm is read as a missing label.
    if re.match(r"^not label", low) and re.search(r"\bctrl\b|0ctrl|control", low):
        return "Control / healthy"
    if re.match(r"^(not reported|not stated|not specified|unknown|not labelled|"
                r"not labeled|no disease recorded|no label)", low):
        return "Not reported"
    if re.match(r"^(healthy|control|normal|hc\b|ctrl\b|non-?case)", low):
        return "Control / healthy"
    if re.match(r"^(case|disease|patient|infected|challenged|affected)", low):
        return "Case / disease"
    if re.search(r"patient|cancer|carcinoma|tumou?r|diabet|malaria|infect|challeng|diagnos|"
                 r"syndrome|\bitis\b|disorder|\bad\b|\bibd\b|\bibs\b|positive|severe", low):
        return "Case / disease"
    if re.search(r"healthy|control|normal", low):
        return "Control / healthy"
    return "Other label"


# ICTV gives every viral rank a mandatory suffix, so a viral name states its own
# rank regardless of which column it was filed under. Some per-paper integrations
# split the lineage string positionally and treated the leading "Viruses" label as
# a rank, shifting every viral name one rank too deep (Fan2023: the `genus` column
# holds Myoviridae, which is a FAMILY). Trusting the suffix repairs that without
# touching bacterial names, which match none of these.
ICTV_SUFFIX_RANK = (
    ("viricetes", "class"),
    ("viricota", "phylum"),
    ("virales", "order"),
    ("viridae", "family"),
    ("virinae", "family"),      # subfamily; nearest Linnaean slot we keep
    ("virae", "kingdom"),
    ("viria", "kingdom"),       # realm, which sits above kingdom
    # NOT "-virus": ICTV genus names end in -virus, but so do many species names
    # ("Moumouvirus", "Fowlpox virus"), so the suffix cannot tell them apart.
    # The six suffixes above are reserved for exactly one rank each.
)


def ictv_rank(name: str) -> str | None:
    """The rank an ICTV name declares through its suffix, or None if not viral."""
    low = (name or "").strip().lower()
    for suffix, rank in ICTV_SUFFIX_RANK:
        if low.endswith(suffix):
            return rank
    return None


def clean_taxon(value: str) -> str:
    """Strip rank prefixes and confidence suffixes: 'k__Bacteria(99)' -> 'Bacteria'."""
    v = (value or "").strip()
    v = re.sub(r"^[a-z]?_{1,2}", "", v)                   # k__ / d__ / p__ / __
    v = re.sub(r"^D_\d+__", "", v)                        # D_0__
    v = re.sub(r"\(\d+(\.\d+)?\)\s*$", "", v)             # (99)
    return v.strip()


def classify_domain(value: str) -> str | None:
    """The domain a per-feature cell states, or None when it states none.

    Two traps this has to survive, both seen in real files:
      * negation - Lutz2025 labels 372 features "non-viral microbial feature",
        which a plain keyword match reads as viral;
      * several domains at once - Mirauta2026 labels every feature
        "bacteria/archaea (UHGG)", a table-level note that does not resolve the
        individual feature. Returning one of them would invent a fact, so the
        caller falls back to the paper level instead.
    """
    v = clean_taxon(value).lower()
    if not v or v in {"unassigned", "unclassified", "unknown", "none", "na", "norank"}:
        return None
    v = re.sub(r"\bno[nt]?-?\s*(viral|virus|viruses|fungal|fungi|bacterial|bacteria|"
               r"archaeal|archaea|eukaryotic|eukaryota)\b", " ", v)
    hits = []
    for label, pattern in DOMAIN_RULES:
        if re.search(pattern, v) and label not in hits:
            hits.append(label)
    return hits[0] if len(hits) == 1 else None


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def folders_in(path: Path) -> set[str]:
    if not path.exists():
        print(f"  ! {path.name} not found - nothing excluded from it", file=sys.stderr)
        return set()
    return {(r.get("paper_folder") or "").strip() for r in read_csv(path)} - {""}


def as_int(value, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def pct(part: int, whole: int) -> str:
    return f"{100.0 * part / whole:.1f}%" if whole else "-"


def median(values: list[float]) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    mid = len(s) // 2
    return float(s[mid]) if len(s) % 2 else (s[mid - 1] + s[mid]) / 2.0


def write_table(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def counter_table(path: Path, counter: Counter, key_name: str, value_name: str) -> list[tuple]:
    total = sum(counter.values())
    items = counter.most_common()
    write_table(path, [key_name, value_name, "share"],
                [[k, v, pct(v, total)] for k, v in items])
    return items


# --------------------------------------------------------------------------
# figures - self-contained SVG, no dependencies
#
# Palette: validated categorical slots (blue, orange, aqua, yellow, magenta,
# violet) on a light surface. Pies are capped at 6 slices and EVERY slice is
# directly labelled, so identity never rests on colour alone; the long-tail
# distributions use single-hue bars, which sidesteps the question entirely.
# --------------------------------------------------------------------------

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e4e3df"
MAX_PIE_SLICES = 6
FONT = ("-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif")


def _esc(text) -> str:
    return html.escape(str(text), quote=True)


def _fold_small(items: list[tuple], max_slices: int = MAX_PIE_SLICES) -> list[tuple]:
    """Keep the largest slices; everything past the cap becomes one 'Other'."""
    if len(items) <= max_slices:
        return list(items)
    head = list(items[: max_slices - 1])
    tail_total = sum(v for _, v in items[max_slices - 1:])
    head.append((f"Other ({len(items) - max_slices + 1} categories)", tail_total))
    return head


def svg_pie(path: Path, title: str, items: list[tuple], subtitle: str = "") -> None:
    items = [(k, v) for k, v in _fold_small(items) if v > 0]
    if not items:
        return
    total = sum(v for _, v in items)
    width, height = 940, 470
    cx, cy, r = 250, 250, 155

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="{FONT}">',
        f'<rect width="{width}" height="{height}" fill="{SURFACE}"/>',
        f'<text x="32" y="42" font-size="20" font-weight="600" fill="{INK}">{_esc(title)}</text>',
    ]
    if subtitle:
        parts.append(f'<text x="32" y="66" font-size="13" fill="{INK_SOFT}">{_esc(subtitle)}</text>')

    angle = -math.pi / 2                                   # start at 12 o'clock
    for i, (label, value) in enumerate(items):
        sweep = 2 * math.pi * value / total
        end = angle + sweep
        x0, y0 = cx + r * math.cos(angle), cy + r * math.sin(angle)
        x1, y1 = cx + r * math.cos(end), cy + r * math.sin(end)
        large = 1 if sweep > math.pi else 0
        if len(items) == 1:                                # a full circle has no arc
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{SERIES[0]}" '
                         f'stroke="{SURFACE}" stroke-width="2"/>')
        else:
            parts.append(
                f'<path d="M {cx:.2f} {cy:.2f} L {x0:.2f} {y0:.2f} '
                f'A {r} {r} 0 {large} 1 {x1:.2f} {y1:.2f} Z" '
                f'fill="{SERIES[i % len(SERIES)]}" stroke="{SURFACE}" stroke-width="2"/>')
        angle = end

    # direct labels, one row per slice, to the right of the pie
    y = 128
    for i, (label, value) in enumerate(items):
        parts.append(f'<rect x="470" y="{y - 11}" width="13" height="13" rx="3" '
                     f'fill="{SERIES[i % len(SERIES)]}"/>')
        parts.append(f'<text x="494" y="{y}" font-size="14" fill="{INK}">{_esc(label)}</text>')
        parts.append(f'<text x="{width - 32}" y="{y}" font-size="14" fill="{INK_SOFT}" '
                     f'text-anchor="end">{value:,}  ({100.0 * value / total:.1f}%)</text>')
        y += 30

    parts.append(f'<text x="32" y="{height - 20}" font-size="12" fill="{INK_SOFT}">'
                 f'total = {total:,}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_barh(path: Path, title: str, items: list[tuple], subtitle: str = "",
             max_bars: int = 18) -> None:
    items = [(k, v) for k, v in items[:max_bars] if v > 0]
    if not items:
        return
    top = max(v for _, v in items)
    label_w, bar_x, bar_w = 310, 330, 470
    row_h, pad_top = 28, 86
    width = bar_x + bar_w + 120
    height = pad_top + row_h * len(items) + 44

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="{FONT}">',
        f'<rect width="{width}" height="{height}" fill="{SURFACE}"/>',
        f'<text x="32" y="42" font-size="20" font-weight="600" fill="{INK}">{_esc(title)}</text>',
    ]
    if subtitle:
        parts.append(f'<text x="32" y="66" font-size="13" fill="{INK_SOFT}">{_esc(subtitle)}</text>')

    for frac in (0.25, 0.5, 0.75, 1.0):                    # recessive gridlines
        x = bar_x + bar_w * frac
        parts.append(f'<line x1="{x:.1f}" y1="{pad_top - 10}" x2="{x:.1f}" '
                     f'y2="{pad_top + row_h * len(items)}" stroke="{GRID}" stroke-width="1"/>')

    for i, (label, value) in enumerate(items):
        y = pad_top + i * row_h
        w = max(2.0, bar_w * value / top)
        text = label if len(str(label)) <= 44 else str(label)[:42] + "..."
        parts.append(f'<text x="{label_w}" y="{y + 15}" font-size="13" fill="{INK}" '
                     f'text-anchor="end">{_esc(text)}</text>')
        parts.append(f'<rect x="{bar_x}" y="{y + 4}" width="{w:.1f}" height="16" rx="4" '
                     f'fill="{SERIES[0]}"/>')
        parts.append(f'<text x="{bar_x + w + 10:.1f}" y="{y + 17}" font-size="13" '
                     f'fill="{INK_SOFT}">{value:,}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_scatter(path: Path, title: str, points: list[tuple], xlabel: str, ylabel: str,
                subtitle: str = "") -> None:
    """points = [(x, y, label)] on log-log axes."""
    pts = [(x, y, lab) for x, y, lab in points if x > 0 and y > 0]
    if not pts:
        return
    width, height = 820, 580
    left, right, top, bottom = 90, 40, 92, 70
    pw, ph = width - left - right, height - top - bottom

    lx = [math.log10(x) for x, _, _ in pts]
    ly = [math.log10(y) for _, y, _ in pts]
    x0, x1 = math.floor(min(lx)), math.ceil(max(lx))
    y0, y1 = math.floor(min(ly)), math.ceil(max(ly))
    x1 = max(x1, x0 + 1)
    y1 = max(y1, y0 + 1)

    def px(v: float) -> float:
        return left + pw * (math.log10(v) - x0) / (x1 - x0)

    def py(v: float) -> float:
        return top + ph - ph * (math.log10(v) - y0) / (y1 - y0)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="{FONT}">',
        f'<rect width="{width}" height="{height}" fill="{SURFACE}"/>',
        f'<text x="32" y="42" font-size="20" font-weight="600" fill="{INK}">{_esc(title)}</text>',
    ]
    if subtitle:
        parts.append(f'<text x="32" y="66" font-size="13" fill="{INK_SOFT}">{_esc(subtitle)}</text>')

    for e in range(x0, x1 + 1):
        x = px(10 ** e)
        parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + ph}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
        parts.append(f'<text x="{x:.1f}" y="{top + ph + 22}" font-size="12" fill="{INK_SOFT}" '
                     f'text-anchor="middle">10^{e}</text>')
    for e in range(y0, y1 + 1):
        y = py(10 ** e)
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + pw}" y2="{y:.1f}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
        parts.append(f'<text x="{left - 12}" y="{y + 4:.1f}" font-size="12" fill="{INK_SOFT}" '
                     f'text-anchor="end">10^{e}</text>')

    for x, y, _ in pts:
        parts.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="5" fill="{SERIES[0]}" '
                     f'fill-opacity="0.62" stroke="{SURFACE}" stroke-width="2"/>')

    parts.append(f'<text x="{left + pw / 2:.0f}" y="{height - 18}" font-size="13" '
                 f'fill="{INK_SOFT}" text-anchor="middle">{_esc(xlabel)}</text>')
    parts.append(f'<text x="22" y="{top + ph / 2:.0f}" font-size="13" fill="{INK_SOFT}" '
                 f'text-anchor="middle" transform="rotate(-90 22 {top + ph / 2:.0f})">'
                 f'{_esc(ylabel)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def maybe_png(fig_dir: Path, pies: list[tuple], bars: list[tuple]) -> None:
    """Mirror the figures as PNG when matplotlib is importable. Optional by design."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed - SVG only; `pip install matplotlib` to add PNGs)")
        return

    plt.rcParams.update({"font.size": 10, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})
    for name, title, items in pies:
        items = [(k, v) for k, v in _fold_small(items) if v > 0]
        if not items:
            continue
        fig, ax = plt.subplots(figsize=(9.4, 4.7))
        total = sum(v for _, v in items)
        ax.pie([v for _, v in items], startangle=90, counterclock=False,
               colors=SERIES[: len(items)],
               wedgeprops={"edgecolor": SURFACE, "linewidth": 2},
               labels=[f"{k}\n{v:,} ({100.0 * v / total:.1f}%)" for k, v in items],
               textprops={"fontsize": 9, "color": INK})
        ax.set_title(title, fontsize=13, fontweight="600", color=INK, loc="left")
        ax.axis("equal")
        fig.tight_layout()
        fig.savefig(fig_dir / f"{name}.png", dpi=200)
        plt.close(fig)

    for name, title, items in bars:
        items = [(k, v) for k, v in items[:18] if v > 0]
        if not items:
            continue
        fig, ax = plt.subplots(figsize=(9.0, 0.38 * len(items) + 1.9))
        ys = range(len(items))[::-1]
        ax.barh(list(ys), [v for _, v in items], color=SERIES[0], height=0.62)
        ax.set_yticks(list(ys))
        ax.set_yticklabels([k if len(str(k)) <= 44 else str(k)[:42] + "..." for k, _ in items])
        for y, (_, v) in zip(ys, items):
            ax.text(v, y, f" {v:,}", va="center", fontsize=9, color=INK_SOFT)
        ax.set_title(title, fontsize=13, fontweight="600", color=INK, loc="left")
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.xaxis.grid(True, color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        ax.margins(x=0.14)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{name}.png", dpi=200)
        plt.close(fig)
    print(f"  PNG copies written alongside the SVGs")


# --------------------------------------------------------------------------
# per-paper scanning
# --------------------------------------------------------------------------

def scan_samples(folder: Path, paper_host_hint: str = "") -> dict:
    """Read one sample_metadata.csv and derive per-paper sample facts."""
    out = {
        "n_meta_rows": 0, "sample_types": Counter(), "hosts": Counter(),
        "disease": Counter(), "methods": Counter(),
        "n_disease_labelled": 0, "n_group_labelled": 0, "n_subject": 0,
        "n_timepoint": 0, "n_sex": 0, "n_age": 0,
        "n_subjects_distinct": 0, "longitudinal": False,
        "raw_sample_types": Counter(), "raw_hosts": Counter(),
        "diseases": Counter(), "disease_variants": defaultdict(set),
    }
    path = folder / "sample_metadata.csv"
    if not path.exists():
        return out
    try:
        rows = read_csv(path)
    except Exception as exc:                               # keep one bad file from killing the run
        print(f"  ! {folder.name}/sample_metadata.csv unreadable: {exc}", file=sys.stderr)
        return out
    if not rows:
        return out

    cols = {(c or "").strip().lower(): c for c in rows[0].keys() if c}

    def get(row: dict, *names: str) -> str:
        for n in names:
            if n in cols:
                v = row.get(cols[n])
                if isinstance(v, list):
                    v = " ".join(map(str, v))
                v = (v or "").strip()
                if v:
                    return v
        return ""

    subject_counts: Counter = Counter()
    out["n_meta_rows"] = len(rows)
    for row in rows:
        stype = get(row, "sample_type", "sample_type_detail", "body_site")
        host = get(row, "host", "host_species", "host_common_name", "organism")
        dstat = get(row, "disease_status", "disease_state")
        method = get(row, "measurement_method", "assay", "sequencing_method")

        out["sample_types"][bucket(stype, SAMPLE_TYPE_RULES, "Other / unspecified")] += 1
        out["hosts"][host_of(host, stype, paper_host_hint)] += 1
        out["disease"][disease_bucket(dstat)] += 1
        out["methods"][bucket(f"{method} {stype}", ASSAY_RULES)] += 1
        dname, _how = normalise_disease(get(row, "disease", "disease_name"))
        if dname:
            out["diseases"][dname] += 1
            out["disease_variants"][dname].add(get(row, "disease", "disease_name")[:60])
        if stype:
            out["raw_sample_types"][stype[:80]] += 1
        if host:
            out["raw_hosts"][host[:80]] += 1

        if dstat and not re.match(r"^(not reported|not applicable|unknown|not stated)", dstat.lower()):
            out["n_disease_labelled"] += 1
        if get(row, "experimental_group", "group", "treatment"):
            out["n_group_labelled"] += 1
        subj = get(row, "subject_id", "subject", "individual_id", "animal_id")
        if subj:
            out["n_subject"] += 1
            subject_counts[subj] += 1
        if get(row, "timepoint", "time_point", "timepoint_description", "day"):
            out["n_timepoint"] += 1
        if get(row, "sex", "gender"):
            out["n_sex"] += 1
        if get(row, "age_years", "age", "age_group", "age_bin"):
            out["n_age"] += 1

    out["n_subjects_distinct"] = len(subject_counts)
    out["longitudinal"] = any(c > 1 for c in subject_counts.values())
    return out


def scan_features(folder: Path, paper_domain: str) -> dict:
    """Stream one feature_metadata.csv; classify every feature into a domain.

    Streaming matters: the largest file here is ~327 MB / 785k rows.
    """
    out = {
        "n_rows": 0, "domains": Counter(), "domain_source": "none",
        "phyla": set(), "genera": set(), "species": set(),
        "n_with_taxonomy": 0,
        "ictv_corrected": Counter(),     # column rank -> true rank, for viral names
        "row_ranks": Counter(),          # what the matrix rows ARE
        "lineage_filled": Counter(),     # how far up the lineage can be aggregated
        "lineage_cols": [],
    }
    path = folder / "feature_metadata.csv"
    if not path.exists():
        return out

    try:
        fh = path.open(encoding="utf-8-sig", newline="")
    except Exception as exc:
        print(f"  ! {folder.name}/feature_metadata.csv unreadable: {exc}", file=sys.stderr)
        return out

    with fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return out

        def find(pattern: str) -> list[int]:
            """Column indices matching `pattern`, excluding host-annotation columns.

            Virome tables carry the phage's bacterial HOST taxonomy next to the
            viral taxonomy (`host_taxonomy`, `host_phylum`, ...). Matching those
            would file every phage whose own taxonomy is blank under Bacteria.
            """
            return [i for i, c in enumerate(header)
                    if c and re.search(pattern, c.strip(), re.I)
                    and not re.search(r"(^|_)host(_|$)", c.strip(), re.I)]

        # a column named Kingdom / Domain / superkingdom, possibly prefixed
        # (silva138_Domain, gtdb_kingdom, ...)
        # Not anchored at the end: real files use domain_note,
        # organism_domain_note, kingdom_group, Kingdom_s ...
        king_idx = find(r"(^|_)(realm|kingdom|domain|superkingdom)(_|$)")
        lineage_idx = find(r"taxonomy|lineage|taxon_path|classification|taxonomy_original|"
                           r"taxonomy_path|full_lineage")
        name_idx = find(r"^(cleaned_name|original_name|feature_id)$")
        rank_idx = find(r"^(taxonomic_rank|feature_level)$")

        # one column index list per Linnaean rank, e.g. 'silva138_Genus' -> genus
        rank_cols = {r: find(rf"(^|_)({r}|{r[0]})$" if r != "kingdom"
                             else r"(^|_)(realm|kingdom|domain|superkingdom)(_|$)")
                     for r in LINEAGE_RANKS}
        out["lineage_cols"] = [r for r in LINEAGE_RANKS if rank_cols[r]]
        phylum_idx, genus_idx, species_idx = (rank_cols["phylum"], rank_cols["genus"],
                                              rank_cols["species"])

        sources_used: Counter = Counter()
        for row in reader:
            out["n_rows"] += 1
            domain = None
            for idx in king_idx:
                if idx < len(row):
                    domain = classify_domain(row[idx])
                    if domain:
                        sources_used["kingdom column"] += 1
                        break
            if not domain:
                for idx in lineage_idx:
                    if idx < len(row):
                        domain = classify_domain(row[idx])
                        if domain:
                            sources_used["lineage string"] += 1
                            break
            if not domain:
                for idx in name_idx:
                    if idx < len(row) and (";" in row[idx] or "__" in row[idx]):
                        domain = classify_domain(row[idx])
                        if domain:
                            sources_used["feature name"] += 1
                            break
            if domain:
                out["n_with_taxonomy"] += 1
            else:
                domain = paper_domain                      # paper-level fallback
                sources_used["paper-level fallback"] += 1
            out["domains"][domain] += 1

            # What this matrix row IS. Only recorded when the paper actually
            # ships a per-feature rank column; otherwise row_ranks stays empty
            # and the caller falls back to the paper-level feature_level text,
            # so there is exactly one answer to this question per paper.
            if rank_idx:
                for i in rank_idx:
                    if i < len(row) and row[i].strip():
                        out["row_ranks"][bucket(row[i], RANK_RULES, "Other / mixed")] += 1
                        break
                else:
                    # the column exists but this row left it blank: still counted,
                    # otherwise the feature silently vanishes from the totals
                    out["row_ranks"]["Other / mixed"] += 1

            # how far up the lineage this row can be rolled up
            for rank_name, idx in rank_cols.items():
                for i in idx:
                    if i < len(row):
                        name = clean_taxon(row[i])
                        if name and not re.match(
                                r"^(unassigned|unclassified|unknown|norank|uncultured|"
                                r"na|none|other|others)$", name, re.I):
                            # a viral name overrides the column it sat in
                            true_rank = ictv_rank(name) or rank_name
                            if true_rank != rank_name:
                                out["ictv_corrected"][f"{rank_name}->{true_rank}"] += 1
                            out["lineage_filled"][true_rank] += 1
                            if true_rank == "phylum":
                                out["phyla"].add(name)
                            elif true_rank == "genus":
                                out["genera"].add(name)
                            elif true_rank == "species":
                                out["species"].add(name)
                            break

        if sources_used:
            out["domain_source"] = sources_used.most_common(1)[0][0]
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=out_dir)
    parser.add_argument("--no-features", action="store_true",
                        help="skip the per-feature scan (~1.2 GB of CSV); "
                             "domain counts then come from the paper-level fallback only")
    args = parser.parse_args()

    dest = args.out_dir
    fig_dir = dest / "figures"
    dest.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ---- cohort -----------------------------------------------------------
    print(f"reading {status_integrated}")
    papers = read_csv(status_integrated)
    excluded = {
        "manual_review": folders_in(manual_review_queue),
        "missing_sample_data": folders_in(missing_sample),
        "missing_feature_data": folders_in(missing_feature),
    }
    # Papers the triage excluded that are nonetheless in active use downstream:
    # the viral dataset is built on all three. Their exclusion reasons concern
    # sample LABELS (no subject/timepoint map; grouping awaiting confirmation),
    # not whether the matrix is usable, so they belong in the cohort - the label
    # gaps then show up honestly in 11_metadata_completeness.csv instead of the
    # papers vanishing from every count.
    RESCUED = {"Garmaeva2024", "Grazioli2022"}   # Pande2026 confirmed 2026-10-14
    drop = set().union(*excluded.values()) - RESCUED
    for name in excluded:
        excluded[name] = excluded[name] - RESCUED
    cohort = [p for p in papers if (p.get("paper_folder") or "").strip() not in drop]

    print(f"  {len(papers)} integrated papers")
    for name, folders in excluded.items():
        print(f"  - {len(folders):>3} excluded ({name})")
    print(f"  + {len(RESCUED):>3} kept despite triage (in use by the viral dataset): "
          f"{', '.join(sorted(RESCUED))}")
    print(f"  = {len(cohort)} papers in the analysis cohort")

    # What the exclusion costs, per paper. Without this the headline domain
    # counts read as database-wide when they are cohort-wide: the triage rules
    # drop papers whose SAMPLES lack labels, which is right for association work
    # but removed Garmaeva2024 - 102,170 VLP viral features - from the viral
    # picture entirely. Anyone reading 03_feature_domain.csv needs this beside it.
    by_folder = {r["paper_folder"].strip(): r for r in papers}
    exc_rows = []
    for why, folders in excluded.items():
        for f in sorted(folders, key=str.lower):
            src = by_folder.get(f, {})
            flev, unit = src.get("feature_level", ""), src.get("main_matrix_unit", "")
            pdom = paper_domain_of(f"{flev} {unit}")
            fm = ({} if args.no_features else scan_features(ROOT / f, pdom))
            dom = fm.get("domains", Counter())
            exc_rows.append([f, why, src.get("n_samples", ""), src.get("n_features", ""),
                             dom.get("Bacteria", 0), dom.get("Viruses", 0),
                             dom.get("Archaea", 0), dom.get("Fungi", 0),
                             dom.get("Eukaryota", 0), unit[:80]])
    exc_rows.sort(key=lambda r: (r[1], -as_int(r[3])))
    write_table(dest / "12_excluded_papers.csv",
                ["paper_folder", "excluded_by", "n_samples", "n_features",
                 "n_features_bacteria", "n_features_viruses", "n_features_archaea",
                 "n_features_fungi", "n_features_eukaryota", "main_matrix_unit"],
                exc_rows)
    exc_viral = sum(as_int(r[5]) for r in exc_rows)
    exc_feat = sum(as_int(r[3]) for r in exc_rows)

    if not cohort:
        print("cohort is empty - nothing to summarise", file=sys.stderr)
        return 1

    # ---- per-paper scan ---------------------------------------------------
    per_paper: list[dict] = []
    tot_samples = 0
    domain_totals: Counter = Counter()
    assay_counter: Counter = Counter()
    unit_counter: Counter = Counter()
    status_counter: Counter = Counter()
    sample_type_counter: Counter = Counter()
    host_counter: Counter = Counter()
    disease_counter: Counter = Counter()
    raw_stypes: Counter = Counter()
    raw_hosts: Counter = Counter()
    phylum_papers: Counter = Counter()
    row_rank_features: Counter = Counter()      # features by matrix-row rank
    row_rank_papers: Counter = Counter()        # papers by dominant matrix-row rank
    lineage_feat: Counter = Counter()           # features carrying each lineage rank
    lineage_paper: Counter = Counter()          # papers carrying each lineage rank
    disease_samples: Counter = Counter()        # samples per named disease
    disease_papers: Counter = Counter()         # papers per named disease
    disease_variants: defaultdict = defaultdict(set)
    ictv_fixes: Counter = Counter()             # viral names re-ranked by suffix
    genus_papers: Counter = Counter()
    species_papers: Counter = Counter()

    for i, p in enumerate(cohort, 1):
        name = (p.get("paper_folder") or "").strip()
        folder = ROOT / name
        n_samples = as_int(p.get("n_samples"))
        n_features = as_int(p.get("n_features"))
        flevel = p.get("feature_level") or ""
        unit = p.get("main_matrix_unit") or ""
        reason = p.get("reason") or ""

        # Only the fields that describe THIS matrix. `reason` is excluded: it
        # describes provenance and alternative tables, so including it made a
        # bacterial paper read as viral (Mirauta2026).
        # Only the fields that describe THIS matrix. `reason` is excluded: it
        # describes provenance and alternative tables, so including it made a
        # bacterial paper read as viral (Mirauta2026).
        paper_domain = paper_domain_of(f"{flevel} {unit}")
        rank = bucket(flevel, RANK_RULES, "Mixed / other")
        unit_b = unit_bucket(unit)

        # paper-level host hint: the per-sample metadata often omits `host`
        # entirely (only 89 of 265 papers have the column), but the title and
        # the integration note almost always name the organism.
        paper_host_hint = bucket(f"{p.get('title', '')} {reason} {flevel}", HOST_RULES, "")
        if paper_host_hint == "Not recorded":
            paper_host_hint = ""
        sm = scan_samples(folder, paper_host_hint)
        assay = (sm["methods"].most_common(1)[0][0] if sm["methods"]
                 else bucket(flevel, ASSAY_RULES))

        if args.no_features:
            fm = {"n_rows": 0, "domains": Counter({paper_domain: n_features}),
                  "domain_source": "paper-level fallback", "phyla": set(),
                  "genera": set(), "species": set(), "n_with_taxonomy": 0,
                  "ictv_corrected": Counter(),
                  "row_ranks": Counter({rank: n_features}),
                  "lineage_filled": Counter(), "lineage_cols": []}
        else:
            fm = scan_features(folder, paper_domain)
            if fm["n_rows"] == 0:                          # no feature_metadata on disk
                fm["domains"] = Counter({paper_domain: n_features})
                fm["domain_source"] = "paper-level fallback"

        print(f"  [{i:>3}/{len(cohort)}] {name:<28} "
              f"{n_samples:>6} samples x {n_features:>7} features", flush=True)

        tot_samples += n_samples
        domain_totals.update(fm["domains"])
        assay_counter[assay] += 1
        unit_counter[unit_b] += 1
        status_counter[p.get("status") or "unknown"] += 1
        sample_type_counter.update(sm["sample_types"])
        host_counter.update(sm["hosts"])
        disease_counter.update(sm["disease"])
        raw_stypes.update(sm["raw_sample_types"])
        raw_hosts.update(sm["raw_hosts"])
        ictv_fixes.update(fm["ictv_corrected"])
        disease_samples.update(sm["diseases"])
        for dname in sm["diseases"]:
            disease_papers[dname] += 1
            disease_variants[dname] |= sm["disease_variants"][dname]
        # Papers without a per-feature rank column fall back to the paper-level
        # rank for BOTH columns, so n_features still sums to the database total.
        if fm["row_ranks"]:
            row_rank_features.update(fm["row_ranks"])
        else:
            row_rank_features[rank] += n_features
        lineage_feat.update(fm["lineage_filled"])
        for r in fm["lineage_cols"]:
            lineage_paper[r] += 1
        matrix_row_rank = (fm["row_ranks"].most_common(1)[0][0] if fm["row_ranks"] else rank)
        row_rank_papers[matrix_row_rank] += 1
        deepest = ([r for r in LINEAGE_RANKS if fm["lineage_filled"].get(r)] or [""])[-1]
        for taxon in fm["phyla"]:
            phylum_papers[taxon] += 1
        for taxon in fm["genera"]:
            genus_papers[taxon] += 1
        for taxon in fm["species"]:
            species_papers[taxon] += 1

        dominant_type = (sm["sample_types"].most_common(1)[0][0]
                         if sm["sample_types"] else "Not recorded")
        dominant_host = (sm["hosts"].most_common(1)[0][0] if sm["hosts"] else "Not recorded")

        per_paper.append({
            "paper_folder": name,
            "doi": p.get("doi", ""),
            "title": (p.get("title") or "")[:200],
            "status": p.get("status", ""),
            "n_samples": n_samples,
            "n_features": n_features,
            "n_cells": n_samples * n_features,
            "feature_rank": rank,
            "matrix_row_rank": matrix_row_rank,
            "lineage_ranks_available": "|".join(fm["lineage_cols"]),
            "deepest_lineage_rank": deepest,
            "can_rollup_to_genus": "genus" in fm["lineage_cols"],
            "can_rollup_to_phylum": "phylum" in fm["lineage_cols"],
            "assay": assay,
            "abundance_unit_bucket": unit_b,
            "main_matrix_unit": unit,
            "is_relative": unit_b.startswith("Relative"),
            "is_count": unit_b.startswith("Read count"),
            "dominant_sample_type": dominant_type,
            "dominant_host": dominant_host,
            "has_disease_info": p.get("has_disease_info", ""),
            "diseases_named": "; ".join(sorted(sm["diseases"])) [:300],
            "n_samples_disease_labelled": sm["n_disease_labelled"],
            "n_samples_group_labelled": sm["n_group_labelled"],
            "n_subjects_distinct": sm["n_subjects_distinct"],
            "is_longitudinal": sm["longitudinal"],
            "n_samples_with_sex": sm["n_sex"],
            "n_samples_with_age": sm["n_age"],
            "n_samples_with_timepoint": sm["n_timepoint"],
            "n_features_bacteria": fm["domains"].get("Bacteria", 0),
            "n_features_archaea": fm["domains"].get("Archaea", 0),
            "n_features_viruses": fm["domains"].get("Viruses", 0),
            "n_features_fungi": fm["domains"].get("Fungi", 0),
            "n_features_eukaryota": fm["domains"].get("Eukaryota", 0),
            "n_features_unclassified": fm["domains"].get("Unclassified", 0),
            "domain_evidence": fm["domain_source"],
            "n_viral_names_rerankeed_by_ictv": sum(fm["ictv_corrected"].values()),
            "n_features_with_real_taxonomy": fm["n_with_taxonomy"],
            "share_features_with_real_taxonomy": pct(fm["n_with_taxonomy"],
                                                     fm["n_rows"] or n_features),
            "n_distinct_phyla": len(fm["phyla"]),
            "n_distinct_genera": len(fm["genera"]),
            "n_distinct_species": len(fm["species"]),
            "has_raw_accession": bool((p.get("raw_data_accessions") or "").strip()),
            "n_samples_label_issue": as_int(p.get("n_samples_label_issue")),
            "n_cells_ambiguous_missing": as_int(p.get("n_cells_ambiguous_missing")),
        })

    # ---- tables -----------------------------------------------------------
    print("\nwriting tables ...")
    fields = list(per_paper[0].keys())
    with (dest / "01_papers_included.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(per_paper)

    tot_features = sum(r["n_features"] for r in per_paper)
    tot_cells = sum(r["n_cells"] for r in per_paper)
    n_disease_papers = sum(1 for r in per_paper if str(r["has_disease_info"]).lower() == "true")
    n_longitudinal = sum(1 for r in per_paper if r["is_longitudinal"])
    n_relative = sum(1 for r in per_paper if r["is_relative"])
    n_count = sum(1 for r in per_paper if r["is_count"])
    n_with_accession = sum(1 for r in per_paper if r["has_raw_accession"])
    samples_disease_labelled = sum(r["n_samples_disease_labelled"] for r in per_paper)
    samples_group_labelled = sum(r["n_samples_group_labelled"] for r in per_paper)
    subjects = sum(r["n_subjects_distinct"] for r in per_paper)

    tot_dis = sum(disease_samples.values())
    headline = [
        ["papers in cohort", len(per_paper)],
        ["papers excluded (3 triage files)", len(drop)],
        ["papers kept despite triage (used downstream)", len(RESCUED)],
        ["features in excluded papers (NOT in any table below)", exc_feat],
        ["viral features in excluded papers", exc_viral],
        ["total samples", tot_samples],
        ["total features (sum over papers)", tot_features],
        ["total matrix cells", tot_cells],
        ["distinct subjects (sum over papers)", subjects],
        ["median samples per paper", int(median([r["n_samples"] for r in per_paper]))],
        ["median features per paper", int(median([r["n_features"] for r in per_paper]))],
        ["largest study (samples)", max(r["n_samples"] for r in per_paper)],
        ["largest study (features)", max(r["n_features"] for r in per_paper)],
        ["papers flagged has_disease_info", n_disease_papers],
        ["papers with longitudinal subjects", n_longitudinal],
        ["papers with relative-abundance matrix", n_relative],
        ["papers with count matrix", n_count],
        ["papers with a raw-read accession", n_with_accession],
        ["samples carrying a real disease label", samples_disease_labelled],
        ["samples carrying an experimental_group", samples_group_labelled],
        ["bacterial features", domain_totals.get("Bacteria", 0)],
        ["archaeal features", domain_totals.get("Archaea", 0)],
        ["viral features", domain_totals.get("Viruses", 0)],
        ["fungal features", domain_totals.get("Fungi", 0)],
        ["other eukaryotic features", domain_totals.get("Eukaryota", 0)],
        ["unclassified features", domain_totals.get("Unclassified", 0)],
        ["distinct named diseases", len(disease_samples)],
        ["samples with a named disease", tot_dis],
        ["features with per-feature taxonomy", sum(
            r["n_features_with_real_taxonomy"] for r in per_paper)],
        ["features whose domain is a paper-level assertion", sum(
            r["n_features"] - r["n_features_with_real_taxonomy"] for r in per_paper)],
        ["papers whose matrix can roll up to genus", sum(
            1 for r in per_paper if r["can_rollup_to_genus"])],
        ["papers whose matrix can roll up to phylum", sum(
            1 for r in per_paper if r["can_rollup_to_phylum"])],
        ["distinct phylum names across DB", len(phylum_papers)],
        ["distinct genus names across DB", len(genus_papers)],
        ["distinct species names across DB", len(species_papers)],
    ]
    write_table(dest / "02_summary_headline.csv", ["metric", "value"], headline)

    domain_items = counter_table(dest / "03_feature_domain.csv", domain_totals,
                                 "domain", "n_features")
    stype_items = counter_table(dest / "05_sample_type.csv", sample_type_counter,
                                "sample_type", "n_samples")
    disease_items = counter_table(dest / "06_disease_status.csv", disease_counter,
                                  "disease_status", "n_samples")
    host_items = counter_table(dest / "07_host.csv", host_counter, "host", "n_samples")
    assay_items = counter_table(dest / "08_assay.csv", assay_counter, "assay", "n_papers")
    unit_items = counter_table(dest / "09_abundance_unit.csv", unit_counter,
                               "abundance_unit", "n_papers")
    counter_table(dest / "10_paper_status.csv", status_counter, "status", "n_papers")

    n = len(per_paper)
    write_table(dest / "11_metadata_completeness.csv",
                ["field", "n_papers_with_any", "share_of_papers",
                 "n_samples_with_value", "share_of_samples"],
                [[label,
                  sum(1 for r in per_paper if r[key] > 0),
                  pct(sum(1 for r in per_paper if r[key] > 0), n),
                  sum(r[key] for r in per_paper),
                  pct(sum(r[key] for r in per_paper), tot_samples)]
                 for label, key in [
                     ("disease label", "n_samples_disease_labelled"),
                     ("experimental_group", "n_samples_group_labelled"),
                     ("timepoint", "n_samples_with_timepoint"),
                     ("sex", "n_samples_with_sex"),
                     ("age", "n_samples_with_age"),
                 ]])

    # 16: what the matrix rows ARE, at both feature and paper level. Ordered
    # coarse -> fine rather than by size, so the resolution gradient is readable.
    tot_rows = sum(row_rank_features.values())
    rank_rows = [[r, row_rank_features.get(r, 0), pct(row_rank_features.get(r, 0), tot_rows),
                  row_rank_papers.get(r, 0), pct(row_rank_papers.get(r, 0), len(per_paper))]
                 for r in MATRIX_ROW_RANKS
                 if row_rank_features.get(r, 0) or row_rank_papers.get(r, 0)]
    write_table(dest / "16_matrix_row_rank.csv",
                ["matrix_row_rank", "n_features", "share_of_features",
                 "n_papers_dominant", "share_of_papers"], rank_rows)

    # 17: how far up each rank can be rolled up. This is the cross-study merge
    # feasibility table: a rank is usable only for the papers that carry it.
    write_table(dest / "17_lineage_depth.csv",
                ["lineage_rank", "n_papers_with_column", "share_of_papers",
                 "n_features_with_value", "share_of_features"],
                [[r, lineage_paper.get(r, 0), pct(lineage_paper.get(r, 0), len(per_paper)),
                  lineage_feat.get(r, 0), pct(lineage_feat.get(r, 0), tot_rows)]
                 for r in LINEAGE_RANKS])

    # 18: which diseases, with the raw spellings that were merged kept visible
    write_table(dest / "18_disease_catalog.csv",
                ["disease", "category", "n_samples", "share_of_disease_samples",
                 "n_papers", "raw_variants_merged"],
                [[d, bucket(d, DISEASE_CATEGORY_RULES, "Other / unspecified"), n,
                  pct(n, tot_dis), disease_papers[d],
                  " | ".join(sorted(disease_variants[d]))]
                 for d, n in disease_samples.most_common()])

    dis_cat: Counter = Counter()
    for d, n in disease_samples.items():
        dis_cat[bucket(d, DISEASE_CATEGORY_RULES, "Other / unspecified")] += n
    dis_cat_items = counter_table(dest / "19_disease_category.csv", dis_cat,
                                  "disease_category", "n_samples")

    for fname, counter, key in [("13_taxon_catalog_phylum.csv", phylum_papers, "phylum"),
                                ("14_taxon_catalog_genus.csv", genus_papers, "genus"),
                                ("15_taxon_catalog_species.csv", species_papers, "species")]:
        write_table(dest / fname, [key, "n_papers_containing_it"], counter.most_common())

    write_table(dest / "20_ictv_rank_corrections.csv",
                ["column_rank -> true_rank (by ICTV suffix)", "n_names"],
                ictv_fixes.most_common())

    write_table(dest / "99_unmapped_values.csv", ["field", "raw_value", "n_samples"],
                [["sample_type", k, v] for k, v in raw_stypes.most_common(400)] +
                [["host", k, v] for k, v in raw_hosts.most_common(400)])

    # ---- figures ----------------------------------------------------------
    print("writing figures ...")
    pies = [
        ("fig01_feature_domain_pie", "Features by domain", domain_items),
        ("fig02_assay_pie", "Papers by assay type", assay_items),
        ("fig03_unit_pie", "Papers by abundance unit", unit_items),
        ("fig04_disease_pie", "Samples by disease status", disease_items),
    ]
    bars = [
        ("fig05_sample_type_bar", "Samples by sample type", stype_items),
        ("fig06_host_bar", "Samples by host", host_items),
        ("fig08_top_genera_bar", "Most widely shared genera (papers containing them)",
         genus_papers.most_common(18)),
        ("fig10_matrix_row_rank_bar", "Features by matrix row rank (coarse -> fine)",
         [(r[0], r[1]) for r in rank_rows]),
        ("fig11_lineage_depth_bar", "Papers carrying each lineage rank (roll-up depth)",
         [(r, lineage_paper.get(r, 0)) for r in LINEAGE_RANKS]),
        ("fig12_top_diseases_bar", "Samples by named disease",
         disease_samples.most_common(18)),
        ("fig13_disease_category_bar", "Samples by disease category", dis_cat_items),
    ]
    for name, title, items in pies:
        svg_pie(fig_dir / f"{name}.svg", title, items,
                subtitle=f"{len(per_paper)} papers, {tot_samples:,} samples")
    for name, title, items in bars:
        svg_barh(fig_dir / f"{name}.svg", title, items,
                 subtitle=f"{len(per_paper)} papers, {tot_samples:,} samples")
    svg_scatter(fig_dir / "fig09_size_scatter.svg", "Study size: samples vs features",
                [(r["n_samples"], r["n_features"], r["paper_folder"]) for r in per_paper],
                "samples per paper", "features per paper",
                subtitle="one dot = one paper (log-log)")
    maybe_png(fig_dir, pies, bars)

    # ---- report -----------------------------------------------------------
    lines = [
        "# Abundance database - summary statistics",
        "",
        f"Cohort: **{len(per_paper)} papers** "
        f"({len(papers)} integrated minus {len(drop)} parked in the three triage files).",
        "",
        "| metric | value |",
        "|---|---:|",
    ]
    lines += [f"| {k} | {v:,} |" if isinstance(v, int) else f"| {k} | {v} |"
              for k, v in headline]

    def section(title: str, items, key_name: str, total: int) -> list[str]:
        rows = [f"", f"## {title}", "", f"| {key_name} | n | share |", "|---|---:|---:|"]
        rows += [f"| {k} | {v:,} | {pct(v, total)} |" for k, v in items]
        return rows

    lines += section("Features by domain", domain_items, "domain", sum(domain_totals.values()))
    lines += section("Papers by assay", assay_items, "assay", len(per_paper))
    lines += section("Papers by abundance unit", unit_items, "unit", len(per_paper))
    lines += ["", "## Matrix row rank (what one row of the matrix IS)", "",
              "| rank | features | share | papers (dominant) |", "|---|---:|---:|---:|"]
    lines += [f"| {r[0]} | {r[1]:,} | {r[2]} | {r[3]} |" for r in rank_rows]
    lines += ["", "## Roll-up depth (how far up a matrix can be aggregated)", "",
              "A rank is usable for cross-study work only in the papers that ship "
              "that lineage column. Aggregation only ever goes coarser: a "
              "genus-level matrix can never be split back into species.", "",
              "| lineage rank | papers with the column | features with a value |",
              "|---|---:|---:|"]
    lines += [f"| {r} | {lineage_paper.get(r, 0)} | {lineage_feat.get(r, 0):,} |"
              for r in LINEAGE_RANKS]
    lines += section("Samples by sample type", stype_items, "sample type", tot_samples)
    lines += section("Samples by host", host_items, "host", tot_samples)
    lines += section("Samples by disease status", disease_items, "disease status", tot_samples)
    lines += section("Samples by disease category", dis_cat_items, "category", tot_dis)
    lines += ["", "## Named diseases (top 25)", "",
              "| disease | category | samples | papers |", "|---|---|---:|---:|"]
    lines += [f"| {d} | {bucket(d, DISEASE_CATEGORY_RULES, 'Other / unspecified')} "
              f"| {n:,} | {disease_papers[d]} |"
              for d, n in disease_samples.most_common(25)]
    lines += [
        "",
        "## Caveats",
        "",
        "- `sample_type` / `host` / `disease_status` are free text in the source "
        "metadata; the buckets come from the keyword rules at the top of this "
        "script. Values that matched nothing are in `99_unmapped_values.csv` - "
        "extend the rules there and re-run.",
        "- Feature domains come from per-feature taxonomy where the paper shipped "
        "it, and from the matrix's own unit description otherwise. "
        "`domain_evidence` in `01_papers_included.csv` says which, per paper, and "
        "`share_features_with_real_taxonomy` says how much of each paper rests on "
        "per-feature evidence.",
        "- A per-feature cell naming several domains at once "
        "(\"bacteria/archaea\") resolves nothing for that feature, so it falls "
        "back to the paper level; a matrix description naming three or more "
        "domains resolves nothing either and lands in Unclassified. Neither is "
        "collapsed to whichever domain a rule happens to list first.",
        "- Viral names are re-ranked by ICTV suffix wherever a file filed them "
        "under the wrong column; `20_ictv_rank_corrections.csv` is the audit "
        "trail.",
        "- Every table here describes the COHORT, not the whole database. The "
        "39 excluded papers are listed with their own feature counts in "
        "`12_excluded_papers.csv`; the viral totals in particular are much "
        "larger once they are counted, because one excluded paper is a "
        "102,170-feature VLP virome whose samples simply lack labels.",
        "- Feature counts are summed across papers, so a genus present in 50 "
        "papers is counted 50 times. For the de-duplicated view use the taxon "
        "catalogues (`13`-`15`).",
    ]
    (dest / "00_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\ndone -> {dest}")
    print(f"  {len(per_paper)} papers | {tot_samples:,} samples | {tot_features:,} features")
    print(f"  bacteria {domain_totals.get('Bacteria', 0):,} | "
          f"viruses {domain_totals.get('Viruses', 0):,} | "
          f"archaea {domain_totals.get('Archaea', 0):,} | "
          f"fungi {domain_totals.get('Fungi', 0):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
