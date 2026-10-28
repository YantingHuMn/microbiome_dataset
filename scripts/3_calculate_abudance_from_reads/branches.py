#!/usr/bin/env python3
"""What kind of sequencing a library is, and which analysis it may enter.

Imported by step 2 (which classifies from ENA metadata) and step 3 (which
refines the amplicon marker from the reads themselves), so the two can never
disagree about what a branch name means.

Two things are kept apart on purpose:

    the ORGANISM a paper studies   - a paper about bacteria may hold ITS reads,
                                     a paper about viruses may hold bulk shotgun
    the ASSAY the reads came from  - which is what decides the analysis

`AMPLICON` in ENA does not mean 16S: it means a PCR product of something. `WGS`
does not mean an analysable microbial metagenome: it is also how an isolate
genome or a human whole-genome library is deposited. Both are resolved further
before any branch is allowed to run, and anything still ambiguous blocks.
"""

# branch -> (what it is, what step 3 does with it)
BRANCHES = {
    "amplicon_16s": ("16S/16S-like rRNA amplicon (bacteria + archaea)",
                     "DADA2 ASVs + SILVA, step3_16s.R"),
    "amplicon_16s_single": ("16S amplicon sequenced single-end",
                            "DADA2 on R1 only, no merge, step3_16s_single.R"),
    "amplicon_its": ("ITS amplicon (fungi)", "not implemented"),
    "amplicon_18s": ("18S rRNA amplicon (eukaryotes)", "not implemented"),
    "amplicon_other": ("an amplicon of some other marker", "not implemented"),
    "dna_shotgun": ("bulk DNA shotgun metagenome",
                    "sylph vs GTDB for bacteria/archaea; optionally also the "
                    "viral branch, step3_shotgun.py"),
    "dna_virome": ("VLP-enriched / virus-targeted DNA library",
                   "read mapping to a viral catalogue, step3_virome.py"),
    "rna": ("RNA-seq / metatranscriptome", "not implemented"),
    "unsupported_isolate_genome": ("a single isolate or host genome, not a community",
                                   "nothing; this is not a microbiome sample"),
    "unresolved": ("could not be established from the available evidence",
                   "nothing; blocked until resolved"),
}

#: Branches step 3 can actually run today.
IMPLEMENTED = {"amplicon_16s", "amplicon_16s_single", "dna_shotgun",
               "dna_virome"}

#: Amplicon markers, resolved from the reads in step 3 rather than assumed.
MARKER_TO_BRANCH = {
    "16S_V1V3": "amplicon_16s", "16S_V3V4": "amplicon_16s", "16S_V4": "amplicon_16s",
    "ITS1": "amplicon_its", "ITS2": "amplicon_its",
    "18S_V4": "amplicon_18s", "18S_V9": "amplicon_18s",
}

# Wording that marks a library as virus-targeted rather than bulk. ENA has no
# field for "we filtered to VLPs", so this is a hint, never a conclusion: the
# paper facts decide, and an unsupported guess blocks instead of defaulting.
VIRAL_SELECTION_HINTS = ("VIRAL", "SIZE FRACTIONATION", "VIRUS")


def branch_from_ena(strategy, source, selection=None):
    """Classify from ENA fields alone. Returns (branch, why).

    Deliberately coarse: the amplicon marker and VLP enrichment cannot be read
    out of these three fields, so they resolve to `amplicon` and to a hint that
    step 3 must confirm.
    """
    s = (strategy or "").upper().strip()
    src = (source or "").upper().strip()
    sel = (selection or "").upper().strip()

    if s == "AMPLICON" or "AMPLICON" in sel:
        # The marker is NOT decided here. 16S, ITS and 18S are all AMPLICON.
        return "amplicon", (f"ENA library_strategy={s or '-'}, selection={sel or '-'}: "
                            "a PCR product, marker still unknown")
    if (s == "RNA-SEQ" or "TRANSCRIPTOMIC" in src or "RNA" in src.split()):
        # Checked before the WGS rules: "VIRAL RNA" is an RNA library whatever
        # the strategy field says, and must not land in a DNA branch.
        return "rna", f"ENA library_source={src}, strategy={s}: an RNA library"
    if s in ("WGS", "WGA"):
        if "VIRAL" in src or any(h in sel for h in VIRAL_SELECTION_HINTS):
            return "dna_virome", (f"ENA library_source={src}, selection={sel or '-'}: "
                                  "virus-targeted")
        if "METAGENOMIC" in src:
            return "dna_shotgun", f"ENA library_source={src}, strategy={s}"
        if "GENOMIC" in src:
            # A single organism's genome, not a community. Counting taxa in it
            # would produce a matrix that means nothing.
            return "unsupported_isolate_genome", (
                f"ENA library_source={src}: an isolate or host genome, not a community")
    return "unresolved", (f"ENA strategy={s or '-'}, source={src or '-'}, "
                          f"selection={sel or '-'} do not identify the assay")


def refine_amplicon(marker):
    """Map a marker measured from the reads to its branch."""
    if not marker:
        return "unresolved"
    if marker.startswith("primer_like_"):
        marker = marker[len("primer_like_"):]
    return MARKER_TO_BRANCH.get(marker, "amplicon_other")


def describe(branch):
    what, does = BRANCHES.get(branch, ("unknown branch", "nothing"))
    return f"{branch}: {what} -> {does}"
