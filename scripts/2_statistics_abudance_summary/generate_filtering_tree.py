#!/usr/bin/env python3
"""Draw the paper/sample/feature filtering tree for abundance integration.

The diagram is calculated from the current CSV files, never from hard-coded
counts.  Samples and features are reported once a paper has an integrated
matrix; before that point they are shown as N/A because these quantities do not
exist consistently for papers with raw reads, group summaries, or no data.

Outputs
-------
database_filtering_tree.svg
database_filtering_tree.png      (when matplotlib is installed)

Usage
-----
python generate_filtering_tree.py
python generate_filtering_tree.py --out-dir /path/to/output
"""

from __future__ import annotations

import argparse
import csv
import html
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


DEFAULT_ROOT = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/data/abudance_integration"
)
DEFAULT_OUT = Path(
    "/hickory/proj/didonglab/dataset/virus/Database/results/"
    "2_statistics_abudance_summary/figures"
)


@dataclass
class Counts:
    papers: int
    samples: int | None
    features: int | None


@dataclass
class Node:
    node_id: str
    parent_id: str
    label: str
    kind: str
    counts: Counts
    note: str = ""


# Why a paper produced no matrix. The category is the status FILE a paper is
# filed in, not the `status` column inside the row: rows were re-filed by hand
# and their original triage verdict was deliberately left untouched, so the
# column is a record of the first pass and the file is the current answer.
# Reading the column instead would, for example, show 10 pipeline errors where
# 31 papers are actually filed under one.
NONINTEGRATED_FILES = [
    ("1_2_status_raw_reads_only.csv", "Only raw reads"),
    ("1_7_status_download_failed.csv", "Download failed"),
    ("1_6_status_pipeline_error.csv", "Pipeline error"),
    ("1_3_status_unsuitable_for_abundance_integration.csv",
     "Unsuitable for abundance integration"),
    ("1_5_status_group_level_only.csv", "Only group-level summaries"),
    ("1_4_status_controlled_access.csv", "Controlled access"),
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(value: str | None) -> int:
    try:
        return int(float((value or "").strip()))
    except ValueError:
        return 0


def folders(rows: list[dict[str, str]]) -> set[str]:
    return {(r.get("paper_folder") or "").strip() for r in rows}


def sum_counts(rows: list[dict[str, str]]) -> Counts:
    return Counts(
        papers=len(rows),
        samples=sum(number(r.get("n_samples")) for r in rows),
        features=sum(number(r.get("n_features")) for r in rows),
    )


def subset(rows: list[dict[str, str]], names: set[str]) -> list[dict[str, str]]:
    return [r for r in rows if (r.get("paper_folder") or "").strip() in names]


def assert_unique_papers(rows: list[dict[str, str]], source: str) -> None:
    names = [(r.get("paper_folder") or "").strip() for r in rows]
    blank = sum(not name for name in names)
    duplicates = sorted({name for name in names if name and names.count(name) > 1})
    if blank or duplicates:
        raise ValueError(
            f"{source}: paper_folder must be non-empty and unique; "
            f"blank={blank}, duplicates={duplicates}"
        )


def assert_partition(parent: Counts, kept: Counts, removed: Counts, stage: str) -> None:
    """Require exact paper/sample/feature conservation at one filtering step."""
    problems = []
    for field in ("papers", "samples", "features"):
        p = getattr(parent, field)
        k = getattr(kept, field)
        r = getattr(removed, field)
        if p is None or k is None or r is None:
            problems.append(f"{field}: cannot verify because one value is N/A")
        elif p != k + r:
            problems.append(f"{field}: parent={p:,}, kept={k:,}, removed={r:,}")
    if problems:
        raise ValueError(
            f"Filtering hierarchy failed strict conservation at '{stage}': "
            + "; ".join(problems)
        )


def build_nodes(root: Path) -> list[Node]:
    overview = read_rows(root / "1_overview.csv")
    integrated = read_rows(root / "1_1_status_integrated.csv")
    assert_unique_papers(overview, "1_overview.csv")
    assert_unique_papers(integrated, "1_1_status_integrated.csv")
    integrated_names = folders(integrated)
    overview_integrated_names = {
        (r.get("paper_folder") or "").strip()
        for r in overview
        if (r.get("status") or "").strip().startswith("integrated")
    }
    if integrated_names != overview_integrated_names:
        raise ValueError(
            "Integrated-paper membership differs between 1_overview.csv and "
            "1_1_status_integrated.csv; only_in_overview="
            f"{sorted(overview_integrated_names - integrated_names)}, "
            f"only_in_status_integrated={sorted(integrated_names - overview_integrated_names)}"
        )

    triage_specs = [
        ("manual", "1_1_1_manual_review_queue.csv", "Manual review required",
         "Sample/group identity or another mapping needs a human decision"),
        ("sample", "1_1_2_status_missing_sample_data.csv", "Missing sample information",
         "Matrix exists, but usable per-sample metadata are missing"),
        ("feature", "1_1_3_status_missing_feature_data.csv", "Missing feature information",
         "Features are opaque IDs: neither taxonomy nor sequence is available"),
    ]
    retained_labels = {
        "manual": "After manual-review exclusion",
        "sample": "After missing-sample-info exclusion",
        "feature": "Final: after missing-feature-info exclusion",
    }
    triage: list[tuple[str, str, str, set[str]]] = []
    seen: set[str] = set()
    for key, filename, label, note in triage_specs:
        path = root / filename
        queue_rows = read_rows(path) if path.exists() else []
        assert_unique_papers(queue_rows, filename)
        names = folders(queue_rows)
        outside = names - integrated_names
        if outside:
            raise ValueError(f"{filename} contains non-integrated papers: {sorted(outside)}")
        overlap = names & seen
        if overlap:
            raise ValueError(
                f"Triage queues overlap, so a sequential tree is ambiguous: {sorted(overlap)}"
            )
        seen |= names
        triage.append((key, label, note, names))

    nodes = [
        Node("all", "", "All candidate papers", "root",
             Counts(len(overview), None, None),
             "Sample/feature totals are N/A before a sample-level matrix exists"),
        Node("integrated", "all", "Sample-level abundance matrix produced", "keep",
             sum_counts(integrated)),
    ]

    nonintegrated_names = {
        (r.get("paper_folder") or "").strip() for r in overview
    } - integrated_names
    status_groups: list[tuple[str, set[str]]] = []
    filed: set[str] = set()
    for filename, label in NONINTEGRATED_FILES:
        path = root / filename
        rows = read_rows(path) if path.exists() else []
        assert_unique_papers(rows, filename)
        names = folders(rows)
        stray = names & integrated_names
        if stray:
            raise ValueError(f"{filename} lists integrated papers: {sorted(stray)}")
        clash = names & filed
        if clash:
            raise ValueError(
                f"{filename} shares papers with an earlier status file, so the "
                f"categories are not a partition: {sorted(clash)}"
            )
        filed |= names
        status_groups.append((label, names))
    # Anything not filed anywhere is shown as exactly that, never folded into a
    # neighbouring category to make the arithmetic work.
    unfiled = nonintegrated_names - filed
    if unfiled:
        status_groups.append(("Not yet categorised", unfiled))
    nonintegrated_n = sum(len(v) for _, v in status_groups)
    if nonintegrated_n != len(nonintegrated_names):
        raise ValueError(
            "Non-integrated branches do not sum to the overview total: "
            f"branches={nonintegrated_n}, expected={len(nonintegrated_names)}"
        )
    nodes.append(Node(
        "not_integrated", "all", "No sample-level abundance matrix", "drop",
        Counts(nonintegrated_n, None, None),
        "Samples/features are not comparable or unavailable at this stage",
    ))
    for label, names in sorted(status_groups, key=lambda x: (-len(x[1]), x[0])):
        nodes.append(Node(
            f"status_{label.lower().replace(' ', '_')}", "not_integrated", label,
            "detail", Counts(len(names), None, None),
        ))

    remaining = integrated
    remaining_counts = sum_counts(remaining)
    parent = "integrated"
    for i, (key, label, note, names) in enumerate(triage, start=1):
        removed = subset(remaining, names)
        removed_counts = sum_counts(removed)
        nodes.append(Node(f"drop_{key}", parent, label, "drop", removed_counts, note))
        remaining = [r for r in remaining if (r.get("paper_folder") or "").strip() not in names]
        kept_counts = sum_counts(remaining)
        assert_partition(remaining_counts, kept_counts, removed_counts, label)
        keep_id = f"after_{key}"
        nodes.append(Node(
            keep_id, parent, retained_labels[key], "keep", kept_counts
        ))
        parent = keep_id
        remaining_counts = kept_counts

    if len(overview) != len(integrated) + nonintegrated_n:
        raise ValueError("overview totals do not reconcile with integrated/non-integrated branches")
    return nodes


def fmt_int(value: int | None) -> str:
    return "N/A" if value is None else f"{value:,}"


def node_lines(node: Node) -> list[str]:
    c = node.counts
    return [node.label, f"Papers: {c.papers:,}",
            f"Samples: {fmt_int(c.samples)}", f"Features: {fmt_int(c.features)}"]


COLORS = {
    "root": ("#17324D", "#FFFFFF"),
    "keep": ("#DDF3E4", "#153C27"),
    "drop": ("#FCE3E1", "#70251F"),
    "detail": ("#FFF4D6", "#5D4611"),
}


def layout(nodes: list[Node]) -> tuple[dict[str, tuple[float, float]], float, float]:
    # Three lanes: retained cohort (left), sequential exclusions (middle), and
    # reasons no matrix was produced (right).  Keeping the latter two separate
    # prevents an availability reason from looking like a post-integration QC
    # filter.
    pos: dict[str, tuple[float, float]] = {
        "all": (700, 70), "integrated": (350, 245), "not_integrated": (1450, 245),
    }
    details = [n for n in nodes if n.parent_id == "not_integrated"]
    for i, n in enumerate(details):
        col, row = i % 2, i // 2
        pos[n.node_id] = (1250 + col * 400, 430 + row * 145)

    y = 430
    current = "integrated"
    for key in ("manual", "sample", "feature"):
        pos[f"drop_{key}"] = (800, y)
        pos[f"after_{key}"] = (360, y + 180)
        current = f"after_{key}"
        y += 180
    detail_bottom = 430 + ((len(details) + 1) // 2) * 145
    return pos, 1850, max(y + 160, detail_bottom + 80)


def write_svg(nodes: list[Node], path: Path) -> None:
    pos, width, height = layout(nodes)
    box_w, box_h = 330, 112

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FAFCFE"/>',
        '<style>text{font-family:Arial,"Noto Sans",sans-serif}.title{font-size:24px;font-weight:700}'
        '.label{font-size:15px;font-weight:700}.count{font-size:13px}.edge{stroke:#728197;stroke-width:2;fill:none}'
        '.dash{stroke-dasharray:7 5}</style>',
        '<text x="40" y="34" class="title" fill="#17324D">Abundance database filtering flow</text>',
    ]
    for n in nodes:
        if not n.parent_id:
            continue
        x, y = pos[n.node_id]
        px, py = pos[n.parent_id]
        dash = ' class="edge dash"' if n.kind in {"drop", "detail"} else ' class="edge"'
        out.append(
            f'<path{dash} d="M {px} {py + box_h/2} C {px} {(py+y)/2}, {x} {(py+y)/2}, {x} {y-box_h/2}"/>'
        )

    for n in nodes:
        x, y = pos[n.node_id]
        fill, text_color = COLORS[n.kind]
        out.append(
            f'<rect x="{x-box_w/2}" y="{y-box_h/2}" width="{box_w}" height="{box_h}" '
            f'rx="13" fill="{fill}" stroke="{text_color}" stroke-width="1.6"/>'
        )
        lines = node_lines(n)
        for i, line in enumerate(lines):
            cls = "label" if i == 0 else "count"
            yy = y - 30 + i * 23
            out.append(
                f'<text x="{x}" y="{yy}" text-anchor="middle" class="{cls}" '
                f'fill="{text_color}">{html.escape(line)}</text>'
            )

    out += [
        f'<text x="40" y="{height-30}" font-size="12" fill="#566474">'
        'Solid path = retained cohort; dashed branches = exclusions or reason breakdown. '
        'Counts are recalculated from the current CSV files.</text>',
        '</svg>',
    ]
    path.write_text("\n".join(out), encoding="utf-8")


def write_png(nodes: list[Node], path: Path) -> bool:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
        from matplotlib.path import Path as MplPath
        from matplotlib.patches import PathPatch
    except ImportError:
        # The project environment may not have matplotlib.  librsvg is a
        # lightweight and faithful fallback because the SVG is always written
        # first by main().
        converter = shutil.which("rsvg-convert")
        svg_path = path.with_suffix(".svg")
        if converter and svg_path.exists():
            subprocess.run(
                [converter, str(svg_path), "-o", str(path)],
                check=True,
            )
            return True
        return False

    pos, width, height = layout(nodes)
    box_w, box_h = 330, 112
    fig, ax = plt.subplots(figsize=(width / 120, height / 120), dpi=180)
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")
    fig.patch.set_facecolor("#FAFCFE")
    ax.set_facecolor("#FAFCFE")
    ax.text(40, 34, "Abundance database filtering flow", fontsize=17,
            fontweight="bold", color="#17324D")

    for n in nodes:
        if not n.parent_id:
            continue
        x, y = pos[n.node_id]
        px, py = pos[n.parent_id]
        verts = [(px, py + box_h/2), (px, (py+y)/2), (x, (py+y)/2), (x, y-box_h/2)]
        patch = PathPatch(MplPath(verts, [1, 4, 4, 4]), facecolor="none",
                          edgecolor="#728197", lw=1.5,
                          linestyle="--" if n.kind in {"drop", "detail"} else "-")
        ax.add_patch(patch)

    for n in nodes:
        x, y = pos[n.node_id]
        fill, text_color = COLORS[n.kind]
        ax.add_patch(FancyBboxPatch(
            (x-box_w/2, y-box_h/2), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=13", facecolor=fill,
            edgecolor=text_color, linewidth=1.2,
        ))
        for i, line in enumerate(node_lines(n)):
            ax.text(x, y - 30 + i*23, line, ha="center", va="center",
                    fontsize=9.5 if i == 0 else 8.5,
                    fontweight="bold" if i == 0 else "normal", color=text_color)
    ax.text(40, height-30,
            "Solid path = retained cohort; dashed branches = exclusions or reason breakdown. "
            "Counts are recalculated from the current CSV files.",
            fontsize=7.5, color="#566474")
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                        help="directory containing overview/status CSV files")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes = build_nodes(args.root)
    stem = args.out_dir / "database_filtering_tree"
    write_svg(nodes, stem.with_suffix(".svg"))
    png_written = write_png(nodes, stem.with_suffix(".png"))

    final = next(n for n in nodes if n.node_id == "after_feature")
    print(f"Wrote {stem.with_suffix('.svg')}")
    if png_written:
        print(f"Wrote {stem.with_suffix('.png')}")
    else:
        print("matplotlib unavailable; skipped PNG")
    print("Final cohort: "
          f"{final.counts.papers:,} papers, {final.counts.samples:,} samples, "
          f"{final.counts.features:,} features")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
