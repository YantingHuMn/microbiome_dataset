
from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import io
import json
import re
import shutil
import sys
import threading
import time
import unicodedata
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

csv.field_size_limit(sys.maxsize)

# Reuse stage6's already-tested download/parse machinery -- no reimplementation.
sys.path.insert(0, str((Path(__file__).resolve().parent.parent / "find_papers")))
from stage6_verify_abundance import (          # noqa: E402
    list_deposit_files, stream_download, iter_tables, sniff_frame,
    to_long_rows, check_dependencies, retain_raw_file, append_raw_manifest, safe_pivot,
    detect_taxonomy_lookup_table, build_id_to_lineage_map, split_taxonomy, classify_domain,
)
from _netutil import GLOBAL_THROTTLE           # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "myproject"))
from find_data_sources import find_data_sources
from find_candidate_files import find_candidate_files

paper_id = "10.1371/journal.pone.0348120"

sources = find_data_sources(paper_id)

print(json.dumps(sources, ensure_ascii=False, indent=2) if sources
      else "can't find public data source")

out_dir = Path("/hickory/proj/didonglab/dataset/virus/Database/data/abundance_matrix_precise")
result = find_candidate_files(
    paper_id=paper_id,
    output_dir=out_dir,
    data_sources=sources,
)

print(json.dumps(result, ensure_ascii=False, indent=2))