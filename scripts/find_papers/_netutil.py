"""Shared per-host rate limiting for stage3/4/6.

Parallelising these stages only pays off if concurrent workers don't all
hammer the same API at once -- NCBI/EBI/figshare/Zenodo will rate-limit or
temporarily ban an IP that does that, which is slower than staying serial.
Throttle enforces one minimum interval PER HOSTNAME, shared across every
worker thread, so N threads hitting N different hosts run fully in
parallel, while N threads hitting the SAME host queue up to that host's
own pace -- concurrency hides latency, it doesn't remove the courtesy limit.
"""
from __future__ import annotations

import threading
import time
import urllib.parse

# Seconds between requests to the same host. Conservative by default;
# NCBI raises its own limit with an API key (see NCBI_API_KEY below).
DEFAULT_MIN_INTERVAL = 0.34          # ~3 req/s, safe for most APIs
HOST_MIN_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.34,  # 3/s anonymous, overridden below if keyed
    "www.ncbi.nlm.nih.gov": 0.34,
    # www.ebi.ac.uk serves BOTH the Europe PMC search API (heavy,
    # pageSize=1000/resultType=core payloads) and the per-article
    # fullTextXML/supplementaryFiles endpoints (light, single-article).
    # Measured empirically: 8 concurrent workers at 0.10s against the
    # SEARCH endpoint triggered repeated HTTP 503s from EBI (absorbed by
    # stage1's retry/backoff, but that's borrowed time, not a safe rate).
    # 0.5s is conservative enough to avoid that while still letting
    # --workers overlap latency across different queries/articles.
    "www.ebi.ac.uk": 0.5,
    "api.figshare.com": 0.10,
    "ndownloader.figshare.com": 0.05,
    "zenodo.org": 0.20,
    "datadryad.org": 0.20,
    "api.osf.io": 0.20,
    "api.github.com": 0.80,
    "raw.githubusercontent.com": 0.05,
}


class Throttle:
    """Thread-safe: call wait(url) immediately before every request."""

    def __init__(self, min_interval: dict[str, float] | None = None, ncbi_api_key: str = ""):
        self._last: dict[str, float] = {}
        self._lock = threading.Lock()
        self._intervals = dict(HOST_MIN_INTERVAL)
        if min_interval:
            self._intervals.update(min_interval)
        if ncbi_api_key:
            # 10 req/s with a key instead of 3 req/s anonymous
            for h in ("eutils.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"):
                self._intervals[h] = 0.11

    def wait(self, url: str) -> None:
        host = urllib.parse.urlparse(url).netloc
        gap = self._intervals.get(host, DEFAULT_MIN_INTERVAL)
        if gap <= 0:
            return
        with self._lock:
            now = time.monotonic()
            due = self._last.get(host, 0.0) + gap
            sleep_for = due - now
            if sleep_for > 0:
                time.sleep(sleep_for)
                now = time.monotonic()
            self._last[host] = now


# One process-wide instance is fine: stage scripts are single-process,
# multi-threaded (ThreadPoolExecutor), so a module-level Throttle is shared
# by every worker thread automatically.
GLOBAL_THROTTLE = Throttle()
