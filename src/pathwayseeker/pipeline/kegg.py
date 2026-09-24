"""KEGG REST access shared by the pipeline: retries with backoff, a disk cache, and a
record of requests that failed so a build never drops data silently.

The cache (default ``~/.cache/pathwayseeker/kegg``, override with
``PATHWAYSEEKER_KEGG_CACHE``) also freezes KEGG responses: a rebuild from the same cache
gives the same graph even after KEGG changes.
"""

import os
import re
import time
from pathlib import Path
from typing import List, Optional, Tuple

import requests

BASE_URL = "https://rest.kegg.jp"
MIN_INTERVAL = 0.35  # KEGG asks for no more than about 3 requests per second

failures: List[Tuple[str, str]] = []
_session = requests.Session()
_last_request = [0.0]


def cache_dir() -> Path:
    return Path(os.environ.get("PATHWAYSEEKER_KEGG_CACHE",
                               Path.home() / ".cache" / "pathwayseeker" / "kegg"))


def _cache_file(path: str) -> Path:
    return cache_dir() / (re.sub(r"[^A-Za-z0-9._-]+", "_", path) + ".txt")


def _throttle():
    wait = MIN_INTERVAL - (time.monotonic() - _last_request[0])
    if wait > 0:
        time.sleep(wait)
    _last_request[0] = time.monotonic()


def kegg_rest(path: str, retries: int = 5, backoff: float = 2.0) -> Optional[str]:
    """GET ``https://rest.kegg.jp/<path>``.

    Returns the response text, ``""`` when KEGG has no such entry, or ``None`` when every
    attempt failed (rate limiting, server or network errors). Failures are appended to
    :data:`failures`.
    """
    cf = _cache_file(path)
    if cf.exists():
        return cf.read_text()
    err = ""
    for attempt in range(retries):
        _throttle()
        try:
            r = _session.get(f"{BASE_URL}/{path}", timeout=30)
        except requests.RequestException as e:
            err = f"{type(e).__name__}: {e}"
        else:
            if r.status_code in (200, 404):
                text = r.text if r.status_code == 200 else ""
                cf.parent.mkdir(parents=True, exist_ok=True)
                cf.write_text(text)
                return text
            err = f"HTTP {r.status_code}"
        time.sleep(backoff * (2 ** attempt))
    failures.append((path, err))
    return None


def report_failures(since: int, step: str) -> List[Tuple[str, str]]:
    """Print and return the failures recorded after index ``since``."""
    new = failures[since:]
    if new:
        print(f"  WARNING: {step}: {len(new)} KEGG request(s) failed after retries; their data is "
              f"missing. Rerun the build to retry them (successful responses are cached).")
        for path, err in new[:10]:
            print(f"    {path}: {err}")
    return new


def entry_field(text: str, field: str) -> List[str]:
    """Lines of a KEGG flat-file field (the field line plus its continuation lines)."""
    out, capture = [], False
    for line in (text or "").split("\n"):
        if line.startswith(field):
            out.append(line[len(field):].strip())
            capture = True
        elif capture and line.startswith(" "):
            out.append(line.strip())
        elif capture:
            break
    return out
