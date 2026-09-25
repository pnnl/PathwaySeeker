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


def kegg_rest(path: str, retries: int = 5, backoff: float = 2.0, cache: bool = True) -> Optional[str]:
    """GET ``https://rest.kegg.jp/<path>``.

    Returns the response text, ``""`` when KEGG has no such entry, or ``None`` when every
    attempt failed (rate limiting, server or network errors). Failures are appended to
    :data:`failures`.
    """
    cf = _cache_file(path)
    if cache and cf.exists():
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
                if cache:
                    _store(path, text)
                return text
            err = f"HTTP {r.status_code}"
        time.sleep(backoff * (2 ** attempt))
    failures.append((path, err))
    return None


def _store(path: str, text: str) -> None:
    cf = _cache_file(path)
    cf.parent.mkdir(parents=True, exist_ok=True)
    cf.write_text(text)


def _chunks(items: List[str], n: int):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def prefetch_entries(db: str, ids, batch: int = 10) -> None:
    """Fetch many ``get/<db>:<id>`` entries, ``batch`` per request, into the cache.

    Later ``kegg_rest(f"get/{db}:{id}")`` calls are then cache hits. Entries KEGG does not
    return are left uncached, so single requests still handle them.
    """
    todo = [i for i in dict.fromkeys(str(x) for x in ids if x) if not _cache_file(f"get/{db}:{i}").exists()]
    for chunk in _chunks(todo, batch):
        mark = len(failures)
        text = kegg_rest("get/" + "+".join(f"{db}:{i}" for i in chunk), cache=False, retries=2)
        del failures[mark:]  # single requests retry these IDs later
        if not text:
            continue
        wanted = set(chunk)
        for entry in text.split("\n///"):
            m = re.match(r"\s*ENTRY\s+(\S+)", entry)
            if m and m.group(1) in wanted:
                _store(f"get/{db}:{m.group(1)}", entry.strip("\n") + "\n///\n")


def prefetch_links(target: str, db: str, ids, batch: int = 10) -> None:
    """Fetch ``link/<target>/<db>:<id>`` for many IDs, ``batch`` per request, into the cache."""
    todo = [i for i in dict.fromkeys(str(x) for x in ids if x)
            if not _cache_file(f"link/{target}/{db}:{i}").exists()]
    for chunk in _chunks(todo, batch):
        mark = len(failures)
        text = kegg_rest(f"link/{target}/" + "+".join(f"{db}:{i}" for i in chunk), cache=False, retries=2)
        del failures[mark:]  # single requests retry these IDs later
        if text is None:
            continue
        lines = {i: [] for i in chunk}
        for line in text.strip().split("\n"):
            src = line.split("\t")[0].split(":", 1)[-1]
            if src in lines:
                lines[src].append(line)
        for i, ls in lines.items():
            _store(f"link/{target}/{db}:{i}", "\n".join(ls) + ("\n" if ls else ""))


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
