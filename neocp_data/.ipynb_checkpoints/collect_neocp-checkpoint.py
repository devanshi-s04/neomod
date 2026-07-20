#!/usr/bin/env python3
"""Collect MPC NEOCP list snapshots and ephemerides.

This script is intentionally self-contained and stdlib-only so cron can run it
without needing the notebook/conda environment. It stores raw MPC responses plus
append-only CSV tables for later analysis.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import html
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


LIST_URL = "https://minorplanetcenter.net/iau/NEO/neocp.txt"
HTML_URL = "https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html"
EPHEM_URL = "https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi"
USER_AGENT = "NEOMOD-NEOCP-collector/0.1 (local research archive)"

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "raw"
SNAPSHOT_DIR = BASE_DIR / "snapshots"
EPHEM_DIR = BASE_DIR / "ephemerides"
TABLE_DIR = BASE_DIR / "tables"
LOG_DIR = BASE_DIR / "logs"

OBJECT_HISTORY = TABLE_DIR / "neocp_objects_history.csv"
EPHEM_HISTORY = TABLE_DIR / "neocp_ephemerides_history.csv"
RUN_HISTORY = TABLE_DIR / "neocp_runs.csv"
LOCK_FILE = BASE_DIR / ".collect_neocp.lock"


NEOCP_LINE_RE = re.compile(
    r"^(?P<designation>\S+)\s+"
    r"(?P<score>\d+)\s+"
    r"(?P<disc_year>\d{4})\s+(?P<disc_month>\d{2})\s+(?P<disc_day>\d+(?:\.\d+)?)\s+"
    r"(?P<ra_hours>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<dec_deg>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<v_mag>\d+(?:\.\d+)?)\s+"
    r"(?P<status>Added|Updated)\s+"
    r"(?P<updated_month>[A-Za-z]+)\s+"
    r"(?P<updated_day>\d+(?:\.\d+)?)\s+UT\s+"
    r"(?:(?P<note>[A-Z])\s+)?"
    r"(?P<nobs>\d+)\s+"
    r"(?P<arc_days>\d+(?:\.\d+)?)\s+"
    r"(?P<h_mag>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<not_seen_days>\d+(?:\.\d+)?)\s*$"
)

EPHEM_OBJECT_RE = re.compile(
    r"<p><b>(?P<object>[^<]+)</b>.*?<pre>(?P<pre>.*?)</pre>",
    re.IGNORECASE | re.DOTALL,
)
TAG_RE = re.compile(r"<[^>]+>")
EPHEM_LINE_RE = re.compile(
    r"^\s*(?P<year>\d{4})\s+(?P<month>\d{2})\s+(?P<day>\d{2})\s+"
    r"(?P<hour>\d{2})\s+"
    r"(?P<ra_h>\d{2})\s+(?P<ra_m>\d{2})\s+(?P<ra_s>\d{2}(?:\.\d+)?)\s+"
    r"(?P<dec_sign>[+-])(?P<dec_d>\d{2})\s+(?P<dec_m>\d{2})\s+(?P<dec_s>\d{2})\s+"
    r"(?P<elong_deg>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<v_mag>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<motion_arcsec_per_min>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<pa_deg>[+-]?\d+(?:\.\d+)?)\s*"
    r"(?P<uncertainty>.*)$"
)


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def timestamp_id(now: datetime) -> str:
    return now.strftime("%Y%m%dT%H%M%SZ")


def ensure_dirs() -> None:
    for path in [RAW_DIR, SNAPSHOT_DIR, EPHEM_DIR, TABLE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def fetch_text(url: str, *, data: list[tuple[str, str]] | None = None) -> str:
    body = urlencode(data).encode("ascii") if data is not None else None
    req = Request(
        url,
        data=body,
        headers={
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST" if data is not None else "GET",
    )
    with urlopen(req, timeout=45) as response:
        raw = response.read()
    return raw.decode("iso-8859-1", errors="replace")


def append_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    if not rows:
        return
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    return float(value)


def to_int(value: str | None) -> int | None:
    if value in (None, ""):
        return None
    return int(value)


def parse_neocp_txt(text: str, *, snapshot_id: str, fetched_at: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        match = NEOCP_LINE_RE.match(line)
        if not match:
            print(f"WARN: could not parse NEOCP line: {line}", file=sys.stderr)
            continue
        item = match.groupdict()
        discovery_date = (
            f"{item['disc_year']}-{item['disc_month']}-{float(item['disc_day']):05.2f}"
        )
        rows.append(
            {
                "snapshot_id": snapshot_id,
                "fetched_at_utc": fetched_at,
                "designation": item["designation"],
                "score": to_int(item["score"]),
                "score_0_1": to_float(item["score"]) / 100.0,
                "discovery_date_mpc": discovery_date,
                "ra_hours": to_float(item["ra_hours"]),
                "ra_deg": to_float(item["ra_hours"]) * 15.0,
                "dec_deg": to_float(item["dec_deg"]),
                "v_mag": to_float(item["v_mag"]),
                "status": item["status"],
                "updated_mpc": f"{item['updated_month']} {item['updated_day']} UT",
                "note": item.get("note") or "",
                "nobs": to_int(item["nobs"]),
                "arc_days": to_float(item["arc_days"]),
                "h_mag": to_float(item["h_mag"]),
                "not_seen_days": to_float(item["not_seen_days"]),
                "source_url": LIST_URL,
            }
        )
    return rows


def parse_active_designations(html_text: str) -> list[str]:
    values = re.findall(
        r'<input\s+type="checkbox"\s+name="obj"\s+VALUE="([^"]+)"',
        html_text,
        flags=re.IGNORECASE,
    )
    seen: set[str] = set()
    active: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            active.append(value)
    return active


def hms_to_deg(h: str, m: str, s: str) -> float:
    return (float(h) + float(m) / 60.0 + float(s) / 3600.0) * 15.0


def dms_to_deg(sign: str, d: str, m: str, s: str) -> float:
    value = float(d) + float(m) / 60.0 + float(s) / 3600.0
    return -value if sign == "-" else value


def strip_tags(value: str) -> str:
    return html.unescape(TAG_RE.sub("", value)).strip()


def parse_ephemeris_html(
    html_text: str,
    *,
    snapshot_id: str,
    fetched_at: str,
    batch_id: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for object_match in EPHEM_OBJECT_RE.finditer(html_text):
        designation = strip_tags(object_match.group("object"))
        pre_text = strip_tags(object_match.group("pre"))
        for line in pre_text.splitlines():
            match = EPHEM_LINE_RE.match(line)
            if not match:
                continue
            item = match.groupdict()
            ra_deg = hms_to_deg(item["ra_h"], item["ra_m"], item["ra_s"])
            dec_deg = dms_to_deg(
                item["dec_sign"], item["dec_d"], item["dec_m"], item["dec_s"]
            )
            rows.append(
                {
                    "snapshot_id": snapshot_id,
                    "fetched_at_utc": fetched_at,
                    "batch_id": batch_id,
                    "designation": designation,
                    "ephem_time_utc": (
                        f"{item['year']}-{item['month']}-{item['day']}T"
                        f"{item['hour']}:00:00Z"
                    ),
                    "ra_hms": f"{item['ra_h']} {item['ra_m']} {item['ra_s']}",
                    "ra_deg": ra_deg,
                    "dec_dms": (
                        f"{item['dec_sign']}{item['dec_d']} "
                        f"{item['dec_m']} {item['dec_s']}"
                    ),
                    "dec_deg": dec_deg,
                    "elong_deg": to_float(item["elong_deg"]),
                    "v_mag": to_float(item["v_mag"]),
                    "motion_arcsec_per_min": to_float(item["motion_arcsec_per_min"]),
                    "pa_deg": to_float(item["pa_deg"]),
                    "uncertainty": item["uncertainty"].strip(),
                    "source_url": EPHEM_URL,
                }
            )
    return rows


def chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for idx in range(0, len(items), size):
        yield items[idx : idx + size]


def ephemeris_form_data(objects: list[str], args: argparse.Namespace) -> list[tuple[str, str]]:
    start_hours = f"{args.start_hours:g}"
    min_altitude = f"{args.min_altitude:g}"
    data: list[tuple[str, str]] = [("W", "j")]
    data.extend(("obj", obj) for obj in objects)
    data.extend(
        [
            ("Parallax", args.parallax),
            ("obscode", args.obscode),
            ("int", args.interval),
            ("start", start_hours),
            ("raty", "a"),
            ("mot", "m"),
            ("dmot", "p"),
            ("out", "f"),
            ("sun", args.sun),
            ("oalt", min_altitude),
        ]
    )
    return data


def run(args: argparse.Namespace) -> int:
    ensure_dirs()
    started = utc_now()
    snapshot_id = timestamp_id(started)
    fetched_at = started.isoformat().replace("+00:00", "Z")
    run_status = "ok"
    error_message = ""
    object_count = 0
    active_count = 0
    ephem_count = 0

    with LOCK_FILE.open("w") as lock_handle:
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("Another NEOCP collection run is already active; exiting.")
            return 0

        try:
            list_text = fetch_text(LIST_URL)
            html_text = fetch_text(HTML_URL)
            (RAW_DIR / f"neocp_{snapshot_id}.txt").write_text(list_text, encoding="utf-8")
            (RAW_DIR / f"toconfirm_{snapshot_id}.html").write_text(
                html_text, encoding="utf-8"
            )

            object_rows = parse_neocp_txt(
                list_text, snapshot_id=snapshot_id, fetched_at=fetched_at
            )
            object_count = len(object_rows)
            object_fields = [
                "snapshot_id",
                "fetched_at_utc",
                "designation",
                "score",
                "score_0_1",
                "discovery_date_mpc",
                "ra_hours",
                "ra_deg",
                "dec_deg",
                "v_mag",
                "status",
                "updated_mpc",
                "note",
                "nobs",
                "arc_days",
                "h_mag",
                "not_seen_days",
                "source_url",
            ]
            write_rows(
                SNAPSHOT_DIR / f"neocp_objects_{snapshot_id}.csv",
                object_rows,
                object_fields,
            )
            append_rows(OBJECT_HISTORY, object_rows, object_fields)

            active_objects = parse_active_designations(html_text)
            active_count = len(active_objects)
            if args.max_objects:
                active_objects = active_objects[: args.max_objects]

            ephem_fields = [
                "snapshot_id",
                "fetched_at_utc",
                "batch_id",
                "designation",
                "ephem_time_utc",
                "ra_hms",
                "ra_deg",
                "dec_dms",
                "dec_deg",
                "elong_deg",
                "v_mag",
                "motion_arcsec_per_min",
                "pa_deg",
                "uncertainty",
                "source_url",
            ]
            all_ephem_rows: list[dict[str, object]] = []
            if not args.no_ephemerides and active_objects:
                for batch_index, batch in enumerate(chunks(active_objects, args.batch_size), 1):
                    batch_id = f"{snapshot_id}_batch{batch_index:02d}"
                    eph_html = fetch_text(EPHEM_URL, data=ephemeris_form_data(batch, args))
                    (EPHEM_DIR / f"neocp_ephemerides_{batch_id}.html").write_text(
                        eph_html, encoding="utf-8"
                    )
                    rows = parse_ephemeris_html(
                        eph_html,
                        snapshot_id=snapshot_id,
                        fetched_at=fetched_at,
                        batch_id=batch_id,
                    )
                    all_ephem_rows.extend(rows)
                    if args.request_pause > 0 and batch_index * args.batch_size < len(active_objects):
                        time.sleep(args.request_pause)

            ephem_count = len(all_ephem_rows)
            write_rows(
                EPHEM_DIR / f"neocp_ephemerides_{snapshot_id}.csv",
                all_ephem_rows,
                ephem_fields,
            )
            append_rows(EPHEM_HISTORY, all_ephem_rows, ephem_fields)

        except (HTTPError, URLError, TimeoutError, OSError, ValueError) as exc:
            run_status = "error"
            error_message = f"{type(exc).__name__}: {exc}"
            print(error_message, file=sys.stderr)

        finished = utc_now()
        run_row = {
            "snapshot_id": snapshot_id,
            "started_at_utc": fetched_at,
            "finished_at_utc": finished.isoformat().replace("+00:00", "Z"),
            "status": run_status,
            "object_count": object_count,
            "active_ephemeris_object_count": active_count,
            "ephemeris_row_count": ephem_count,
            "error_message": error_message,
        }
        append_rows(
            RUN_HISTORY,
            [run_row],
            [
                "snapshot_id",
                "started_at_utc",
                "finished_at_utc",
                "status",
                "object_count",
                "active_ephemeris_object_count",
                "ephemeris_row_count",
                "error_message",
            ],
        )
        print(
            f"{snapshot_id}: status={run_status}, objects={object_count}, "
            f"active={active_count}, ephemeris_rows={ephem_count}"
        )
        return 0 if run_status == "ok" else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-ephemerides", action="store_true")
    parser.add_argument("--max-objects", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--request-pause", type=float, default=1.0)
    parser.add_argument(
        "--parallax",
        choices=["0", "1"],
        default="0",
        help="0=geocentric, 1=observatory code. Default is geocentric.",
    )
    parser.add_argument("--obscode", default="500")
    parser.add_argument(
        "--interval",
        choices=["0", "1", "2", "3"],
        default="0",
        help="MPC interval option: 0=1 hour, 1=30 min, 2=10 min, 3=1 min.",
    )
    parser.add_argument("--start-hours", type=float, default=0.0)
    parser.add_argument("--sun", choices=["x", "s", "c", "n", "a"], default="x")
    parser.add_argument("--min-altitude", type=float, default=0.0)
    return parser


def main() -> int:
    return run(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
