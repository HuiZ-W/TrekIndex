#!/usr/bin/env python3
"""Replay saved d-separation workloads through one paper_baselines method.

Each workload file is the `zskeleton --workload-out` format. With --limit N,
only the first N requests are replayed (for per-request index rebuilds).
Output: one JSON line per workload with timing and the raw answer string,
so answers can be checked against zskeleton's afterwards.
"""
import argparse, json, pathlib, subprocess, sys, time

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-baselines", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--workloads", required=True, help="comma-separated workload files")
    ap.add_argument("--limit", type=int, default=0, help="replay only the first N requests (0 = all)")
    ap.add_argument("--tmp", required=True, type=pathlib.Path, help="directory for truncated workload files")
    a = ap.parse_args()
    a.tmp.mkdir(parents=True, exist_ok=True)
    for path in map(pathlib.Path, a.workloads.split(",")):
        lines = path.read_text().splitlines()
        count = len(lines) - 1
        src = path
        if a.limit and a.limit < count:
            count = a.limit
            src = a.tmp / f"{path.name}.first{count}"
            src.write_text("\n".join([str(count)] + lines[1:1 + count]) + "\n")
        row = {"workload": path.name, "method": a.method, "queries": count}
        t = time.time()
        with open(src) as fh:
            p = subprocess.run([a.paper_baselines, "bench", a.method, "17", a.graph], stdin=fh,
                               capture_output=True, text=True)
        row["process_s"] = round(time.time() - t, 3)
        if p.returncode != 0:
            row["status"] = f"error {p.returncode}: {p.stderr[-400:]}"
        else:
            r = json.loads(p.stdout.strip().splitlines()[-1])
            row.update(build_ms=r["build_ms"], query_ms=r["query_ms"], mean_us=1000 * r["query_ms"] / max(count, 1),
                       p50_us=r.get("p50_us"), p99_us=r.get("p99_us"), answers=r["answers"])
        print(json.dumps(row), flush=True)
        print(f"{path.name} {a.method}: {row.get('mean_us', row.get('status'))}", file=sys.stderr, flush=True)

if __name__ == "__main__":
    main()
