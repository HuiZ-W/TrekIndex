#!/usr/bin/env bash
# Run every experiment of the paper and regenerate its tables, plot data and text macros.
# This is days of single-threaded computation on the six graphs; the runs of scripts/paper_jobs.py can also
# be started one by one (see README.md). A run whose directory holds DONE is skipped, and run_jobs.py skips
# finished jobs, so the script can be restarted.
#   usage: scripts/run_paper_experiments.sh DATA RESULTS [CPUS]
#   DATA     the six graphs (<key>.bin) and the replay inputs (see scripts/paper_jobs.py)
#   RESULTS  output directory, one sub-directory per run
#   CPUS     comma-separated CPU ids, one job per CPU at a time (default 1; pinning needs Linux taskset)
set -euo pipefail
cd "$(dirname "$0")/.."
DATA=$1; RES=$2; CPUS=${3:-1}
make all
python3 scripts/paper_jobs.py jobs --data "$DATA" --results "$RES"
RUNS=$(python3 -c 'import sys; sys.path.insert(0, "scripts"); import paper_jobs; print(" ".join(paper_jobs.RUNS))')
for run in $RUNS; do
  [ -f "$RES/$run/DONE" ] && continue
  python3 scripts/run_jobs.py --jobs "$RES/$run/jobs.jsonl" --cores "$CPUS" --cwd . --out "$RES/$run"
  failed=$(grep -L '"exit": 0' "$RES/$run"/*.status || true)
  if [ -n "$failed" ]; then echo "$run: failed jobs: $failed" >&2; exit 1; fi
  touch "$RES/$run/DONE"
done
python3 scripts/paper_jobs.py compare "$RES/classic-v1" "$RES/e2-merged"
python3 scripts/paper_jobs.py compare "$RES/ablation-b" "$RES/e2-merged"
cp "$DATA/causenet-mediation/causenet-mediation.stats.json" "$RES/"
python3 scripts/make_paper_artifacts.py --results "$RES" --out "$RES/../gen"
