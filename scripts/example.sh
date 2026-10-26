#!/usr/bin/env bash
# Small end-to-end example (make example): build TrekIndex on a synthetic DAG with a core and a
# tree-like periphery, answer d-separation and criteria queries with TrekIndex and every baseline,
# and check that all methods give the same answers. Takes well under a minute.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=build/example
mkdir -p "$OUT"

# 1. Graph: a Barabasi-Albert core of 3,000 vertices (scripts/gen_synthetic_dags.py, m = 3) plus
#    12,000 peripheral vertices, each attached by one edge (random direction) to an earlier vertex.
python3 - "$OUT/example.bin" <<'PY'
import random, struct, sys
sys.path.insert(0, 'scripts')
from gen_synthetic_dags import barabasi_albert, orient, write
rng = random.Random('example')
core = 3000
edges = orient(core, barabasi_albert(core, 3, rng), rng)
n = core + 12000
for v in range(core, n):
    u = rng.randrange(v)
    edges.append((u, v) if rng.random() < 0.5 else (v, u))
# the periphery edges go both ways; orient everything by a random topological order of all vertices
rank = list(range(n)); rng.shuffle(rank)
edges = [(u, v) if rank[u] < rank[v] else (v, u) for u, v in edges]
write(sys.argv[1], n, edges)
print(f'example graph: n={n} m={len(edges)}')
PY

check() {   # every method of every JSON row agrees with the first method of the row
  python3 - "$1" <<'PY'
import json, sys
rows = [json.loads(l) for l in open(sys.argv[1]) if l.startswith("{")]
bad = [m["method"] for r in rows for m in r["methods"] if m["mismatches_vs_first"]]
names = [m["method"] for m in rows[0]["methods"]]
print("%s: %d configurations x %d methods (%s), mismatches: %d" % (sys.argv[1], len(rows), len(names), ", ".join(names), len(bad)))
for r in rows[:2]:
    t = " ".join("%s=%.1fus" % (m["method"], m["mean_us"]) for m in r["methods"])
    print("   %s %s k=%s: %s" % (r.get("kind", "d-sep"), r["zmode"], r.get("zsize", r.get("k")), t))
sys.exit(1 if bad else 0)
PY
}

# 2. d-separation: BB (bidirectional) is the reference; ctskel is TrekIndex.
./build/zskeleton bench --graph "$OUT/example.bin" --queries 300 --zsizes 1,4,16 --zmodes random,ancestor \
  --pairmode walk --seed 1 --warmup 1 --warmup-disjoint 1 --progress 0 --emit-answers 1 \
  --methods bidirectional,geiger,bidirectional-anc,rp,cifly,dcollision,trek-bidir,trek-edge,zskel,ctskel \
  --workload-out "$OUT/wl" > "$OUT/dsep.jsonl" 2> "$OUT/dsep.stderr"
check "$OUT/dsep.jsonl"

# 3. The eight criteria (back-door, rules 1-3, front-door, sequential back-door, mediation, selection).
./build/zskeleton causal --graph "$OUT/example.bin" --queries 200 --zsizes 1,4 --zmodes random,ancestor \
  --kinds backdoor,rule1,rule2,rule3,frontdoor,seqbackdoor,mediation,selection --pairmode walk --seed 1 \
  --warmup 1 --warmup-disjoint 1 --progress 0 \
  --methods bidirectional,geiger,bidirectional-anc,rp,cifly,dcollision,ctskel,ctskel-sf > "$OUT/causal.jsonl" 2> "$OUT/causal.stderr"
check "$OUT/causal.jsonl"

# 4. RP-separation through the standalone driver on the saved workloads (as for Exp-1); its answers
#    must equal those of BB in step 2.
python3 scripts/replay_paper_baselines.py --paper-baselines build/paper-baselines-native --method rp-boolean \
  --graph "$OUT/example.bin" --workloads "$(ls "$OUT"/wl.* | paste -sd, -)" --tmp "$OUT/tmp" \
  > "$OUT/rp.jsonl" 2> "$OUT/rp.stderr"
python3 - "$OUT/dsep.jsonl" "$OUT/rp.jsonl" <<'PY'
import json, sys
bench = {(r["zmode"], r["zsize"]): next(m["answers"] for m in r["methods"] if m["method"] == "bidirectional")
         for r in map(json.loads, open(sys.argv[1]))}
bad = 0
rows = [json.loads(l) for l in open(sys.argv[2])]
for r in rows:
    zmode, k = r["workload"].rsplit(".", 1)[1].rsplit("-", 1)
    bad += r["answers"] != bench[(zmode, int(k))]
print("%s: RP replay of %d saved workloads, workloads with answers differing from BB: %d" % (sys.argv[2], len(rows), bad))
sys.exit(1 if bad or not rows else 0)
PY
echo "example: all methods agree"
