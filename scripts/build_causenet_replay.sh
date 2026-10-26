#!/usr/bin/env bash
# Inputs of the literature replay (Exp-8): the CauseNet DAG by evidence-ordered cycle breaking
# (causenet-evidence-dag.bin, 12,186,022 vertices, 10,291,150 edges) and the back-door checks of the
# mediation hypotheses of Toonsi et al. mapped onto its concepts (causenet-mediation/causenet-observed.requests,
# causenet-all.requests, causenet-mediation.stats.json).
#   CAUSENET_BZ2  causenet-full.jsonl.bz2 of CauseNet-Full (Heindorf et al., CIKM 2020; CC BY 4.0)
#   STUDY_DIR     checkout of github.com/bio-ontology-research-group/Mediation-Analysis-using-Causal-Knowledge-Graph
#                 (commit 8175041; only results/*_significant_results.tsv are read)
# usage: scripts/build_causenet_replay.sh CAUSENET_BZ2 STUDY_DIR OUT_DIR     (run `make` first)
set -euo pipefail
cd "$(dirname "$0")/.."
BZ2=$1; STUDY=$2; OUT=$3
mkdir -p "$OUT"
python3 scripts/extract_causenet.py "$BZ2" "$OUT/causenet-full-relations.tsv"
./build/causenet_dag "$OUT/causenet-full-relations.tsv" "$OUT/causenet-evidence-dag"
python3 scripts/build_causenet_mediation.py --names "$OUT/causenet-evidence-dag.names.tsv" \
  --graph "$OUT/causenet-evidence-dag.bin" --study "$STUDY" --output "$OUT/causenet-mediation"
