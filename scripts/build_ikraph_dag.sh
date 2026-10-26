#!/usr/bin/env bash
# Build the iKraph DAG of the paper (ikraph-causal-dag.bin, 1,944,437 vertices, 5,319,797 edges) from
# iKraph_full.tar.gz (iKraph, Zenodo record 14851275, CC BY 4.0):
#   https://zenodo.org/records/14851275  (file iKraph_full.tar.gz)
# Directed PubMed relations of four types (scripts/extract_ikraph_causal.sh), sorted
# (scripts/sort_relations.py), and made acyclic by evidence-ordered cycle breaking (tools/causenet_dag.cpp:
# inside every strongly connected component, remove the edges of least support until no cycle remains).
# usage: scripts/build_ikraph_dag.sh iKraph_full.tar.gz OUT_DIR     (run `make` first)
set -euo pipefail
cd "$(dirname "$0")/.."
TAR=$1; OUT=$2
mkdir -p "$OUT"
tar -xzf "$TAR" -C "$OUT" iKraph_full/PubMedList.json
bash scripts/extract_ikraph_causal.sh "$OUT/iKraph_full/PubMedList.json" "$OUT/ikraph-relations.tsv"
python3 scripts/sort_relations.py "$OUT/ikraph-relations.tsv" "$OUT/ikraph-relations.sorted.tsv"
./build/causenet_dag "$OUT/ikraph-relations.sorted.tsv" "$OUT/ikraph-causal-dag"
rm -f "$OUT/iKraph_full/PubMedList.json"; rmdir "$OUT/iKraph_full" 2>/dev/null || true
cat "$OUT/ikraph-causal-dag.stats.json"; echo
