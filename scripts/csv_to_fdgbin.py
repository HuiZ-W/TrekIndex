#!/usr/bin/env python3
"""Convert a nodes/edges CSV pair (d-collision repository format: nodes have a
`name` column, edges have `source` and `target`) to FDGBIN01, keeping node
order. Duplicate edges are dropped; edges to unknown nodes, self-loops and
cycles are errors. NUL characters inside names (present in CNFDAG) are
replaced identically in both files. Prints a JSON summary."""
import csv, json, struct, sys

def main(nodes_csv, edges_csv, out):
    csv.field_size_limit(sys.maxsize)
    ids = {}
    with open(nodes_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(line.replace("\0", "<NUL>") for line in f):
            if row["name"] in ids: raise SystemExit(f"duplicate node {row['name']!r}")
            ids[row["name"]] = len(ids)
    n = len(ids); edges = set(); dup = 0
    with open(edges_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(line.replace("\0", "<NUL>") for line in f):
            u, v = ids.get(row["source"]), ids.get(row["target"])
            if u is None or v is None: raise SystemExit(f"edge to unknown node: {row}")
            if u == v: raise SystemExit(f"self-loop: {row}")
            if (u, v) in edges: dup += 1
            edges.add((u, v))
    edges = sorted(edges)
    indeg = [0] * n; ch = [[] for _ in range(n)]
    for u, v in edges: indeg[v] += 1; ch[u].append(v)
    stack = [v for v in range(n) if indeg[v] == 0]; seen = 0
    while stack:
        u = stack.pop(); seen += 1
        for v in ch[u]:
            indeg[v] -= 1
            if indeg[v] == 0: stack.append(v)
    if seen != n: raise SystemExit("graph has a cycle")
    with open(out, "wb") as f:
        f.write(b"FDGBIN01"); f.write(struct.pack("<QQ", n, len(edges)))
        for u, v in edges: f.write(struct.pack("<QQ", u, v))
    print(json.dumps({"nodes": n, "edges": len(edges), "duplicate_edges_dropped": dup, "acyclic": True}))

if __name__ == "__main__":
    main(*sys.argv[1:4])
