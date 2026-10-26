#!/usr/bin/env python3
"""Synthetic DAG families of the d-collision experiments (Palmiotto et al., PACMMOD 4(6), 2026; script
Scripts/S_scripts/09_S_DAG_generation.py of their artifact), re-implemented without networkx
and written as FDGBIN01 (magic, u64 n, u64 m, m x (u64 parent, u64 child)).

  BA  Barabasi-Albert, m = 3, edges oriented by a random topological order
  ER  Erdos-Renyi, oriented by a random order, largest weakly connected
      component kept. The original fixes p = 0.0025 (expected degree 25 at
      n = 10,000); larger n keep that expected degree (p = 25 / (n - 1))
      instead of the fixed p, which would give 3.1e8 edges at n = 500,000.
  LF  layered feed-forward: int(sqrt(n)) + 1 layers of random sizes in
      [L/2, 3L/2]; each vertex draws 1 + Binomial(|prev| - 1, p) parents in
      the previous layer, p = 5 / (L/2 - 1); a vertex without children draws
      1 + Binomial(|next| - 1, p) children in the next layer.
  TR  uniformly random labelled rooted tree, edges directed away from the root.

Usage: gen_synthetic_dags.py OUT_DIR [--sizes 10000,100000,500000] [--families BA,ER,LF,TR]
The paper uses BA 500k and BA 2M: gen_synthetic_dags.py OUT_DIR --sizes 500000,2000000 --families BA
(random.Random("BA-500000") and random.Random("BA-2000000"); the output is byte-identical to the paper's files).
"""
import argparse, heapq, json, pathlib, random, struct, sys

def binomial(rng, n, p):
    """Successes in n Bernoulli(p) trials by geometric skipping (O(np + 1))."""
    if n <= 0 or p <= 0: return 0
    if p >= 1: return n
    import math
    log_q = math.log(1.0 - p); k = 0; i = -1
    while True:
        i += int(math.log(1.0 - rng.random()) / log_q) + 1
        if i >= n: return k
        k += 1

def orient(n, edges, rng):
    rank = list(range(n)); rng.shuffle(rank)
    return [(u, v) if rank[u] < rank[v] else (v, u) for u, v in edges]

def barabasi_albert(n, m, rng):
    # networkx: star graph on m + 1 vertices, then each new vertex attaches to
    # m distinct targets drawn from the degree-weighted "repeated nodes" list.
    edges = [(0, i) for i in range(1, m + 1)]
    repeated = [0] * m + list(range(1, m + 1))
    for v in range(m + 1, n):
        targets = set()
        while len(targets) < m: targets.add(rng.choice(repeated))
        for t in targets: edges.append((t, v))
        repeated.extend(targets); repeated.extend([v] * m)
    return edges

def erdos_renyi(n, degree, rng):
    m = round(degree * n / 2)
    keys = set()
    while len(keys) < m:
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v: keys.add(min(u, v) * n + max(u, v))
    return [(k // n, k % n) for k in keys]

def largest_wcc(n, edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for u, v in edges:
        a, b = find(u), find(v)
        if a != b: parent[a] = b
    size = {}
    for v in range(n): r = find(v); size[r] = size.get(r, 0) + 1
    best = max(size, key=size.get)
    ids = {}
    for v in range(n):
        if find(v) == best: ids[v] = len(ids)
    return len(ids), [(ids[u], ids[v]) for u, v in edges if u in ids]

def layered(n, rng):
    L = int(n ** 0.5) + 1; low, high = L // 2, 3 * L // 2
    sizes = [low] * L; remaining = n - L * low
    while remaining > 0:
        i = rng.randrange(L)
        if sizes[i] < high: sizes[i] += 1; remaining -= 1
    rng.shuffle(sizes)
    p = 5 / (low - 1)
    layers, nid = [], 0
    for s in sizes: layers.append(list(range(nid, nid + s))); nid += s
    edges = set()
    for U, V in zip(layers, layers[1:]):
        out = set()
        for v in V:
            for u in rng.sample(U, binomial(rng, len(U) - 1, p) + 1): edges.add((u, v)); out.add(u)
        for u in U:
            if u not in out:
                for v in rng.sample(V, binomial(rng, len(V) - 1, p) + 1): edges.add((u, v))
    return nid, sorted(edges)

def random_tree(n, rng):
    # Uniform labelled tree from a random Pruefer sequence, then a random root.
    seq = [rng.randrange(n) for _ in range(n - 2)]
    deg = [1] * n
    for x in seq: deg[x] += 1
    leaves = [v for v in range(n) if deg[v] == 1]; heapq.heapify(leaves)
    adj = [[] for _ in range(n)]
    for x in seq:
        leaf = heapq.heappop(leaves); adj[leaf].append(x); adj[x].append(leaf)
        deg[x] -= 1
        if deg[x] == 1: heapq.heappush(leaves, x)
    a, b = heapq.heappop(leaves), heapq.heappop(leaves); adj[a].append(b); adj[b].append(a)
    root = rng.randrange(n); seen = [False] * n; seen[root] = True; order = [root]; edges = []
    for u in order:
        for w in adj[u]:
            if not seen[w]: seen[w] = True; order.append(w); edges.append((u, w))
    return edges

def write(path, n, edges):
    with open(path, "wb") as f:
        f.write(b"FDGBIN01"); f.write(struct.pack("<QQ", n, len(edges)))
        f.write(struct.pack(f"<{2 * len(edges)}Q", *[x for e in edges for x in e]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=pathlib.Path)
    ap.add_argument("--sizes", default="10000,100000,500000")
    ap.add_argument("--families", default="BA,ER,LF,TR")
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    stats = []
    for size in map(int, a.sizes.split(",")):
        for fam in a.families.split(","):
            rng = random.Random(f"{fam}-{size}")
            if fam == "BA": n, edges = size, orient(size, barabasi_albert(size, 3, rng), rng)
            elif fam == "ER": n, edges = largest_wcc(size, orient(size, erdos_renyi(size, 0.0025 * (10000 - 1), rng), rng))
            elif fam == "LF": n, edges = layered(size, rng)
            else: n, edges = size, random_tree(size, rng)
            name = f"syn-{fam.lower()}-{size // 1000}k"
            write(a.out / f"{name}.bin", n, edges)
            stats.append({"name": name, "family": fam, "target_n": size, "n": n, "m": len(edges)})
            print(json.dumps(stats[-1]), file=sys.stderr, flush=True)
    p = a.out / "synthetic-stats.json"
    old = json.loads(p.read_text()) if p.exists() else []
    merged = {s["name"]: s for s in old + stats}
    p.write_text(json.dumps(list(merged.values()), indent=1) + "\n")

if __name__ == "__main__":
    main()
