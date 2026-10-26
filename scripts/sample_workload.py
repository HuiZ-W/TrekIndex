"""Uniform sample of k requests of a workload file (first line: count; one request per line), in file order.

usage: sample_workload.py SRC DST K SEED   (PC workload of Exp-8: K = 10000, SEED = 1; prints a JSON summary)
"""
import random, sys, json
src, dst, k, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
with open(src) as f:
    n = int(f.readline()); lines = f.readlines()
assert len(lines) == n, (len(lines), n)
idx = sorted(random.Random(seed).sample(range(n), min(k, n)))
with open(dst, 'w') as f:
    f.write(f'{len(idx)}\n'); f.writelines(lines[i] for i in idx)
print(json.dumps({'source': src, 'tests': n, 'sample': len(idx), 'seed': seed}))
