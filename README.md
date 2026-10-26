# TrekIndex

An index for d-separation and d-separation-based queries on large causal DAGs, with six traversal baselines.
The full version of the paper, with its appendix, is [paper/TrekIndex-full.pdf](paper/TrekIndex-full.pdf).

## Installation

Requires a C++17 compiler (GCC or Clang), GNU make and Python ≥ 3.9 (standard library only).

```bash
make            # build/zskeleton (index, baselines, workloads) and the tools
make test       # invariant tests against independent oracles
make example    # small end-to-end run: every method gives the same answers
```

## Data formats

### Datasets

The paper uses four real-world causal DAGs and two synthetic DAGs. Click a dataset name for its original
source (or generator). The processed graphs, compressed with zstd and with SHA-256 checksums, are on
[OneDrive](https://cuhko365-my.sharepoint.com/:f:/g/personal/224045005_link_cuhk_edu_cn/IgD1p2-LA_WzTaW2TWB1A5tuAUmTz1mlOIynZJrTdjNtfNg?e=9Uxevx).

| Dataset | Nodes | Edges |
|---|---:|---:|
| [BA 500k](scripts/gen_synthetic_dags.py) | 500,000 | 1,499,991 |
| [ASER](https://github.com/HKUST-KnowComp/ASER) | 1,310,917 | 8,657,788 |
| [iKraph](https://zenodo.org/records/14851275) | 1,944,437 | 5,319,797 |
| [BA 2M](scripts/gen_synthetic_dags.py) | 2,000,000 | 5,999,991 |
| [CNFDAG](https://doi.org/10.1145/3856322) | 5,394,394 | 7,423,212 |
| [CauseNet](https://causenet.org/) | 12,186,310 | 10,840,656 |

### Formats

- **Graph (`.bin`, FDGBIN01):** the eight bytes `FDGBIN01`, then `n`, `m` and `m` edges `(parent, child)`
  as little-endian unsigned 64-bit integers. Vertex IDs are zero-based; the graph must be acyclic.
  `scripts/csv_to_fdgbin.py NODES.csv EDGES.csv OUT.bin` converts node and edge CSV files.
- **Workload (`--workload-in`, `--workload-out`):** a count, then one d-separation request
  `k z_1 ... z_k x y` per line.
- **Output:** one JSON line per configuration with, for each method, the mean, p99 and maximum time per
  query, the build time and index size, and `mismatches_vs_first`, the number of answers that differ from
  those of the first method.

## Running

Methods: `ctskel` (TrekIndex), `geiger`, `bidirectional` (BB), `bidirectional-anc` (Anc-BB), `rp`, `cifly`,
`dcollision`; for the ablations, `zskel` (2-hop labels without witness contraction), `trek-bidir` and
`trek-edge` (no index), and `ctskel-sf` (surgery-first exploration).

```bash
# d-separation (Exp-1): random and ancestor Z of 1 to 32 vertices, 1,000 queries each.
./build/zskeleton bench --graph DATA/syn-ba-500k.bin --queries 1000 --zsizes 1,2,4,8,16,32 \
    --zmodes random,ancestor --pairmode walk --seed 1 --warmup 1 --warmup-disjoint 1 \
    --methods bidirectional,geiger,cifly,dcollision,ctskel

# The eight criteria (Exp-2), with sets of 1, 4 and 16 vertices.
./build/zskeleton causal --graph DATA/syn-ba-500k.bin \
    --kinds backdoor,rule1,rule2,rule3,frontdoor,seqbackdoor,mediation,selection \
    --zmodes random,ancestor --queries 1000 --zsizes 1,4,16 --pairmode walk --seed 1 \
    --methods bidirectional,ctskel

# Applications (Exp-8): covariate screening, literature replay, conditional-independence tests of PC.
./build/zskeleton screen --graph DATA/causenet-full-dag.bin --hypotheses 100000 --methods bidirectional,ctskel
./build/zskeleton replay --graph DATA/causenet-evidence-dag.bin \
    --requests DATA/causenet-mediation/causenet-observed.requests --methods bidirectional,ctskel
./build/zskeleton pcgen --graph DATA/syn-ba-500k.bin --seed 1 --regions 20 --min-size 800 \
    --region-size 1000 --workload-out pc.workload

# Every run of the paper, then its tables, figures and numbers.
bash scripts/run_paper_experiments.sh DATA results
python3 scripts/make_paper_artifacts.py --results results --out gen
```

`DATA` is the decompressed dataset folder (`zstd -d *.zst causenet-mediation/*.zst`). Each method builds its
index once per graph; the reported query times exclude construction. The paper's runs were single-threaded
(GCC 13.4, `-O3`, one job per core), and the whole suite takes several days. `scripts/paper_jobs.py` lists
the individual runs, and `make_paper_artifacts.py` checks every claim of the paper against the results.
