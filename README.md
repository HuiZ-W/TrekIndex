# TrekIndex

TrekIndex is an index for d-separation queries and the d-separation-based criteria built on them (back-door,
the three rules of do-calculus, front-door, sequential back-door, mediation, and selection bias) on large
causal DAGs. It reduces every query to forbidden-set reachability tests on a *trek DAG* with four states per
vertex, contracts the tree-like periphery of the graph into shortcut arcs with witnesses, and keeps pruned
2-hop labels on the remaining core; a query explores a small query graph H_Z whose edge tests use per-query
hub bitmasks and label-pruned searches.

This repository contains the code of the paper *TrekIndex: A Single Index for d-Separation-Based Queries on
Large Causal DAGs*. The online appendix (proofs and additional results) is
[`paper/TrekIndex-appendix.pdf`](paper/TrekIndex-appendix.pdf).

## Repository layout

| Path | Content |
|---|---|
| `src/coretree.hpp` | TrekIndex: witness contraction of the periphery and the index on the core (`CoreTreeZSkeleton`) |
| `src/zskeleton.hpp` | trek-DAG states, 2-hop labels (`TrekLabels`), query-graph exploration (`ZSkeleton`), the index-free variants of Exp-4 (`TrekSearch`) |
| `src/causal.hpp` | the criteria on top of an index, with graph surgery at query time |
| `src/polytree.hpp`, `src/universal.hpp`, `src/paper_algorithms.hpp` | the d-separation baselines (see [Baselines](#baselines)) and shared graph code (`src/core.hpp`) |
| `src/surgery_baselines.hpp`, `src/backdoor_baselines.hpp`, `src/reach_labels.hpp` | the baselines for endpoint sets and graph surgery, used for the criteria |
| `src/zskeleton.cpp` | benchmark driver `build/zskeleton`: d-separation, criteria and application workloads for every method |
| `src/paper_baselines.cpp` | standalone RP-separation driver `build/paper-baselines-native` (RP column of Exp-1) |
| `tools/` | `causenet_dag.cpp` (relation list to DAG by evidence-ordered cycle breaking), `markov_workload.cpp` (local-Markov requests), `dump_running_example.cpp` (the index of the paper's running example) |
| `tests/` | invariant tests against independent oracles, and a brute-force check of the reduction |
| `scripts/` | dataset preparation, experiment job lists and runner, and the generator of the paper's tables and figures |
| `paper/TrekIndex-appendix.pdf` | the online appendix (proofs and additional results) |

Scripts (Python 3 and bash):

| Script | Purpose |
|---|---|
| `gen_synthetic_dags.py` | Barabási–Albert DAGs (BA 500k, BA 2M) |
| `build_ikraph_dag.sh` (with `extract_ikraph_causal.sh`, `sort_relations.py`) | the iKraph DAG from the iKraph release |
| `csv_to_fdgbin.py` | node and edge CSV files to FDGBIN01 (CNFDAG) |
| `build_causenet_replay.sh` (with `extract_causenet.py`, `build_causenet_mediation.py`) | the inputs of the literature replay |
| `paper_jobs.py`, `run_jobs.py`, `run_paper_experiments.sh` | job lists of the paper's runs and their execution |
| `replay_paper_baselines.py` | RP-separation on saved workloads through `build/paper-baselines-native` |
| `sample_workload.py` | uniform sample of a workload (PC tests) |
| `make_paper_artifacts.py` | tables, plot data and text macros of the paper from the results |
| `example.sh` | the small end-to-end example of `make example` |

In the code, the query graph H_Z is called the *Z-skeleton* and witness contraction the *core–tree
decomposition*; the method names of the driver are listed under [Methods](#methods).

## Build

Requirements: a C++17 compiler (GCC or Clang), GNU make, bash, and Python 3.9 or later (standard library only).
Each binary is compiled from a single source file with the headers of `src/`; there are no other dependencies.

    make            # build/zskeleton, build/paper-baselines-native and the tools of tools/
    make test       # invariant tests and the brute-force check of the reduction (about one minute)
    make example    # build TrekIndex on a small synthetic DAG and check its answers against every baseline

The default flags are `-std=c++17 -O3 -DNDEBUG`; set `CXX` to choose the compiler (`make CXX=g++-13`).
The paper's experiments used GCC 13.4 on Linux (x86-64), single-threaded. `make test-sanitize` runs the
invariant tests under AddressSanitizer and UndefinedBehaviorSanitizer.

## Input format

Graphs are binary edge lists in the format FDGBIN01: the 8 ASCII bytes `FDGBIN01`, the number of vertices `n`
and of edges `m` as little-endian unsigned 64-bit integers, then `m` pairs `(parent, child)` of little-endian
unsigned 64-bit vertex IDs in `0..n-1` (file size `24 + 16 m`). The graph must be acyclic.
`scripts/csv_to_fdgbin.py NODES.csv EDGES.csv OUT.bin` converts a node list (column `name`) and an edge list
(columns `source`, `target`), and `build/causenet_dag RELATIONS.tsv OUT_PREFIX` turns a list of weighted
relations (`cause<TAB>effect<TAB>support<TAB>...`) into a DAG by evidence-ordered cycle breaking.

A workload file (`--workload-in`, `--workload-out`) holds a count on its first line, then one d-separation
request `k z_1 ... z_k x y` per line.

## Usage

    ./build/zskeleton bench --graph G.bin --queries 1000 --zsizes 1,4,16 --zmodes random,ancestor \
        --pairmode walk --seed 1 --warmup 1 --warmup-disjoint 1 --methods bidirectional,ctskel

builds every listed method once, generates the requests of each configuration (or reads them with
`--workload-in F1,F2,...`), runs one untimed warm-up pass (`--warmup-disjoint 1`: on different requests),
times every method on the same requests, and prints one JSON line per configuration with, for each method,
the mean, median, p90, p99 and maximum time per query, the build time and index statistics, and the number of
answers that differ from those of the first method (`mismatches_vs_first`, 0 in every run of the paper). Subcommands:

| Subcommand | Workload |
|---|---|
| `bench` | d-separation: `x` uniform, `y` at the end of a random walk of 1 to 6 steps (`--pairmode walk`), `Z` of `--zsizes` vertices drawn uniformly (`random`) or as ends of upward walks (`ancestor`) |
| `causal --kinds backdoor,rule1,rule2,rule3,frontdoor,seqbackdoor,mediation,selection` | the eight criteria (`--zmodes` also accepts `parents` for back-door and `children` for front-door) |
| `screen --hypotheses 100000` | covariate screening (Exp-8) |
| `replay --requests FILE` | replay of back-door checks (literature replay, Exp-8) |
| `pcgen` | the conditional-independence tests of the PC algorithm on ancestral sets (Exp-8) |
| `stats` | structural statistics of TrekIndex (shortcut arcs, label sizes; not timed) |

### Methods

| Paper | `--methods` name | Code |
|---|---|---|
| TrekIndex | `ctskel` | `CoreTreeZSkeleton` (`src/coretree.hpp`) |
| TrekIndex with / without surgery-first exploration (Exp-6) | `ctskel-sf` / `ctskel` with `--surgery-first 0` | (surgery-first is on by default) |
| 2-hop labels on the full trek DAG, without witness contraction (Exp-5) | `zskel` | `ZSkeleton` (`src/zskeleton.hpp`) |
| Trek traversal, trek search (Exp-4) | `trek-bidir`, `trek-edge` | `TrekSearch` (`src/zskeleton.hpp`) |
| Geiger | `geiger` | `GeigerReachable` (`src/polytree.hpp`) |
| Bidirectional Bayes-Ball (BB) | `bidirectional` | `BidirectionalCIfly` (`src/polytree.hpp`) |
| Ancestral Bayes-Ball (Anc-BB) | `bidirectional-anc` | `AncestralBidirectional` (`src/polytree.hpp`) |
| RP-separation (RP) | `rp` (or `build/paper-baselines-native`) | `RPReachable` (`src/paper_algorithms.hpp`) |
| CIfly | `cifly` | `StreamingCIfly` (`src/universal.hpp`) |
| d-collision | `dcollision` | `DCollision` (`src/polytree.hpp`) |

For the criteria, each baseline runs with endpoint sets and graph surgery (`src/surgery_baselines.hpp`);
`traversal` is BB with surgery and BFS side conditions, the reference of the older runs `round2` and `round3`.
The driver also accepts development variants that the paper does not report (`zskel-adaptive`, `darwiche`,
`trek-sweep`, `dcollision-full`, `ctskel-noclo`, and the `-r` variants with reachability labels).

### Baselines

All baselines are implemented in C++ in this repository; they stop as soon as the two ends of a test are
connected, reuse their workspace across queries, and skip the ignored edges.

- Geiger: the procedure of D. Geiger, T. Verma, J. Pearl, *Identifying Independence in Bayesian Networks*,
  Networks 20(5), 1990.
- Bayes-Ball: R. D. Shachter, *Bayes-Ball: The Rational Pastime*, UAI 1998, searched from both ends.
- Ancestral Bayes-Ball: Bayes-Ball restricted to the ancestral set, after S. L. Lauritzen, A. P. Dawid,
  B. N. Larsen, H.-G. Leimer, *Independence Properties of Directed Markov Fields*, Networks 20(5), 1990.
- RP-separation: C. J. Butz, A. E. dos Santos, J. S. Oliveira, *Relevant Path Separation: A Faster Method for
  Testing Independencies in Bayesian Networks*, PGM 2016 (Algorithm 1).
- CIfly: our C++ port of the rule-table approach of M. Wienöbst, S. Weichwald, L. Henckel, *Linear-Time
  Primitives for Algorithm Development in Graphical Causal Inference*, arXiv:2506.15758, 2025
  (original implementation: https://github.com/mwien/CIfly), with the Bayes-Ball rules for d-separation.
- d-collision: our C++ port of the Neo4j user-defined function of the d-collision graph method of
  M. Palmiotto, A. Bonifati, A. Mauri, *Computing d-Separation on Causal Property Graphs*, PACMMOD 4(6), 2026.

The ports are ours; credit for the methods belongs to their authors. The core labels of TrekIndex are 2-hop
reachability labels (E. Cohen, E. Halperin, H. Kaplan, U. Zwick, SIAM J. Comput. 32(5), 2003) built by pruned
search in the hub order of distribution labeling (R. Jin, G. Wang, PVLDB 6(14), 2013).

## Datasets

The six processed graphs of the paper (Table 3) and the inputs of the literature replay are available at
**https://cuhko365-my.sharepoint.com/:f:/g/personal/224045005_link_cuhk_edu_cn/IgD1p2-LA_WzTaW2TWB1A5tuAUmTz1mlOIynZJrTdjNtfNg?e=9Uxevx**
(zstd-compressed, with SHA-256 checksums and a README on their sources and licenses). After `zstd -d` the
folder is the `DATA` directory of the scripts below.

| File | Graph | \|V\| | \|E\| | Source | Preparation |
|---|---|---:|---:|---|---|
| `syn-ba-500k.bin` | BA 500k | 500,000 | 1,499,991 | Barabási–Albert, 3 edges per vertex | `scripts/gen_synthetic_dags.py OUT --sizes 500000,2000000 --families BA` |
| `aser-concept-dag.bin` | ASER | 1,310,917 | 8,657,788 | ASER concept graph (Zhang et al., WWW 2020) | processed file only (see the dataset README) |
| `ikraph-causal-dag.bin` | iKraph | 1,944,437 | 5,319,797 | iKraph, Zenodo record 14851275 (CC BY 4.0) | `scripts/build_ikraph_dag.sh iKraph_full.tar.gz OUT` |
| `syn-ba-2000k.bin` | BA 2M | 2,000,000 | 5,999,991 | Barabási–Albert, 3 edges per vertex | as BA 500k |
| `dcol-cnfdag.bin` | CNFDAG | 5,394,394 | 7,423,212 | CNFDAG of the d-collision work (from CauseNet, CC BY 4.0) | `scripts/csv_to_fdgbin.py nodes_CNFDAG.csv edges_CNFDAG.csv OUT.bin` |
| `causenet-full-dag.bin` | CauseNet | 12,186,310 | 10,840,656 | CauseNet-Full (Heindorf et al., CIKM 2020; CC BY 4.0) | processed file only (see the dataset README) |
| `causenet-evidence-dag.bin`, `causenet-mediation/` | literature replay (Exp-8) | 12,186,022 | 10,291,150 | CauseNet-Full and the results of Toonsi et al. | `scripts/build_causenet_replay.sh causenet-full.jsonl.bz2 STUDY_DIR OUT` |

The folder also holds the vertex names of iKraph and of the replay graph (`*.names.tsv`). The scripts
reproduce the distributed files byte for byte, except ASER and CauseNet, which were made
acyclic by an earlier data pipeline that is not part of this repository (the rule is described in the
dataset README).

## Running the experiments

`scripts/paper_jobs.py` writes the job list of every run of the paper (one directory per run, with the
commands of the paper), and `scripts/run_jobs.py` executes a list, one single-threaded job per CPU at a time
(pinned with `taskset` on Linux). Every job leaves `NAME.stdout` (JSON lines), `NAME.stderr` and `NAME.status`
in its run directory.

    make
    python3 scripts/paper_jobs.py jobs --data DATA --results results
    python3 scripts/run_jobs.py --jobs results/e2-merged/jobs.jsonl --cores 1 --cwd . --out results/e2-merged

`scripts/run_paper_experiments.sh DATA results [CPUS]` runs all runs in order, writes the files the generator
checks (`DONE`, `workloads-vs-e2.txt`), and regenerates the tables. The complete suite takes several days of
single-threaded computation. In the paper, every job ran single-threaded and pinned to its own core of a
server with two Intel Xeon Gold 6226R CPUs (2.9 GHz) and 512 GB of memory.

| Experiment (paper) | Runs |
|---|---|
| Table 3, Exp-3: index construction (mean of three builds) | `round2` (`round2-extra` for iKraph): `sweep-s{1,2,3}-G`; `index-stats` |
| Exp-1: d-separation | `e2-merged` (BB, Anc-BB, CIfly, TrekIndex; RP through `paper-baselines-native` on the same requests), `classic-v1` (Geiger, d-collision) |
| Exp-2: criteria | `criteria-all` (every baseline and TrekIndex), `classic-v1` (Geiger, d-collision), `round2`/`round2-extra` `causal-G` and `round3` (three seeds) |
| Exp-4: effect of the index | `ablation-b` |
| Exp-5: witness contraction | `round2` sweeps (construction), `decomp-order{,-r2,-r3}` (query time, both method orders, three repetitions) |
| Exp-6: surgery-first exploration | `surgery-order{,-r2,-r3}` |
| Exp-7: where edge tests are decided | `edge-tests`, checked against `ablation-final` |
| Exp-8: applications | `apps-all` (screening, replay), `pc-v2` (PC tests), `classic-v1` (Geiger, d-collision) |
| inputs and checks only | `markov-v1` (local-Markov requests), `exp5-markov`, `closure-ablation{,-rev}{,-r2,-r3}` |

When two methods share a process, the one timed second gains from warm caches; the runs that compare two
variants of TrekIndex (`decomp-order`, `surgery-order`, `closure-ablation`) therefore use both method orders
and three repetitions, and the generator reports the geometric mean over the orders of the median of the
repetitions.

## Regenerating the paper's tables and figures

    python3 scripts/make_paper_artifacts.py --results results --out gen

reads the run directories above (the full list of files is in the docstring of the script), checks every
claim of the text against the data, and writes `numbers.tex` (text macros), `tab-*.tex` (tables), `*.dat`
(plot data) and `paper-data.json` (every number with the files it came from). It needs every run listed in the
table above, and the replay statistics `causenet-mediation.stats.json` (option `--replay-stats`, default
`results/causenet-mediation.stats.json`). A missing or failed job stops it with a message naming the file.
The experiments of the paper ran on intermediate commits of the development repository, which the comments of
the script name; this repository contains the final version of that code.
