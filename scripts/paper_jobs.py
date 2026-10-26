#!/usr/bin/env python3
"""Job lists of the paper's experiments, one run directory per list, for scripts/run_jobs.py.

The commands and run names are those of the paper's runs; scripts/make_paper_artifacts.py reads the
run directories under these names (see its docstring for which file feeds which table or figure).

  python3 scripts/paper_jobs.py jobs --data DATA --results RESULTS [--runs RUN,...]
      writes RESULTS/<run>/jobs.jsonl for every run (default: all, in the order of RUNS below).
  python3 scripts/paper_jobs.py compare RESULTS/<run> RESULTS/e2-merged
      writes RESULTS/<run>/workloads-vs-e2.txt: whether the d-separation requests a run wrote
      (--workload-out) are byte-identical to those of e2-merged (checked by the generator for
      classic-v1 and ablation-b).

DATA holds the six graphs as <key>.bin (keys below) and, for the literature replay (Exp-8),
causenet-evidence-dag.bin and causenet-mediation/causenet-observed.requests (causenet-all.requests is
optional; that workload is not in the paper). scripts/run_paper_experiments.sh runs everything in order.
"""
import argparse, filecmp, json, pathlib, sys

GRAPHS = ['syn-ba-500k', 'aser-concept-dag', 'ikraph-causal-dag', 'syn-ba-2000k', 'dcol-cnfdag', 'causenet-full-dag']
EST = {'syn-ba-500k': 10, 'aser-concept-dag': 30, 'ikraph-causal-dag': 20, 'syn-ba-2000k': 40, 'dcol-cnfdag': 50,
       'causenet-full-dag': 60}   # relative job lengths; run_jobs.py starts the longest first
ZB = './build/zskeleton'
PB = './build/paper-baselines-native'
COMMON = '--warmup 1 --warmup-disjoint 1 --progress 0'
DSEP = f'--queries 1000 --zsizes 1,2,4,8,16,32 --zmodes random,ancestor --pairmode walk --seed 1 {COMMON}'
CRIT = '--queries 1000 --zsizes 1,4,16 --pairmode walk --seed 1'
ALL_BASELINES = 'bidirectional,bidirectional-anc,cifly,rp,ctskel'
CLASSIC = 'bidirectional,geiger,dcollision'
REPS = ('', '-r2', '-r3')
RUNS = (['markov-v1', 'round2', 'round2-extra', 'e2-merged', 'criteria-all', 'round3', 'pc-v2', 'classic-v1', 'apps-all',
         'exp5-markov', 'ablation-final', 'ablation-b', 'edge-tests']
        + [f'decomp-order{r}' for r in REPS] + [f'surgery-order{r}' for r in REPS]
        + [f'closure-ablation{v}{r}' for v in ('', '-rev') for r in REPS] + ['index-stats'])


def wait_for(status):
    """Shell prefix that waits until another job of the same run has finished with exit 0."""
    return (f'while ! grep -q \'"exit": 0\' {status} 2>/dev/null; do '
            f'if [ -f {status} ]; then echo "{status}: failed" >&2; exit 1; fi; sleep 10; done; ')


def markov_files(w, g):
    return f'$(ls {w}/{g}.random-* {w}/{g}.walk-* | paste -sd, -)'


def run_jobs(run, data, res):
    out = res / run
    G = lambda g: data / f'{g}.bin'
    jobs = []
    job = lambda name, est, cmd: jobs.append({'name': name, 'est_min': est, 'cmd': cmd})
    mw = res / 'markov-v1' / 'workloads'
    if run == 'markov-v1':   # local-Markov requests (16 buckets of 1,000 per graph), input of round2*/markov-G and exp5-markov
        for g in GRAPHS:
            job(f'mw-{g}', EST[g], f'mkdir -p {mw} && ./build/markov_workload {G(g)} {mw}/{g} 1000 1')
    elif run in ('round2', 'round2-extra'):   # Table 3 / Exp-3 (three builds), Exp-5 (plain 2-hop), Exp-2 (TrekIndex vs BB)
        for g in GRAPHS:
            if (g == 'ikraph-causal-dag') != (run == 'round2-extra'):
                continue
            for s in (1, 2, 3):
                job(f'sweep-s{s}-{g}', 3 * EST[g] + s, f'{ZB} bench --graph {G(g)} --queries 1000 --zsizes 1,2,4,8,16,32 '
                    f'--zmodes random,ancestor --pairmode walk --seed {s} --methods bidirectional,zskel,ctskel {COMMON}')
            job(f'causal-{g}', 2 * EST[g], f'{ZB} causal --graph {G(g)} {CRIT} --zmodes parents,random,ancestor '
                f'--kinds backdoor,rule1,rule2,rule3 --methods traversal,zskel,ctskel {COMMON}')
            job(f'markov-{g}', EST[g], f'{ZB} bench --graph {G(g)} --workload-in {markov_files(mw, g)} '
                '--methods bidirectional,ctskel --warmup 1 --progress 0')
    elif run == 'e2-merged':   # Exp-1: BB, Anc-BB, CIfly, TrekIndex; RP replays the same requests with the standalone driver
        for g in GRAPHS:
            job(f'bench-s1-{g}', 1000 + EST[g], f'mkdir -p {out}/workloads && {ZB} bench --graph {G(g)} {DSEP} --methods bidirectional,bidirectional-anc,cifly,ctskel '
                f'--workload-out {out}/workloads/{g}-s1')
            job(f'rp-s1-{g}', EST[g], wait_for(f'{out}/bench-s1-{g}.status') +
                f'python3 scripts/replay_paper_baselines.py --paper-baselines {PB} --method rp-boolean --graph {G(g)} '
                f'--workloads $(ls {out}/workloads/{g}-s1.random-* {out}/workloads/{g}-s1.ancestor-* | paste -sd, -) --tmp {out}/tmp')
    elif run == 'criteria-all':   # Exp-2: every baseline and TrekIndex
        for g in GRAPHS:
            for name, kinds, zmodes in (('backdoor', 'backdoor', 'parents,random,ancestor'), ('rules', 'rule1,rule2,rule3', 'random,ancestor'),
                                        ('frontdoor', 'frontdoor', 'children,random,ancestor'),
                                        ('further', 'seqbackdoor,mediation,selection', 'random,ancestor')):
                job(f'{name}-{g}', EST[g], f'{ZB} causal --graph {G(g)} --kinds {kinds} --zmodes {zmodes} {CRIT} '
                    f'--methods {ALL_BASELINES} {COMMON}')
    elif run == 'round3':   # Exp-2 (TrekIndex vs BB, three seeds); surgery-first on (ctskel-sf) and off (ctskel)
        for g in GRAPHS:
            for s in (1, 2, 3):
                job(f'rules-s{s}-{g}', EST[g], f'{ZB} causal --graph {G(g)} --zsizes 1,4,16 --zmodes random,ancestor --kinds rule2,rule3 '
                    f'--methods traversal,ctskel,ctskel-sf --surgery-first 0 --seed {s} --queries 1000 --pairmode walk {COMMON}')
                job(f'further-s{s}-{g}', EST[g], f'{ZB} causal --graph {G(g)} --zsizes 1,4,16 --zmodes children,random,ancestor '
                    f'--kinds frontdoor,seqbackdoor,mediation,selection --methods traversal,ctskel --seed {s} --queries 1000 '
                    f'--pairmode walk {COMMON}')
    elif run == 'pc-v2':   # Exp-8, PC tests: level 1 on 20 ancestral sets of 800-1,000 vertices, uniform sample of 10,000 tests
        for g in GRAPHS:
            job(f'pc-{g}', EST[g], f'mkdir -p {out}/workloads && {ZB} pcgen --graph {G(g)} --seed 1 --regions 20 --min-size 800 '
                f'--region-size 1000 --workload-out {out}/workloads/{g}.pc > {out}/pcgen-{g}.json && '
                f'python3 scripts/sample_workload.py {out}/workloads/{g}.pc {out}/workloads/{g}.sample.pc 10000 1 > {out}/sample-{g}.json && '
                f'{ZB} bench --graph {G(g)} --workload-in {out}/workloads/{g}.sample.pc --methods ctskel,bidirectional-anc,rp,bidirectional,cifly '
                '--warmup 1 --progress 0')
    elif run == 'classic-v1':   # Geiger and d-collision on every workload of the paper, BB as the answer reference
        for g in GRAPHS:
            job(f'dsep-{g}', 2 * EST[g], f'mkdir -p {out}/workloads && {ZB} bench --graph {G(g)} {DSEP} --methods {CLASSIC} --workload-out {out}/workloads/{g}-s1')
            for name, kinds, zmodes in (('backdoor', 'backdoor', 'parents,random,ancestor'), ('rules', 'rule1,rule2,rule3', 'random,ancestor'),
                                        ('frontdoor', 'frontdoor', 'children,random,ancestor'),
                                        ('further', 'seqbackdoor,mediation,selection', 'random,ancestor')):
                job(f'{name}-{g}', EST[g], f'{ZB} causal --graph {G(g)} --kinds {kinds} --zmodes {zmodes} {CRIT} {COMMON} --methods {CLASSIC}')
            job(f'pc-{g}', EST[g] // 2, f'{ZB} bench --graph {G(g)} --workload-in {res}/pc-v2/workloads/{g}.sample.pc '
                f'--methods {CLASSIC} --warmup 1 --progress 0')
        for g in ('causenet-full-dag', 'aser-concept-dag'):
            job(f'screen-{g}', 3 * EST[g], f'{ZB} screen --graph {G(g)} --hypotheses 100000 --methods {CLASSIC} --progress 0')
        job('replay-causenet-observed', 200, f'{ZB} replay --graph {data}/causenet-evidence-dag.bin '
            f'--requests {data}/causenet-mediation/causenet-observed.requests --methods {CLASSIC} --warmup 1 --progress 0')
    elif run == 'apps-all':   # Exp-8: screening and literature replay, every baseline in back-door mode
        M = 'traversal,ctskel,bidirectional,bidirectional-anc,cifly,rp'
        for g in ('causenet-full-dag', 'aser-concept-dag'):
            job(f'screen-{g}', EST[g], f'{ZB} screen --graph {G(g)} --hypotheses 100000 --methods {M} --progress 0')
        for w in ('observed', 'all'):
            req = data / 'causenet-mediation' / f'causenet-{w}.requests'
            if w == 'all' and not req.exists():
                continue   # optional: not in the paper
            job(f'med-causenet-{w}', 50, f'{ZB} replay --graph {data}/causenet-evidence-dag.bin --requests {req} --methods {M} '
                '--warmup 1 --progress 0')
    elif run == 'exp5-markov':   # local-Markov tests (not in the paper; read by the generator's checks)
        for g in GRAPHS:
            job(f'markov-{g}', 100 + EST[g], f'{ZB} bench --graph {G(g)} --workload-in {markov_files(mw, g)} '
                '--methods bidirectional,bidirectional-anc,cifly,ctskel --warmup 1 --progress 0')
            job(f'rp-{g}', EST[g], f'python3 scripts/replay_paper_baselines.py --paper-baselines {PB} --method rp-boolean '
                f'--graph {G(g)} --workloads {markov_files(mw, g)} --tmp {out}/tmp')
    elif run == 'ablation-final':   # reference counters for Exp-6/7 (default, without up-closures, plain 2-hop labels)
        for g in GRAPHS:
            for name, flags in (('default', '--methods ctskel'), ('noclosure', '--methods ctskel --closure 0'), ('fulllabels', '--methods zskel')):
                job(f'{name}-{g}', EST[g], f'{ZB} bench --graph {G(g)} {DSEP} {flags}')
    elif run == 'ablation-b':   # Exp-4: BB, TrekIndex, trek traversal (trek-bidir), trek search (trek-edge) in one process
        for g in GRAPHS:
            job(f'abl-{g}', EST[g], f'mkdir -p {out}/workloads && {ZB} bench --graph {G(g)} {DSEP} --methods bidirectional,ctskel,trek-bidir,trek-edge '
                f'--workload-out {out}/workloads/{g}-s1')
    elif run == 'edge-tests':   # Exp-7: how TrekIndex decides each edge test of H_Z (counters et_*)
        for g in GRAPHS:
            job(f'et-{g}', EST[g], f'{ZB} bench --graph {G(g)} {DSEP} --methods ctskel')
    elif run.startswith('decomp-order'):   # Exp-5: plain 2-hop labels (zskel) and TrekIndex, both orders
        for g in GRAPHS:
            job(f'dc-ab-{g}', EST[g], f'{ZB} bench --graph {G(g)} {DSEP} --methods zskel,ctskel')
            job(f'dc-ba-{g}', EST[g], f'{ZB} bench --graph {G(g)} {DSEP} --methods ctskel,zskel')
    elif run.startswith('surgery-order'):   # Exp-6: surgery-first off (ctskel) and on (ctskel-sf), both orders
        for g in GRAPHS:
            for name, methods in (('sf-ab', 'ctskel,ctskel-sf'), ('sf-ba', 'ctskel-sf,ctskel')):
                job(f'{name}-{g}', EST[g], f'{ZB} causal --graph {G(g)} --zsizes 1,4,16 --zmodes random,ancestor --kinds rule2,rule3 '
                    f'--surgery-first 0 --seed 1 --queries 1000 --pairmode walk {COMMON} --methods {methods}')
    elif run.startswith('closure-ablation'):   # up-closures on and off, both orders (paper-data.json only)
        methods = 'ctskel-noclo,ctskel' if '-rev' in run else 'ctskel,ctskel-noclo'
        for g in GRAPHS:
            job(f'clo-{g}', EST[g], f'{ZB} bench --graph {G(g)} {DSEP} --methods {methods}')
    elif run == 'index-stats':   # shortcut-arc statistics of TrekIndex (full version, Exp-3); not timed
        cmd = ' && '.join(f'{ZB} stats --graph {G(g)} --progress 0 >> {out}/stats.jsonl' for g in GRAPHS)
        job('stats', 1, f'rm -f {out}/stats.jsonl && {cmd}')
    else:
        sys.exit(f'unknown run {run}')
    return jobs


def compare(run, ref):
    """workloads-vs-e2.txt: one line 'same NAME' or 'differ NAME' per d-separation workload file of the run."""
    run, ref = pathlib.Path(run), pathlib.Path(ref)
    lines = []
    for f in sorted((run / 'workloads').glob('*-s1.*')):
        r = ref / 'workloads' / f.name
        lines.append(f'{"same" if r.exists() and filecmp.cmp(f, r, shallow=False) else "differ"} {f.name}')
    (run / 'workloads-vs-e2.txt').write_text('\n'.join(lines) + '\n')
    print(f'{run}/workloads-vs-e2.txt: {sum(l.startswith("same") for l in lines)} of {len(lines)} workloads identical')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    j = sub.add_parser('jobs')
    j.add_argument('--data', required=True, type=pathlib.Path)
    j.add_argument('--results', required=True, type=pathlib.Path)
    j.add_argument('--runs', default=','.join(RUNS))
    c = sub.add_parser('compare')
    c.add_argument('run'); c.add_argument('ref')
    a = ap.parse_args()
    if a.cmd == 'compare':
        return compare(a.run, a.ref)
    data, res = a.data.resolve(), a.results.resolve()
    missing = [g for g in GRAPHS if not (data / f'{g}.bin').exists()]
    if missing:
        print(f'warning: {data} lacks {", ".join(m + ".bin" for m in missing)}', file=sys.stderr)
    for run in a.runs.split(','):
        jobs = run_jobs(run, data, res)
        (res / run).mkdir(parents=True, exist_ok=True)
        (res / run / 'jobs.jsonl').write_text(''.join(json.dumps(x) + '\n' for x in jobs))
        print(f'{res / run}/jobs.jsonl: {len(jobs)} jobs')


if __name__ == '__main__':
    main()
