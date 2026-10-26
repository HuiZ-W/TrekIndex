#!/usr/bin/env python3
"""Generate every experiment number of the paper (tables, plot data, text macros) from the result files.

usage: python3 scripts/make_paper_artifacts.py [--results DIR] [--out DIR] [--replay-stats FILE]

  --results       directory with one sub-directory per experiment run (default: results/, or $TREKINDEX_RESULTS)
  --out           directory for the generated files (default: gen/, or $TREKINDEX_GEN_OUT)
  --replay-stats  causenet-mediation.stats.json written by scripts/build_causenet_mediation.py
                  (default: <results>/causenet-mediation.stats.json)

Every job of a run writes <name>.stdout (JSON lines), <name>.stderr and <name>.status (exit code), as
scripts/run_jobs.py does; scripts/paper_jobs.py writes the job lists. G is one of the six graph keys below.
Inputs, by experiment of the paper (Section 7):
  Table 3, Exp-3   round2/sweep-s{1,2,3}-G (round2-extra/ for iKraph): n, m, build time (mean of the three
                   builds), label entries, index size, width, core and peeled vertices; the plain 2-hop index
                   ("zskel") from the same jobs; round2/markov-G (round2-extra/ for iKraph; peak RSS, written
                   to paper-data.json only); index-stats/stats.jsonl (shortcut arcs, appendix).
  Exp-1            e2-merged/bench-s1-G (BB, Anc-BB, CIfly, TrekIndex) and e2-merged/rp-s1-G (RP), seed 1;
                   classic-v1/dsep-G (Geiger, d-collision) and classic-v1/workloads-vs-e2.txt (same requests).
  Exp-2            criteria-all/{backdoor,rules,frontdoor,further}-G (every baseline and TrekIndex),
                   classic-v1/{backdoor,rules,frontdoor,further}-G (Geiger, d-collision); round2/causal-G
                   (round2-extra/ for iKraph) and round3/{rules,further}-s{1,2,3}-G (TrekIndex and BB).
  Exp-4            ablation-b/abl-G, ablation-b/workloads-vs-e2.txt, ablation-b/DONE.
  Exp-5            decomp-order{,-r2,-r3}/dc-{ab,ba}-G (plain 2-hop labels and TrekIndex, both orders).
  Exp-6            surgery-order{,-r2,-r3}/sf-{ab,ba}-G (surgery-first off and on, both orders);
                   closure-ablation{,-rev}{,-r2,-r3}/clo-G (up-closures; paper-data.json only).
  Exp-7            edge-tests/et-G, checked against ablation-final/{default,noclosure,fulllabels}-G.
  Exp-8            apps-all/{screen-G,med-causenet-observed,med-causenet-all}, classic-v1/{screen-G,
                   replay-causenet-observed,pc-G}, pc-v2/{pc-G.stdout,pcgen-G.json,sample-G.json}, the
                   replay statistics (--replay-stats), and exp5-markov/{markov,rp}-G (local-Markov tests,
                   not reported in the paper but required by the checks of the text macros).
"git <hash>" in the comments below names the commit of the development repository that produced a run.
The comments call the experiments by their working names (E1-E6, Exp-4 to Exp-7), which predate the
numbering of the paper.

Outputs (in --out): numbers.tex (text macros), tab-*.tex (tables), *.dat (plot data) and paper-data.json
(every normalized number used by the tables, plot data and macros, and the list of source files read).

Conventions: per cell, the mean of the per-seed means and the mean of the per-seed p99; speedup =
baseline mean / index mean; every cell on its own, no medians, no averages over cells (min/max over
cells only in text macros). A missing or failed job aborts the run, except the application macros,
which become a red ?? and are listed.
"""
import argparse, json, math, os, pathlib, re, statistics, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = pathlib.Path(os.environ.get('TREKINDEX_RESULTS', ROOT / 'results'))
OUT = pathlib.Path(os.environ.get('TREKINDEX_GEN_OUT', ROOT / 'gen'))
REPLAY_STATS = None   # --replay-stats; default RES / 'causenet-mediation.stats.json'
SCRIPT = 'scripts/make_paper_artifacts.py'

# file key, table name, text name, macro key, kind
GRAPHS = [
    ('ceg-dag', 'CEG', 'CEG', 'CEG', 'causal'),
    ('syn-ba-500k', 'BA 500k', 'BA 500k', 'BAfive', 'synthetic'),
    ('aser-concept-dag', 'ASER', 'ASER', 'ASER', 'causal'),
    ('ikraph-causal-dag', 'iKraph', 'iKraph', 'iKraph', 'causal'),
    ('syn-ba-2000k', 'BA 2M', 'BA 2M', 'BAtwo', 'synthetic'),
    ('dcol-cnfdag', 'CNFDAG', 'CNFDAG', 'CNFDAG', 'causal'),
    ('causenet-full-dag', 'CauseNet', 'CauseNet', 'CauseNet', 'causal'),
    ('twitter', 'Twitter (s)', 'Twitter', 'Twitter', 'structural'),
]
# CEG and Twitter (s) were benchmarked during development but are not in the paper; they are not read.
DROPPED = ('ceg-dag', 'twitter')
GRAPHS = [g for g in GRAPHS if g[0] not in DROPPED]
KEYS = [g[0] for g in GRAPHS]
TAB = {g[0]: g[1] for g in GRAPHS}
TXT = {g[0]: g[2] for g in GRAPHS}
MAC = {g[0]: g[3] for g in GRAPHS}
KIND = {g[0]: g[4] for g in GRAPHS}
# domain column of the dataset table (CNFDAG and CauseNet both come from the CauseNet claims)
DOMAIN = {'ceg-dag': 'Lexical', 'syn-ba-500k': 'Random graph', 'aser-concept-dag': 'Event',
          'ikraph-causal-dag': 'Biomedical', 'syn-ba-2000k': 'Random graph', 'dcol-cnfdag': 'Web claim',
          'causenet-full-dag': 'Web claim', 'twitter': 'Social'}

# which round2 directory holds a graph's jobs
SWEEP_DIR = {'ikraph-causal-dag': 'round2-extra', 'twitter': 'round2-extra'}
CAUSAL_DIR = {'ikraph-causal-dag': 'round2-extra', 'twitter': 'round2d'}
MARKOV_DIR = {'ikraph-causal-dag': 'round2-extra'}

SEEDS = (1, 2, 3)
E2_SEEDS = (1,)   # e2-merged also ran seeds 2 and 3, which were partly slowed by other load on the machine; seed 1 is used
E2_SEEDS_TEXT = 'seed 1'
ZMODES = ('random', 'ancestor')
ZSIZES = (1, 2, 4, 8, 16, 32)
# The six traversal baselines. Geiger and d-collision (the d-collision graph method, Neo4j UDF variant, ported to
# C++) were added last; their times come from results/classic-v1 (every workload of the paper, Bidirectional
# first in each process as the answer reference, 12 jobs at a time).
# In the order of the Competitors list of the paper (classical procedure first), which the table columns follow.
BASELINES = [('geiger', 'Geiger'), ('bidirectional', 'Bi'), ('bidirectional-anc', 'AnBi'), ('rp', 'RP'),
             ('cifly', 'CIfly'), ('dcollision', 'DCol')]
CLASSIC = ('geiger', 'dcollision')
CLASSIC_RUN = 'classic-v1'
OLD_BASELINES = [b for b, _ in BASELINES if b not in CLASSIC]   # the baselines of the paper's own runs
BASE_TEXT = {'bidirectional': 'bidirectional Bayes-Ball', 'bidirectional-anc': 'ancestral Bayes-Ball',
             'cifly': 'CIfly', 'rp': 'RP-separation', 'geiger': 'Geiger', 'dcollision': 'd-collision'}
COL = {'bidirectional': 'bidir', 'bidirectional-anc': 'ancbidir', 'cifly': 'cifly', 'rp': 'rp',
       'geiger': 'geiger', 'dcollision': 'dcol'}   # columns of the critbar data files
# classic-v1 job of each application workload of apps-all
CLASSIC_APP = {'screen-causenet-full-dag': 'screen-causenet-full-dag', 'screen-aser-concept-dag': 'screen-aser-concept-dag',
               'med-causenet-observed': 'replay-causenet-observed'}
# classic-v1 jobs still running: the methods not yet finished, with how long (s) the pending method had run without
# finishing its warm-up pass. The timed pass repeats the warm-up pass, so a method that has run longer than twice the
# fastest finished baseline cannot be the fastest; its cell is printed as a red ?? until the job finishes, and the
# speedup is over the finished baselines (checked in load_apps). Remove an entry once its job has a status file.
CLASSIC_PENDING = {}   # empty: every classic-v1 job has finished
PC_RUN = 'pc-v2'   # Exp-7 PC workload: 20 ancestral regions of 800-1000 vertices, uniform sample of 10,000 tests
APP_TABLE = (('screen', ('causenet-full-dag', 'aser-concept-dag')), ('replay-observed', ('causenet-full-dag',)))
CRITERIA = [  # kind, macro key, text name, table name, Z modes, source
    ('backdoor', 'Backdoor', 'back-door', 'Back-door', ('parents', 'random', 'ancestor'), 'round2'),
    ('rule1', 'RuleOne', 'rule 1', 'Rule 1', ('random', 'ancestor'), 'round2'),
    ('rule2', 'RuleTwo', 'rule 2', 'Rule 2', ('random', 'ancestor'), 'round3-rules'),
    ('rule3', 'RuleThree', 'rule 3', 'Rule 3', ('random', 'ancestor'), 'round3-rules'),
    ('frontdoor', 'Frontdoor', 'front-door', 'Front-door', ('children', 'random', 'ancestor'), 'round3-further'),
    ('seqbackdoor', 'SeqBackdoor', 'sequential back-door', 'Seq.\\ back-door', ('random', 'ancestor'), 'round3-further'),
    ('mediation', 'Mediation', 'mediation', 'Mediation', ('random', 'ancestor'), 'round3-further'),
    ('selection', 'Selection', 'selection', 'Selection', ('random', 'ancestor'), 'round3-further'),
]
CRIT_K = (1, 4, 16)
# Figure 6 (bar chart, random sets, mean query time over the sizes CRIT_K) and the appendix figures
# (ancestor sets, and every size separately): every graph and every method. Taken from criteria-all
# (every baseline, git 3d16980) where the graph's job has finished, otherwise from the older runs
# above, which time Bidirectional ("traversal") and the index only. Exp-2 reports the mean over the
# sizes; the per-size cells stay in the appendix and in CritCell* macros.
CRIT_ALL_JOB = {'backdoor': 'backdoor', 'rule1': 'rules', 'rule2': 'rules', 'rule3': 'rules', 'frontdoor': 'frontdoor',
                'seqbackdoor': 'further', 'mediation': 'further', 'selection': 'further'}
QUERIES_PER_RUN = 1000
RED = '\\textcolor{red}{??}'


class SourceError(Exception):
    pass


# ------------------------------------------------------------------ reading
USED = set()   # every source file read, for paper-data.json


def rel(p):
    p = pathlib.Path(p)
    return 'results/' + str(p.relative_to(RES)) if p.is_relative_to(RES) else str(p)


def rows(path):
    """JSON lines of a job's stdout; the job's .status must exist and report exit 0."""
    path = pathlib.Path(path)
    status = path.with_suffix('.status')
    if not path.exists():
        raise SourceError(f'missing source {rel(path)}')
    if not status.exists():
        raise SourceError(f'missing status {rel(status)}')
    st = json.loads(status.read_text().strip().splitlines()[-1])
    if st.get('exit') != 0:
        raise SourceError(f'{rel(status)}: exit {st.get("exit")}')
    out = [json.loads(re.sub(r':(-?nan|-?inf)\b', ':null', l)) for l in path.read_text().splitlines() if l.startswith('{')]
    if not out:
        raise SourceError(f'{rel(path)}: no result rows')
    USED.add(rel(path)); USED.add(rel(status))
    return out


def meth(r):
    return {m['method']: m for m in r['methods']}


def mean(v):
    return statistics.mean(list(v))


def check(cond, msg):
    if not cond:
        raise SourceError(msg)


# --------------------------------------------------------------- formatting
def num(x):
    """Integer with LaTeX thousands separators."""
    return f'{int(math.floor(x + 0.5)):,}'.replace(',', '{,}')


def sig(x):
    """3 significant digits below 100, integer with separators from 100 on."""
    if x >= 99.95:
        return num(x)
    if x == 0:
        return '0'
    r = float(f'{x:.3g}')
    return f'{r:.{max(0, 2 - math.floor(math.log10(r)))}f}'


def spd(x):
    """Speedup for tables: bold below 1."""
    return f'\\textbf{{{sig(x)}}}' if x < 1 else sig(x)


def tus(us):
    """Microseconds -> LaTeX time with unit, 3 significant digits."""
    r = float(f'{us:.3g}')
    if r >= 1e6:
        return f'{sig(us / 1e6)}\\,s'
    if r >= 1e3:
        return f'{sig(us / 1e3)}\\,ms'
    return f'{sig(us)}\\,\\textmu s'


def tsec(s):
    return tus(s * 1e6)


def nbytes(b):
    r = b / 2**30
    return f'{sig(r)}\\,GiB' if r >= 0.9995 else f'{sig(b / 2**20)}\\,MiB'


def count(x):
    """Large count -> compact (k / M)."""
    if x >= 1e6:
        return f'{sig(x / 1e6)}\\,M'
    if x >= 1e4:
        return f'{sig(x / 1e3)}\\,k'
    return num(x)


def words(x):
    if x >= 1e9:
        return f'{sig(x / 1e9)} billion'
    if x >= 1e6:
        return f'{sig(x / 1e6)} million'
    return num(x)


def header(sources):
    return f'% Generated by {SCRIPT} from {sources}; do not edit.'


# ======================================================================= E1
def load_e1():
    e1 = {}
    for g in KEYS:
        d = RES / SWEEP_DIR.get(g, 'round2')
        seeds = [rows(d / f'sweep-s{s}-{g}.stdout') for s in SEEDS]
        ct = [meth(r[0])['ctskel'] for r in seeds]
        zs = [meth(r[0]).get('zskel') for r in seeds]
        for f in ('label_entries', 'index_bytes', 'fringe_bytes', 'width', 'core_n', 'peeled'):
            check(len({c[f] for c in ct}) == 1, f'{g}: ctskel {f} differs across seeds')
        n, m = seeds[0][0]['n'], seeds[0][0]['m']
        c = ct[0]
        check(c['core_n'] + c['peeled'] == n, f'{g}: core + peeled != n')
        rec = {'n': n, 'm': m, 'kind': KIND[g], 'source': rel(d / f'sweep-s[123]-{g}.stdout'),
               'build_s_per_seed': [x['build_ms'] / 1e3 for x in ct],
               'build_s': mean(x['build_ms'] / 1e3 for x in ct),
               'peel_s': mean(x['peel_ms'] / 1e3 for x in ct),
               'label_entries': c['label_entries'], 'label_bytes': c['index_bytes'], 'fringe_bytes': c['fringe_bytes'],
               'index_bytes': c['index_bytes'] + c['fringe_bytes'], 'width': c['width'], 'core_n': c['core_n'],
               'core_m': c['core_m'], 'peeled': c['peeled'], 'peeled_frac': c['peeled'] / n, 'shortcuts': c['shortcuts'],
               'peak_rss_with_plain_bytes': max(r['peak_rss_bytes'] for s in seeds for r in s)}
        if all(zs):
            check(len({z['label_entries'] for z in zs}) == 1, f'{g}: zskel entries differ across seeds')
            rec['plain'] = {'build_s_per_seed': [z['build_ms'] / 1e3 for z in zs], 'build_s': mean(z['build_ms'] / 1e3 for z in zs),
                            'label_entries': zs[0]['label_entries'], 'label_bytes': zs[0]['index_bytes']}
        else:
            check(not any(zs), f'{g}: zskel present in some seeds only')
            # round2 has no plain 2-hop run for Twitter; the ablation run (git 72da81c, one build) has
            ap = RES / 'ablation-final' / f'fulllabels-{g}.stdout'
            ar = rows(ap)
            az = [meth(r)['zskel'] for r in ar]
            check(all(r['n'] == n and r['m'] == m for r in ar), f'{g}: n/m differ in {rel(ap)}')
            check(len({z['build_ms'] for z in az}) == 1 and len({z['label_entries'] for z in az}) == 1,
                  f'{g}: {rel(ap)} built more than once')
            rec['plain'] = {'build_s_per_seed': [az[0]['build_ms'] / 1e3], 'build_s': az[0]['build_ms'] / 1e3,
                            'label_entries': az[0]['label_entries'], 'label_bytes': az[0]['index_bytes'],
                            'source': rel(ap)}
        # peak RSS of a process with the graph, bidirectional traversal and the core-tree index only
        if g == 'twitter':
            mp, mrows = d / f'sweep-s[123]-{g}.stdout', [r for s in seeds for r in s]
        else:
            mp = RES / MARKOV_DIR.get(g, 'round2') / f'markov-{g}.stdout'
            mrows = rows(mp)
        check(all(set(meth(r)) == {'bidirectional', 'ctskel'} for r in mrows), f'{g}: peak-RSS job has other methods')
        check(all(r['n'] == n for r in mrows), f'{g}: n differs in the peak-RSS job')
        rec['peak_rss_bytes'] = max(r['peak_rss_bytes'] for r in mrows)
        rec['peak_rss_source'] = rel(mp)
        e1[g] = rec
    return e1


def sweeps(g):
    """round2 d-separation sweep of graph g: {seed: {(zmode, k): row}}, every cell present."""
    d = RES / SWEEP_DIR.get(g, 'round2')
    out = {}
    for s in SEEDS:
        rs = {}
        for r in rows(d / f'sweep-s{s}-{g}.stdout'):
            check(r['seed'] == s and r['queries'] == QUERIES_PER_RUN and r['pairmode'] == 'walk', f'{g} s{s}: unexpected sweep row')
            rs[(r['zmode'], r['zsize'])] = r
        check(sorted(rs) == sorted((zm, k) for zm in ZMODES for k in ZSIZES), f'{g} s{s}: sweep cells {sorted(rs)}')
        out[s] = rs
    return out


# ============================================== E6: core-tree vs plain 2-hop query time
def load_ablation():
    """Exp-4: mean query time of the plain 2-hop index on the full trek DAG ("zskel") / TrekIndex
    ("ctskel") for d-separation. Both answer the same requests in one process, and the method timed
    second gains from warm caches (the round2 sweeps always timed zskel first), so every graph runs in
    both orders, decomp-order{REP}/dc-ab-G (plain first) and dc-ba-G (TrekIndex first), three times each
    (REPS), git 88b068a, seed 1; the ratio is the geometric mean over the two orders of the per-order
    median of the three runs, as in load_surgery. On graphs without a periphery (BA) the two indexes
    coincide, and the ratios above 1 at large ancestor sets come from the up-closures, which only
    TrekIndex has (compare load_closure)."""
    abl = {}
    for g in KEYS:
        runs = {}
        for rep in REPS:
            D = RES / ('decomp-order' + rep)
            for order, job, first in (('plain_first', 'dc-ab', 'zskel'), ('trek_first', 'dc-ba', 'ctskel')):
                for r in rows(D / f'{job}-{g}.stdout'):
                    ms = meth(r)
                    check(list(ms)[0] == first, f'{D.name}/{job}-{g}: unexpected method order {list(ms)}')
                    check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{D.name}/{job}-{g}: mismatches')
                    runs.setdefault((r['zmode'], r['zsize']), {}).setdefault(order, []).append(
                        (ms['zskel']['mean_us'], ms['ctskel']['mean_us']))
        for zm in ZMODES:
            for k in ZSIZES:
                v = runs.get((zm, k), {})
                check(set(v) == {'plain_first', 'trek_first'} and all(len(x) == len(REPS) for x in v.values()),
                      f'{g} {zm}-{k}: missing order or repetition in decomp-order')
                a = statistics.median(z / c for z, c in v['plain_first'])
                b = statistics.median(z / c for z, c in v['trek_first'])
                allz = [z for o in v.values() for z, _ in o]; allc = [c for o in v.values() for _, c in o]
                abl[(g, zm, k)] = {'graph': g, 'zmode': zm, 'k': k, 'ratio': math.sqrt(a * b),
                                   'ratio_plain_first': a, 'ratio_trek_first': b,
                                   'zskel_mean_us': statistics.median(allz), 'ctskel_mean_us': statistics.median(allc),
                                   'zskel_mismatches': 0,
                                   'source': rel(RES / 'decomp-order{,-r2,-r3}' / f'dc-{{ab,ba}}-{g}.stdout')}
    return abl


# ===================================================================== index ablation
# Exp-4, the trek-DAG reduction without an index:
# results/ablation-b, code-ablB = git 15b1bff + ABLATION_B.diff (TrekSearch). Bench methods: "trek-bidir"
# (trek traversal: one bidirectional search per query in T without the states of Z vertices, plus an arc
# E_z -> P_z for each z in Z, i.e. Bayes-Ball on the four states of T) and "trek-edge" (trek search: H_Z explored as
# by TrekIndex, each edge test a bidirectional search in T). One process per graph runs bidirectional, ctskel,
# trek-bidir and trek-edge, under the conditions of e2-merged (Exp-1 protocol, seed 1, jobs one after another on
# CPU 1), and the run checked that its 72 workloads are byte-identical to those of e2-merged. The first run,
# results/ablation-a (six jobs at once, with the unidirectional "trek-sweep"), is superseded.
ABLA_RUN = 'ablation-b'


def load_index_ablation():
    D = RES / ABLA_RUN
    check((D / 'DONE').exists(), f'{rel(D)}: run not finished')
    same = [l for l in (D / 'workloads-vs-e2.txt').read_text().splitlines() if l.strip()]
    check(len(same) == len(KEYS) * len(ZMODES) * len(ZSIZES) and all(l.startswith('same ') for l in same),
          f'{rel(D)}: workloads differ from e2-merged')
    USED.add(rel(D / 'workloads-vs-e2.txt'))
    cells = {}
    for g in KEYS:
        e = {(r['zmode'], r['zsize']): r for r in rows(D / f'abl-{g}.stdout')}
        check(set(e) == {(zm, k) for zm in ZMODES for k in ZSIZES}, f'{rel(D)}: configurations of {g}')
        for key, r in e.items():
            m = meth(r)
            check([x['method'] for x in r['methods']] == ['bidirectional', 'ctskel', 'trek-bidir', 'trek-edge'],
                  f'{g} {key}: methods and order')
            check(all(x['mismatches_vs_first'] == 0 for x in m.values()), f'{g} {key}: answers differ')
            check(r['queries'] == 1000, f'{g} {key}: query count')
            cells[(g,) + key] = {'bi_us': m['bidirectional']['mean_us'], 'idx_us': m['ctskel']['mean_us'],
                                 'trav_us': m['trek-bidir']['mean_us'], 'edge_us': m['trek-edge']['mean_us'],
                                 'tests': m['trek-edge']['trek_searches'] / r['queries']}
    return cells


ABLA_TAB_K = (1, 32)   # |Z| of the configurations in Table 6; the Exp-4 ranges are taken over these, so the table shows them


def tab_index_ablation(ab):
    """Exp-4 table (main text): mean d-separation query time of Bidirectional, trek traversal (trek-bidir), trek search
    (trek-edge) and TrekIndex, measured in the same process, for |Z| = 1 and 32 with random and ancestor sets.
    The fastest method of each row and |Z| is bold on the band of TrekIndex's column in the appendix tables
    (\\trekcell, figures/heat_macros.tex). Appendix: every configuration, with the edge tests of trek search per query."""
    meths = [('bi_us', FULL_ABBR['bidirectional'], 1e3), ('trav_us', 'Trav.', 1e3), ('edge_us', 'Search', 1e3), ('idx_us', '\\Method{}', 1)]
    unit = lambda d: 'ms' if d == 1e3 else '\\textmu s'
    L = [header(f'results/{ABLA_RUN} (abl-*; code-ablB, seed 1)'),
         '\\input{figures/heat_macros}',
         '\\begin{table}[t]\\centering\\footnotesize',
         '\\caption{Mean query time with and without the index.}',
         '\\label{tab:index-abl}', '\\setlength{\\tabcolsep}{2pt}',
         colwidth_tabular('ll' + 'r' * (2 * len(meths))), '\\toprule',
         ' & '.join(['&'] + [f'\\multicolumn{{{len(meths)}}}{{c}}{{$|Z|={k}$}}' for k in ABLA_TAB_K]) + '\\\\',
         f'\\cmidrule(lr){{3-{2 + len(meths)}}}\\cmidrule(lr){{{3 + len(meths)}-{2 + 2 * len(meths)}}}',
         'Graph & $Z$ & ' + ' & '.join([t for _, t, _ in meths] * 2) + '\\\\',
         '& & ' + ' & '.join([unit(d) for _, _, d in meths] * 2) + '\\\\', '\\midrule']
    for gi, g in enumerate(KEYS):
        if gi:
            L.append('\\addlinespace[2pt]')
        for zi, (zm, zt) in enumerate((('random', 'Rand.'), ('ancestor', 'Anc.'))):
            vals = []
            for k in ABLA_TAB_K:
                c = ab[(g, zm, k)]
                best = min(meths, key=lambda m: c[m[0]])[0]   # all four times are in microseconds here
                vals += [bestcell(sig(c[f] / d)) if f == best else sig(c[f] / d) for f, _, d in meths]
            L.append(f'{TAB[g] if zi == 0 else ""} & {zt} & ' + ' & '.join(vals) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-index-abl.tex', L)
    head = [f'& & \\multicolumn{{{len(meths) + 1}}}{{c}}{{Random $Z$}} & \\multicolumn{{{len(meths) + 1}}}{{c}}{{Ancestor $Z$}}\\\\',
            f'\\cmidrule(lr){{3-{3 + len(meths)}}}\\cmidrule(lr){{{4 + len(meths)}-{4 + 2 * len(meths)}}}',
            '& & ' + ' & '.join(([t for _, t, _ in meths] + ['Tests']) * 2) + '\\\\',
            'Graph & $|Z|$ & ' + ' & '.join(([unit(d) for _, _, d in meths] + ['']) * 2) + '\\\\']
    rws = []
    for gi, g in enumerate(KEYS):
        if gi:
            rws.append('\\midrule')
        for i, k in enumerate(ZSIZES):
            vals = []
            for zm in ZMODES:
                c = ab[(g, zm, k)]
                vals += [sig(c[f] / d) for f, _, d in meths] + [f'{c["tests"]:.1f}']
            rws.append(f'{TAB[g] if i == 0 else ""} & {k} & ' + ' & '.join(vals) + '\\\\')
    n = len(meths) + 1
    write('tab-index-abl-full.tex', longtable(header(f'results/{ABLA_RUN} (abl-*; code-ablB, seed 1)'),
          '@{}lr' + 'r' * n + '@{\\hspace{14pt}}' + 'r' * n + '@{}', 2 * n + 2,
          'Mean d-separation query time of BB, trek traversal, trek search, and \\Method{}, and edge tests of '
          'trek search per query.',
          'tab:index-abl-full', head, rws, sep='6pt'))


def index_ablation_macros(ab):
    """Exp-4 text on the index ablation; every claim is checked here. The ranges are taken over the configurations of
    Table 6 (|Z| in ABLA_TAB_K), so the reader can find both ends in the table; the claims on every configuration
    (trek search is slower than BB and its slowdown grows with |Z|, so do the edge tests) are checked on all 72."""
    M = []

    def add(name, val, raw=None, note=''):
        check(re.fullmatch(r'[A-Za-z]+', name), f'bad macro name {name}')
        M.append((name, val, raw, note))
    v = [c for key, c in ab.items() if key[2] in ABLA_TAB_K]
    rng = lambda f: (min(f(c) for c in v), max(f(c) for c in v))
    lo, hi = rng(lambda c: c['trav_us'] / c['bi_us'])
    add('AblaTravMin', sig(lo), lo, 'trek traversal / BB, configurations of Table 6'); add('AblaTravMax', sig(hi), hi)
    check(lo > 0.75 and hi < 2.5, 'Exp-4 text: trek traversal is within a small factor of Bidirectional')
    check(all(c['edge_us'] > c['bi_us'] for c in ab.values()), 'Exp-4 text: trek search is slower than BB in every configuration')
    lo, hi = rng(lambda c: c['edge_us'] / c['bi_us'])
    add('AblaEdgeMin', sig(lo), lo, 'trek search / BB'); add('AblaEdgeMax', sig(hi), hi)
    for zm in ZMODES:
        for g in KEYS:
            seq = [ab[(g, zm, k)]['edge_us'] / ab[(g, zm, k)]['bi_us'] for k in ZSIZES]
            check(all(x < y for x, y in zip(seq, seq[1:])), f'Exp-4 text: trek search over Bidirectional grows with |Z| ({g}, {zm})')
    lo, hi = rng(lambda c: c['trav_us'] / c['idx_us'])
    add('AblaIdxTravMin', sig(lo), lo, 'trek traversal / TrekIndex'); add('AblaIdxTravMax', sig(hi), hi)
    lo, hi = rng(lambda c: c['edge_us'] / c['idx_us'])
    add('AblaIdxEdgeMin', sig(lo), lo, 'trek search / TrekIndex'); add('AblaIdxEdgeMax', sig(hi), hi)
    t1 = [ab[(g, zm, 1)]['tests'] for g in KEYS for zm in ZMODES]
    t32 = [ab[(g, zm, 32)]['tests'] for g in KEYS for zm in ZMODES]
    check(max(t1) < min(t32), 'Exp-4 text: a query performs more edge tests at |Z| = 32 than at |Z| = 1 on every graph')
    return M


# ===================================================================== edge tests
ET_FIG_K = (1, 8, 32)   # |Z| of the bars of Figure 7
ET_CLASSES = ('et_proj', 'et_refuted', 'et_certified', 'et_searched', 'et_closure')


def load_edgetests():
    """How the core-tree index ("ctskel") decides each H_Z edge test, from results/edge-tests
    (git 336c6c6, which only adds the counters et_* to the round2 code), seed 1: the requests of E2 seed 1
    and of ablation-final. edge() classifies each test by its costliest step:
      et_searched:  at least one label-pruned bidirectional search;
      et_certified: no search, a core test certified by the labels;
      et_closure:   no search or certificate, answered by an up-closure lookup (edge);
      et_refuted:   no search, certificate, or lookup; every core test refuted by the labels (no edge);
      et_proj:      no core test: a junction in the projections (edge), no live state, or no usable
                    pair of core states (no edge).
    Checks: all other counters and the answers equal those of ablation-final default-<g> (same requests,
    git 72da81c), and the five classes sum to the edge_checks of the plain 2-hop index (ablation-final
    fulllabels-<g>), which explores H_Z in the same order with exact edge answers."""
    D = RES / 'edge-tests'; A = RES / 'ablation-final'
    out = {}
    for g in KEYS:
        et = {(r['zmode'], r['zsize']): meth(r)['ctskel'] for r in rows(D / f'et-{g}.stdout')}
        de = {(r['zmode'], r['zsize']): meth(r)['ctskel'] for r in rows(A / f'default-{g}.stdout')}
        fl = {(r['zmode'], r['zsize']): meth(r)['zskel'] for r in rows(A / f'fulllabels-{g}.stdout')}
        for zm in ZMODES:
            for k in ZSIZES:
                c, d, z = et[(zm, k)], de[(zm, k)], fl[(zm, k)]
                for f in COUNTERS + ('separated',):
                    check(c[f] == d[f], f'{g} {zm}-{k}: {f} differs from ablation-final default ({c[f]} vs {d[f]})')
                check(z['separated'] == c['separated'], f'{g} {zm}-{k}: plain and core-tree answers differ')
                n = {f: c[f] for f in ET_CLASSES}
                tests = sum(n.values())
                check(tests == z['edge_checks'], f'{g} {zm}-{k}: {tests} classified tests, plain index {z["edge_checks"]}')
                q = QUERIES_PER_RUN
                core = c['neg_reach'] + c['cert_free'] + c['searches']
                out[(g, zm, k)] = {'graph': g, 'zmode': zm, 'k': k, 'queries': q, 'counts': n, 'tests': tests,
                                   'tests_per_query': tests / q, 'share': {f[3:]: v / tests for f, v in n.items()},
                                   'core_tests_per_query': core / q, 'searches_per_query': c['searches'] / q,
                                   'refuted_share_of_core_tests': c['neg_reach'] / core if core else None,
                                   'source': rel(D / f'et-{g}.stdout')}
    return out


def dat_edgetests(et):
    """gen/edgetests-<zmode>.dat for Figure 7: one row per graph and |Z| in ET_FIG_K, at x = 4 * graph index + the
    index of |Z|, with the share of the H_Z edge tests decided by each class. Up-closure lookups (at most 1% of the
    tests in every cell, checked) are counted with the certified tests: both find a trek without a search."""
    for zm in ZMODES:
        L = [header('results/edge-tests (et-*, git 336c6c6, seed 1)'),
             f'% d-separation, {zm} Z: share of the H_Z edge tests by how the index decides them; definitions in',
             '% load_edgetests(). certified includes up-closure lookups. tests: edge tests per query.',
             'x graph k proj refuted certified searched tests']
        for gi, g in enumerate(KEYS):
            for ki, k in enumerate(ET_FIG_K):
                r = et[(g, zm, k)]; sh = r['share']
                check(sh['closure'] <= 0.01, f'{g} {zm}-{k}: up-closure lookups decide {sh["closure"]:.3f} of the tests')
                L.append(f'{4 * gi + ki} {MAC[g]} {k} {sh["proj"]:.5f} {sh["refuted"]:.5f} '
                         f'{sh["certified"] + sh["closure"]:.5f} {sh["searched"]:.5f} {r["tests_per_query"]:.3f}')
        write(f'edgetests-{zm}.dat', L)


COUNTERS = ('edge_checks', 'neg_reach', 'cert_free', 'searches', 'search_found', 'fringe_resolved', 'core_routes',
            'closures', 'h_nodes')


# ======================================================================= E2
def load_e2():
    D = RES / 'e2-merged'
    C = RES / CLASSIC_RUN
    same = [l for l in (C / 'workloads-vs-e2.txt').read_text().splitlines() if l.strip()]
    check(len(same) == len(KEYS) * len(ZMODES) * len(ZSIZES) and all(l.startswith('same ') for l in same),
          f'{rel(C)}: d-separation workloads differ from e2-merged')
    USED.add(rel(C / 'workloads-vs-e2.txt'))
    raw = {}     # (g, zmode, k, method) -> {seed: (mean, p99)}
    runs = mism = 0
    rp_count_agree = rp_runs = 0
    nm = {}
    for g in KEYS:
        for s in E2_SEEDS:
            sep = {}
            for r in rows(D / f'bench-s{s}-{g}.stdout'):
                check(r['seed'] == s and r['queries'] == QUERIES_PER_RUN and r['pairmode'] == 'walk', f'{g} s{s}: unexpected bench row')
                nm[g] = (r['n'], r['m'])
                ms = meth(r)
                check(list(ms) == ['bidirectional', 'bidirectional-anc', 'cifly', 'ctskel'], f'{g}: unexpected E2 methods {list(ms)}')
                for name, m in ms.items():
                    raw.setdefault((g, r['zmode'], r['zsize'], name), {})[s] = (m['mean_us'], m['p99_us'])
                    if name != 'bidirectional':
                        runs += 1; mism += m['mismatches_vs_first']
                sep[(r['zmode'], r['zsize'])] = ms['bidirectional']['separated']
            for r in rows(D / f'rp-s{s}-{g}.stdout'):
                mm = re.match(r'(.*)-s(\d)\.(\w+)-(\d+)$', r['workload'])
                check(mm and mm.group(1) == g and int(mm.group(2)) == s, f'unexpected RP workload {r["workload"]}')
                check(r.get('mean_us') is not None and r['queries'] == QUERIES_PER_RUN, f'RP {r["workload"]}: no result')
                zm, k = mm.group(3), int(mm.group(4))
                raw.setdefault((g, zm, k, 'rp'), {})[s] = (r['mean_us'], r['p99_us'])
                rp_runs += 1; rp_count_agree += r['answers'].count('1') == sep[(zm, k)]
            if s == 1:   # Geiger and d-collision on the same workloads (classic-v1, Bidirectional as reference)
                for r in rows(C / f'dsep-{g}.stdout'):
                    check(r['seed'] == s and r['queries'] == QUERIES_PER_RUN and r['pairmode'] == 'walk', f'{g}: unexpected classic row')
                    ms = meth(r)
                    check(list(ms) == ['bidirectional', *CLASSIC], f'{g}: unexpected classic methods {list(ms)}')
                    check(ms['bidirectional']['separated'] == sep[(r['zmode'], r['zsize'])], f'{g}: classic answers differ from e2-merged')
                    for name in CLASSIC:
                        raw.setdefault((g, r['zmode'], r['zsize'], name), {})[s] = (ms[name]['mean_us'], ms[name]['p99_us'])
                        runs += 1; mism += ms[name]['mismatches_vs_first']
    cells = {}
    for g in KEYS:
        for zm in ZMODES:
            for k in ZSIZES:
                c = {'graph': g, 'zmode': zm, 'k': k}
                for name in [b for b, _ in BASELINES] + ['ctskel']:
                    v = raw.get((g, zm, k, name), {})
                    check(sorted(v) == list(E2_SEEDS), f'E2 {g} {zm}-{k} {name}: seeds {sorted(v)}')
                    c[name] = {'mean_us': mean(v[s][0] for s in E2_SEEDS), 'p99_us': mean(v[s][1] for s in E2_SEEDS),
                               'per_seed_mean_us': [v[s][0] for s in E2_SEEDS]}
                best = min((c[b]['mean_us'], i, b) for i, (b, _) in enumerate(BASELINES))
                c['fastest_baseline'] = best[2]
                c['speedup_vs_fastest'] = best[0] / c['ctskel']['mean_us']
                c['speedup_vs_bidirectional'] = c['bidirectional']['mean_us'] / c['ctskel']['mean_us']
                c['p99_ratio_vs_best_baseline_p99'], c['best_p99_baseline'] = p99_ratio(c)
                c['p99_ratio_vs_fastest_by_mean'] = c[best[2]]['p99_us'] / c['ctskel']['p99_us']
                cells[(g, zm, k)] = c
    return {'cells': cells, 'nm': nm, 'runs': runs, 'mismatches': mism, 'rp_runs': rp_runs, 'rp_count_agree': rp_count_agree}


# ======================================================================= E3
def load_e3():
    cells = {}
    runs = mism = 0
    for g in KEYS:
        # back-door and rule 1, seed 1
        p = RES / CAUSAL_DIR.get(g, 'round2') / f'causal-{g}.stdout'
        for r in rows(p):
            if r['kind'] not in ('backdoor', 'rule1'):
                continue
            check(r['queries'] == QUERIES_PER_RUN and r['pairmode'] == 'walk', f'{rel(p)}: unexpected row')
            ms = meth(r); t, c = ms['traversal'], ms['ctskel']
            runs += 1; mism += c['mismatches_vs_first']
            key = (r['kind'], g, r['zmode'], r['k'])
            check(key not in cells, f'duplicate cell {key}')
            cells[key] = {'traversal_mean_us': t['mean_us'], 'index_mean_us': c['mean_us'], 'index_p99_us': c['p99_us'],
                          'speedup': t['mean_us'] / c['mean_us'], 'seeds': [1], 'source': rel(p)}
        # round 3, seeds 1-3
        raw = {}
        for prefix in ('rules', 'further'):
            for s in SEEDS:
                p = RES / 'round3' / f'{prefix}-s{s}-{g}.stdout'
                for r in rows(p):
                    check(r['queries'] == QUERIES_PER_RUN and r['pairmode'] == 'walk', f'{rel(p)}: unexpected row')
                    for name, m in meth(r).items():
                        raw.setdefault((r['kind'], r['zmode'], r['k'], name), {})[s] = (m['mean_us'], m['p99_us'])
                        if name != 'traversal':
                            runs += 1; mism += m['mismatches_vs_first']
        for (kind, zm, k, name), v in raw.items():
            if name == 'traversal':
                continue
            check(sorted(v) == list(SEEDS), f'round3 {g} {kind} {zm}-{k} {name}: seeds {sorted(v)}')
            t = raw[(kind, zm, k, 'traversal')]
            check(sorted(t) == list(SEEDS), f'round3 {g} {kind} {zm}-{k} traversal: seeds {sorted(t)}')
            tm, cm = mean(t[s][0] for s in SEEDS), mean(v[s][0] for s in SEEDS)
            rec = {'traversal_mean_us': tm, 'index_mean_us': cm, 'index_p99_us': mean(v[s][1] for s in SEEDS),
                   'speedup': tm / cm, 'seeds': list(SEEDS),
                   'source': f'results/round3/{"rules" if kind in ("rule2", "rule3") else "further"}-s[123]-{g}.stdout'}
            if kind in ('rule2', 'rule3'):
                cells[(kind + {'ctskel': '-off', 'ctskel-sf': ''}[name], g, zm, k)] = rec
            else:
                check(name == 'ctskel', f'unexpected method {name}')
                cells[(kind, g, zm, k)] = rec
    # every expected cell present
    for kind, _, _, _, zms, _ in CRITERIA:
        for g in KEYS:
            for zm in zms:
                for k in CRIT_K:
                    for v in ((kind, kind + '-off') if kind in ('rule2', 'rule3') else (kind,)):
                        check((v, g, zm, k) in cells, f'missing criteria cell {v} {g} {zm}-{k}')
    return {'cells': cells, 'runs': runs, 'mismatches': mism}


# ======================================================================= E5
def load_apps(missing):
    """Application workloads against every baseline, query time only (construction: Exp-3).
    Screening and replay: apps-all (one job per workload, all methods in one process). Model
    testing: exp5-markov, 16 buckets of 1,000 local-Markov requests per graph; bidirectional,
    ancestral bidirectional, CIfly and the index in one process (markov-G), RP-separation replays
    the same request files (rp-G). The fastest baseline is the one with the smallest total."""
    apps = []
    names = OLD_BASELINES

    def finish(rec, base_ms, index_ms):
        rec['baseline_query_s'] = {b: base_ms[b] / 1e3 for b in base_ms}
        rec['index_query_s'] = index_ms / 1e3
        fb = min(base_ms, key=lambda b: base_ms[b])
        rec['fastest_baseline'] = fb
        rec['query_speedup'] = base_ms[fb] / index_ms
        rec['missing'] = False
        return rec

    def job(name, label, workload, graph, nq_field):
        p = RES / 'apps-all' / f'{name}.stdout'
        try:
            rs = rows(p)
        except SourceError as e:
            missing.append(str(e)); apps.append({'name': name, 'label': label, 'workload': workload, 'graph': graph, 'missing': True})
            return
        check(len(rs) == 1, f'{name}: expected one row')
        r = rs[0]; ms = meth(r)
        check(all(b in ms for b in names + ['ctskel']), f'{name}: methods {list(ms)}')
        check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{name}: mismatches')
        nq = ms['ctskel'][nq_field] if nq_field == 'tests' else r['count']
        if nq_field == 'tests':
            check(all(ms[b]['tests'] == nq for b in names), f'{name}: test counts differ')
        rec = {'name': name, 'label': label, 'workload': workload, 'graph': graph, 'source': rel(p), 'queries': nq,
               'n': r['n'], 'm': r['m'], 'index_build_s': ms['ctskel']['build_ms'] / 1e3,
               'mismatches': sum(m['mismatches_vs_first'] for m in ms.values()),
               'observed_frac': r.get('observed_frac'), 'cap': r.get('cap')}
        base = {b: ms[b]['query_ms'] for b in names}
        if name in CLASSIC_APP:   # Geiger and d-collision: classic-v1, the same workload, Bidirectional as reference
            pc = RES / CLASSIC_RUN / f'{CLASSIC_APP[name]}.stdout'
            pending = CLASSIC_PENDING.get(CLASSIC_APP[name], {})
            try:
                rc = rows(pc)
                check(not pending, f'{rel(pc)}: finished, remove it from CLASSIC_PENDING')
            except SourceError as e:
                if not pending:
                    missing.append(str(e)); apps.append({'name': name, 'label': label, 'workload': workload, 'graph': graph, 'missing': True})
                    return
                # still running: the methods that have finished, from the partial JSON row
                txt = pc.read_text()
                rc = [{'count': int(re.search(r'"count":(\d+)', txt).group(1)),
                       'methods': [json.loads(x) for x in re.findall(r'\{"method":"[^"]+"[^{}]*\}', txt)]}]
                missing.append(f'{rel(pc)}: running, pending {sorted(pending)} (red ?? cells)')
            check(len(rc) == 1, f'{rel(pc)}: expected one row'); mc = meth(rc[0])
            check(list(mc) + sorted(pending) == ['bidirectional', *CLASSIC], f'{rel(pc)}: methods {list(mc)}')
            check(all(m['mismatches_vs_first'] == 0 for m in mc.values()), f'{rel(pc)}: mismatches')
            if nq_field == 'tests':
                check(all(mc[b]['tests'] == nq for b in mc), f'{rel(pc)}: test counts differ from {name}')
                check((mc['bidirectional']['valid_z1'], mc['bidirectional']['valid_z2']) == (ms['bidirectional']['valid_z1'], ms['bidirectional']['valid_z2']),
                      f'{rel(pc)}: answers differ from {name}')
            else:
                check(rc[0]['count'] == nq and mc['bidirectional']['valid'] == ms['bidirectional']['valid'], f'{rel(pc)}: answers differ from {name}')
            base.update({b: mc[b]['query_ms'] for b in CLASSIC if b not in pending}); rec['source_classic'] = rel(pc)
            if pending:
                fastest_s = min(base.values()) / 1e3
                check(all(t > 2 * fastest_s for t in pending.values()), f'{rel(pc)}: a pending method may still be the fastest')
                rec['pending'] = sorted(pending)
        apps.append(finish(rec, base, ms['ctskel']['query_ms']))

    for g in ('causenet-full-dag', 'aser-concept-dag'):
        job(f'screen-{g}', f'Screening {TXT[g]}', 'screen', g, 'tests')
    job('med-causenet-observed', 'Replay CauseNet (observed covariates)', 'replay-observed', 'causenet-full-dag', 'count')
    job('med-causenet-all', 'Replay CauseNet (all covariates)', 'replay-all', 'causenet-full-dag', 'count')
    for g in KEYS:
        p, prp = RES / 'exp5-markov' / f'markov-{g}.stdout', RES / 'exp5-markov' / f'rp-{g}.stdout'
        try:
            rs, rrp = rows(p), rows(prp)
        except SourceError as e:
            missing.append(str(e)); apps.append({'name': f'markov-{g}', 'workload': 'markov', 'graph': g, 'missing': True})
            continue
        # the iKraph job also replays three d-separation workloads (random-1/4/16); keep the Markov buckets
        bucket = re.compile(re.escape(g) + r'\.(random|walk)-(all|pa[0-9+-]+)$')
        rs = [r for r in rs if bucket.match(r['zmode'].rsplit('/', 1)[-1])]
        rrp = [x for x in rrp if bucket.match(x['workload'])]
        check(len(rs) == 16 and len(rrp) == 16, f'{g}: expected 16 Markov buckets, got {len(rs)} and {len(rrp)}')
        ms = [meth(r) for r in rs]
        check(all(list(m) == ['bidirectional', 'bidirectional-anc', 'cifly', 'ctskel'] for m in ms), f'{g}: Markov methods')
        check(all(m[x]['mismatches_vs_first'] == 0 for m in ms for x in m), f'{g}: Markov mismatches')
        check(all(r['queries'] == rs[0]['queries'] for r in rs), f'{g}: Markov buckets differ in size')
        wl = {r['zmode'].rsplit('/', 1)[-1]: r for r in rs}
        check(len(wl) == len(rs) == len(rrp) and all(x['workload'] in wl for x in rrp), f'{g}: RP buckets differ')
        for x in rrp:
            r = wl[x['workload']]
            check(x['queries'] == r['queries'], f'{g} {x["workload"]}: RP query count')
            check(x['answers'].count('1') == meth(r)['bidirectional']['separated'], f'{g} {x["workload"]}: RP answers disagree')
        base = {b: sum(m[b]['query_ms'] for m in ms) for b in names if b != 'rp'}
        base['rp'] = sum(x['query_ms'] for x in rrp)
        rec = {'name': f'markov-{g}', 'label': f'Model testing {TXT[g]}', 'workload': 'markov', 'graph': g,
               'source': rel(p), 'source_rp': rel(prp), 'queries': sum(r['queries'] for r in rs),
               'buckets': len(rs), 'queries_per_bucket': rs[0]['queries'], 'mismatches': 0}
        apps.append(finish(rec, base, sum(m['ctskel']['query_ms'] for m in ms)))
    # PC workload (Exp-7): per graph, the tests of PC-stable at level 1 on 20 ancestral regions of
    # 800-1000 vertices; all methods time the same uniform sample of 10,000 of these tests.
    for g in KEYS:
        p = RES / PC_RUN / f'pc-{g}.stdout'
        if g == 'twitter':
            continue   # Twitter has no ancestral set of 800 to 1,000 vertices
        try:
            rs = rows(p); gen = json.loads((RES / PC_RUN / f'pcgen-{g}.json').read_text().splitlines()[0])
            smp = json.loads((RES / PC_RUN / f'sample-{g}.json').read_text().splitlines()[0])
        except (SourceError, OSError, ValueError, IndexError) as e:
            missing.append(f'{PC_RUN}/pc-{g}: {e}')
            apps.append({'name': f'pc-{g}', 'workload': 'pc', 'graph': g, 'missing': True})
            continue
        check(len(rs) == 1, f'{p}: expected one row')
        r = rs[0]; ms = meth(r)
        check(list(ms) == ['ctskel', 'bidirectional-anc', 'rp', 'bidirectional', 'cifly'], f'{p}: methods {list(ms)}')
        check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{p}: mismatches')
        check(r['queries'] == smp['sample'] and smp['tests'] == gen['tests'], f'{p}: sample sizes')
        pc = RES / CLASSIC_RUN / f'pc-{g}.stdout'
        rc = rows(pc); check(len(rc) == 1, f'{rel(pc)}: expected one row'); mc = meth(rc[0])
        check(list(mc) == ['bidirectional', *CLASSIC] and all(m['mismatches_vs_first'] == 0 for m in mc.values()), f'{rel(pc)}: methods or answers')
        check(rc[0]['queries'] == r['queries'] and mc['bidirectional']['separated'] == ms['bidirectional']['separated'], f'{rel(pc)}: answers differ from {rel(p)}')
        rec = {'name': f'pc-{g}', 'label': f'PC {TXT[g]}', 'workload': 'pc', 'graph': g, 'source': rel(p), 'source_classic': rel(pc),
               'queries': r['queries'], 'pc_tests': gen['tests'], 'pc_separated': gen['separated'],
               'pc_capped_regions': gen['capped_regions'], 'pc_regions': gen['regions'],
               'pc_region_min': gen['min_size'], 'pc_region_max': gen['region_size'], 'pc_level': gen['max_level'],
               'mean_ancestral': ms['bidirectional-anc'].get('mean_ancestral'),
               'mean_us': {**{b: ms[b]['mean_us'] for b in names + ['ctskel']}, **{b: mc[b]['mean_us'] for b in CLASSIC}}, 'mismatches': 0}
        apps.append(finish(rec, {**{b: ms[b]['query_ms'] for b in names}, **{b: mc[b]['query_ms'] for b in CLASSIC}}, ms['ctskel']['query_ms']))
    return apps

# =================================================================== outputs
def write(name, lines):
    (OUT / name).write_text('\n'.join(lines) + '\n')


def colwidth_tabular(spec):
    """Begin a table that spans the column: the first column flush left, the extra space spread
    between the columns (tabular* with \\extracolsep{\\fill})."""
    return '\\begin{tabular*}{\\columnwidth}{@{}' + spec[0] + '@{\\extracolsep{\\fill}}' + spec[1:] + '@{}}'


def tab_datasets(e1):
    L = [header('results/round2 and round2-extra (sweep-s1-*: n, m)'),
         '\\begin{table}[t]\\centering\\small',
         '\\caption{Statistics of the datasets used in our experiments.}',
         '\\label{tab:datasets}', '\\setlength{\\tabcolsep}{4pt}',
         colwidth_tabular('lllrr'), '\\toprule',
         'Graph & Type & Domain & $|V|$ & $|E|$\\\\', '\\midrule']
    for g in KEYS:
        r = e1[g]
        L.append(f'{TAB[g]} & {r["kind"].capitalize()} & {DOMAIN[g]} & {num(r["n"])} & {num(r["m"])}\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-datasets.tex', L)


def tab_build(e1):
    L = [header('results/round2 and round2-extra (sweep-s[123]-*: width, core, peeled, build, entries, bytes)'),
         '\\begin{table}[t]\\centering\\small',
         '\\caption{Index construction of \\Method{} on one thread.}',
         '\\label{tab:build}', '\\setlength{\\tabcolsep}{4pt}',
         colwidth_tabular('lrrrrrr'), '\\toprule',
         'Graph & \\makecell[r]{Width\\\\$w$} & \\makecell[r]{Peeled\\\\vertices} & \\makecell[r]{Core\\\\vertices} & '
         '\\makecell[r]{Build\\\\time} & \\makecell[r]{Label\\\\entries} & \\makecell[r]{Index\\\\size}\\\\', '\\midrule']
    for g in KEYS:
        r = e1[g]
        L.append(f'{TAB[g]} & {r["width"]} & {100 * r["peeled_frac"]:.1f}\\,\\% & {count(r["core_n"])} & '
                 f'{tsec(r["build_s"])} & {count(r["label_entries"])} & {nbytes(r["index_bytes"])}\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-build.tex', L)


def tab_decomp(e1, abl):
    """Exp-4: plain 2-hop on the full trek DAG vs the core-tree index: build time, label entries, and
    the largest plain / core-tree mean query time over the 12 d-separation configurations."""
    L = [header('results/round2 and round2-extra (sweep-s[123]-*: build, entries), decomp-order{,-r2,-r3} '
                '(query time, both method orders)'),
         '\\begin{table}[t]\\centering\\small',
         '\\caption{Effect of the core--tree decomposition.}',
         '\\label{tab:decomp}', '\\setlength{\\tabcolsep}{3pt}',
         colwidth_tabular('lrrrrr'), '\\toprule',
         '& \\multicolumn{2}{c}{Build time} & \\multicolumn{2}{c}{Label entries} & \\makecell[r]{Max.\\ query}\\\\',
         '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}',
         'Graph & Plain & \\Method{} & Plain & \\Method{} & speedup\\\\', '\\midrule']
    for g in KEYS:
        r = e1[g]; p = r['plain']
        rs = [c['ratio'] for c in abl.values() if c['graph'] == g]
        check(len(rs) == len(ZMODES) * len(ZSIZES), f'{g}: ablation cells {len(rs)}')
        L.append(f'{TAB[g]} & {tsec(p["build_s"])} & {tsec(r["build_s"])} & {count(p["label_entries"])} & '
                 f'{count(r["label_entries"])} & {sig(max(rs))}$\\times$\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-decomp.tex', L)


def dat_decomp(e1, abl):
    """gen/decomp.dat for the Exp-4 dumbbell chart (figures/decomp_plot.tex, which replaces tab-decomp in the
    body): one row per graph in increasing order of the share of peeled vertices (row 0 at the bottom; ties keep the
    order of KEYS), with the label entries of plain 2-hop labels and of TrekIndex and the largest query speedup over
    plain labels (the same values as tab-decomp)."""
    L = [header('results/round2 and round2-extra (sweep-s[123]-*: label entries), decomp-order{,-r2,-r3} '
                '(query time, both method orders)'),
         '% label: y tick text; peeled: % of peeled vertices; plainentries/trekentries: label entries;',
         '% maxspeedup: largest plain / TrekIndex mean query time over the 12 d-separation configurations.',
         'row graph label peeled plainentries trekentries maxspeedup']
    for i, g in enumerate(sorted(KEYS, key=lambda g: e1[g]['peeled_frac'])):
        r = e1[g]
        rs = [c['ratio'] for c in abl.values() if c['graph'] == g]
        check(len(rs) == len(ZMODES) * len(ZSIZES), f'{g}: ablation cells {len(rs)}')
        label = TAB[g].replace(' ', '~') + f'~({round(100 * r["peeled_frac"])}\\,\\%)'
        L.append(f'{i} {MAC[g]} {label} {100 * r["peeled_frac"]:.1f} {r["plain"]["label_entries"]} '
                 f'{r["label_entries"]} {max(rs):.4g}')
    write('decomp.dat', L)


def tab_app_plain(e1):
    """Appendix (witness contraction), one column wide: build time, label entries and index size of 2-hop
    labels on the full trek DAG ("zskel" of the round2 sweeps, whose index is its labels alone) and of TrekIndex
    (Table 4). The query times of the two indexes are in Figure 7 (largest speedup) and not repeated here."""
    L = [header('results/round2 and round2-extra (sweep-s[123]-*: build, entries, bytes)'),
         '\\begin{table}[tp]\\centering\\footnotesize',
         '\\caption{Construction of 2-hop labels on the full trek DAG (Plain) and of \\Method{}.}',
         '\\label{tab:app:plain}', '\\setlength{\\tabcolsep}{2pt}',
         colwidth_tabular('lrrrrrr'), '\\toprule',
         '& \\multicolumn{2}{c}{Build time} & \\multicolumn{2}{c}{Label entries} & \\multicolumn{2}{c}{Index size}\\\\',
         '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}',
         'Graph & Plain & \\Method{} & Plain & \\Method{} & Plain & \\Method{}\\\\', '\\midrule']
    for g in KEYS:
        r = e1[g]; p = r['plain']
        L.append(f'{TAB[g]} & {tsec(p["build_s"])} & {tsec(r["build_s"])} & {count(p["label_entries"])} & '
                 f'{count(r["label_entries"])} & {nbytes(p["label_bytes"])} & {nbytes(r["index_bytes"])}\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-app-plain.tex', L)


def tab_dsep(e2):
    cells = e2['cells']; abbr = dict(BASELINES)
    L = [header('results/e2-merged (bench-s1-*, rp-s1-*)'),
         '\\begin{table*}[t]\\centering\\scriptsize',
         '\\caption{d-separation per cell: mean and p99 query time of the index, mean of the fastest baseline '
         '(Bi: bidirectional; AnBi: ancestral bidirectional; CIfly: C++ port; RP: RP-separation) and speedup '
         '= fastest-baseline mean / index mean. Seed 1, 1{,}000 walk-pair queries per cell; '
         '(s): structural, non-causal. Bold: below 1.}',
         '\\label{tab:dsep}', '\\setlength{\\tabcolsep}{3pt}',
         '\\begin{tabular}{lrrrrrrrrr}', '\\toprule',
         '& & \\multicolumn{4}{c}{Random $Z$} & \\multicolumn{4}{c}{Ancestor $Z$}\\\\',
         '\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}',
         'Graph & $|Z|$ & Mean & p99 & Fastest baseline & Speedup & Mean & p99 & Fastest baseline & Speedup\\\\']
    for g in KEYS:
        L.append('\\midrule')
        for i, k in enumerate(ZSIZES):
            parts = []
            for zm in ZMODES:
                c = cells[(g, zm, k)]; b = c['fastest_baseline']
                parts.append(f'{tus(c["ctskel"]["mean_us"])} & {tus(c["ctskel"]["p99_us"])} & '
                             f'{tus(c[b]["mean_us"])} {abbr[b]} & {spd(c["speedup_vs_fastest"])}')
            head = f'\\multirow{{{len(ZSIZES)}}}{{*}}{{{TAB[g]}}}' if i == 0 else ''
            L.append(f'{head} & {k} & ' + ' & '.join(parts) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular}', '\\end{table*}']
    write('tab-dsep.tex', L)


def dat_dsep(e2):
    cols = ['bidirectional', 'bidirectional-anc', 'cifly', 'rp', 'geiger', 'dcollision', 'ctskel']
    names = ['bidir', 'ancbidir', 'cifly', 'rp', 'geiger', 'dcol', 'index']
    for g in KEYS:
        L = [header(f'results/e2-merged (bench-s1-{g}, rp-s1-{g}), {CLASSIC_RUN} (dsep-{g})'),
             f'% {TXT[g]}: mean query time in microseconds ({E2_SEEDS_TEXT}) per method.',
             'zmode k ' + ' '.join(names)]
        for zm in ZMODES:
            for k in ZSIZES:
                c = e2['cells'][(g, zm, k)]
                L.append(f'{zm} {k} ' + ' '.join(f'{c[m]["mean_us"]:.3f}' for m in cols))
        write(f'dsep-{MAC[g]}.dat', L)


def p99_ratio(c):
    """Smallest baseline p99 in the cell / index p99 (the baseline with the best mean need not
    have the best tail, so the comparison is against the best tail)."""
    b = min((c[name]['p99_us'], i, name) for i, (name, _) in enumerate(BASELINES))
    return b[0] / c['ctskel']['p99_us'], b[2]


def tab_p99(e2):
    cells = e2['cells']
    L = [header('results/e2-merged (bench-s1-*)'),
         '\\begin{table*}[t]\\centering\\small',
         '\\caption{Tail latency of the index: p99 query time in \\textmu s per cell, seed 1 '
         '(1{,}000 walk-pair queries per cell). (s): structural, non-causal.}',
         '\\label{tab:p99}', '\\setlength{\\tabcolsep}{4pt}',
         '\\begin{tabular}{l' + 'r' * (2 * len(ZSIZES)) + '}', '\\toprule',
         f'& \\multicolumn{{{len(ZSIZES)}}}{{c}}{{Random $Z$, $|Z|$}} & \\multicolumn{{{len(ZSIZES)}}}{{c}}{{Ancestor $Z$, $|Z|$}}\\\\',
         f'\\cmidrule(lr){{2-{1 + len(ZSIZES)}}}\\cmidrule(lr){{{2 + len(ZSIZES)}-{1 + 2 * len(ZSIZES)}}}',
         'Graph & ' + ' & '.join(str(k) for _ in ZMODES for k in ZSIZES) + '\\\\', '\\midrule']
    for g in KEYS:
        L.append(f'{TAB[g]} & ' + ' & '.join(sig(cells[(g, zm, k)]['ctskel']['p99_us']) for zm in ZMODES for k in ZSIZES) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular}', '\\end{table*}']
    write('tab-p99.tex', L)


def crit_configs(zms):
    """Fixed plot order: random, ancestor, then parents (p) or children (c) where they exist."""
    order = [zm for zm in ('random', 'ancestor') if zm in zms] + [zm for zm in zms if zm not in ('random', 'ancestor')]
    return [(zm, k) for zm in order for k in CRIT_K]


def dat_crit(e3):
    for kind, key, ctext, _, zms, src in CRITERIA:
        srcs = {'round2': 'results/round2 (causal-*), round2-extra (causal-ikraph-causal-dag), round2d (causal-twitter)',
                'round3-rules': 'results/round3 (rules-s[123]-*)',
                'round3-further': 'results/round3 (further-s[123]-*)'}[src]
        L = [header(srcs),
             f'% {ctext}: speedup over traversal (traversal mean / index mean) per configuration'
             + ('; surgery-first on' if kind in ('rule2', 'rule3') else '')
             + '; labels r/a/p/c = random/ancestor/parents/children Z, number = |Z|.',
             'pos label ' + ' '.join(MAC[g] for g in KEYS)]
        for pos, (zm, k) in enumerate(crit_configs(zms)):
            vals = [crit_value(e3, kind, g, zm, k) for g in KEYS]
            L.append(f'{pos} {zm[0]}{k} ' + ' '.join('nan' if v is None else f'{v:.4f}' for v in vals))
        write(f'crit-{key}.dat', L)


def load_crit_bars(e3):
    """{(kind, g, zmode, k): {baseline: mean_us or None, 'ctskel': mean_us, 'source': ...}} for every
    configuration of every criterion. criteria-all (every baseline) where the graph's job has
    finished, otherwise the older runs (Bidirectional = "traversal" and the index only)."""
    out = {}
    for kind, _, _, _, zms, _ in CRITERIA:
        for g in KEYS:
            p = RES / 'criteria-all' / f'{CRIT_ALL_JOB[kind]}-{g}.stdout'
            new = {}
            try:
                for r in rows(p):
                    if r['kind'] != kind:
                        continue
                    ms = meth(r)
                    check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{rel(p)}: mismatches')
                    rec = {b: ms[b]['mean_us'] for b in OLD_BASELINES}
                    rec['ctskel'] = ms['ctskel']['mean_us']; rec['source'] = rel(p)
                    new[(r['zmode'], r['k'])] = rec
            except SourceError:
                pass
            classic = {}
            pc = RES / CLASSIC_RUN / f'{CRIT_ALL_JOB[kind]}-{g}.stdout'
            for r in rows(pc):
                if r['kind'] != kind:
                    continue
                ms = meth(r)
                check(list(ms) == ['bidirectional', *CLASSIC], f'{rel(pc)}: methods {list(ms)}')
                check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{rel(pc)}: mismatches')
                classic[(r['zmode'], r['k'])] = {b: ms[b]['mean_us'] for b in CLASSIC}
            for zm in zms:
                for k in CRIT_K:
                    rec = new.get((zm, k))
                    if rec is None:     # older runs: Bidirectional and the index only
                        c = e3['cells'][(kind, g, zm, k)]
                        rec = {b: None for b, _ in BASELINES}
                        rec['bidirectional'] = c['traversal_mean_us']; rec['ctskel'] = c['index_mean_us']
                        rec['source'] = c['source']
                    check((zm, k) in classic, f'{rel(pc)}: {kind} {zm}-{k} missing')
                    rec.update(classic[(zm, k)])
                    out[(kind, g, zm, k)] = rec
    return out


def dat_crit_bars(cb):
    """gen/critbar-<C>-<z><k>.dat, z = r/a/p/c (random, ancestor, parents, children)."""
    for kind, key, ctext, _, zms, _ in CRITERIA:
        for zm in zms:
            for k in CRIT_K:
                recs = [cb[(kind, g, zm, k)] for g in KEYS]
                srcs = sorted({r['source'] for r in recs})
                L = [header(', '.join(srcs) + f'; Geiger, d-collision: results/{CLASSIC_RUN}'),
                     f'% {ctext}, {zm} sets of {k} vertices: mean query time in seconds per method; nan = not measured.',
                     'graph ' + ' '.join(COL[b] for b, _ in BASELINES) + ' index']
                for g, r in zip(KEYS, recs):
                    vals = [r[b] for b, _ in BASELINES] + [r['ctskel']]
                    L.append(f'{MAC[g]} ' + ' '.join('nan' if v is None else f'{v / 1e6:.6e}' for v in vals))
                write(f'critbar-{key}-{zm[0]}{k}.dat', L)
            # mean query time over the sizes CRIT_K (each size has the same number of queries)
            recs = [[cb[(kind, g, zm, k)] for k in CRIT_K] for g in KEYS]
            srcs = sorted({r['source'] for rs in recs for r in rs})
            L = [header(', '.join(srcs) + f'; Geiger, d-collision: results/{CLASSIC_RUN}'),
                 f'% {ctext}, {zm} sets: mean query time in seconds over sets of '
                 + ', '.join(map(str, CRIT_K)) + ' vertices (the same number of queries each), per method; nan = not measured.',
                 'graph ' + ' '.join(COL[b] for b, _ in BASELINES) + ' index']
            for g, rs in zip(KEYS, recs):
                vals = []
                for col in [b for b, _ in BASELINES] + ['ctskel']:
                    xs = [r[col] for r in rs]
                    vals.append(None if any(x is None for x in xs) else sum(xs) / len(xs))
                L.append(f'{MAC[g]} ' + ' '.join('nan' if v is None else f'{v / 1e6:.6e}' for v in vals))
            write(f'critbar-{key}-{zm[0]}m.dat', L)


def dat_ablation(abl):
    gs = [g for g in KEYS if any(c['graph'] == g for c in abl.values())]
    L = [header('results/decomp-order{,-r2,-r3} (dc-ab-*, dc-ba-*; git 88b068a, seed 1)'),
         '% Mean query time (microseconds) of the plain 2-hop index on the full trek DAG (zskel) and the core-tree index',
         '% (ctskel) on the same requests: medians over both method orders and three runs; ratio = geometric mean over',
         '% the two orders of the per-order median of zskel / ctskel (load_ablation).',
         'graph zmode k zskel_mean_us ctskel_mean_us ratio']
    for g in gs:
        for zm in ZMODES:
            for k in ZSIZES:
                c = abl[(g, zm, k)]
                L.append(f'{MAC[g]} {zm} {k} {c["zskel_mean_us"]:.3f} {c["ctskel_mean_us"]:.3f} {c["ratio"]:.4f}')
    write('abl-coretree.dat', L)


def crit_value(e3, kind, g, zm, k):
    c = e3['cells'].get((kind, g, zm, k))
    return c['speedup'] if c else None


def tab_criteria(e3):
    L = [header('results/round2 (causal-*), round2-extra (causal-ikraph-causal-dag), '
                'round2d (causal-twitter), round3 (rules-s[123]-*, further-s[123]-*)'),
         '\\begin{table*}[t]\\centering\\scriptsize',
         '\\caption{Causal criteria: speedup of the index over traversal per cell (traversal mean / index mean). '
         'Back-door and rule 1: seed 1; all others: means over seeds 1--3; 1{,}000 walk-pair queries per seed. '
         'Rules 2 and 3 with surgery-first exploration. First $Z$ group: parents of the treatment (back-door) '
         'or children (front-door). (s): structural, non-causal. Bold: below 1.}',
         '\\label{tab:criteria}', '\\setlength{\\tabcolsep}{3.5pt}',
         '\\begin{tabular}{ll' + 'r' * 9 + '}', '\\toprule',
         '& & \\multicolumn{3}{c}{Parents / children} & \\multicolumn{3}{c}{Random $Z$} & \\multicolumn{3}{c}{Ancestor $Z$}\\\\',
         '\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}',
         'Criterion & Graph & ' + ' & '.join(f'$|Z|={k}$' for _ in range(3) for k in CRIT_K) + '\\\\']
    for kind, _, _, tname, zms, _ in CRITERIA:
        L.append('\\midrule')
        first = zms[0] if zms[0] in ('parents', 'children') else None
        for i, g in enumerate(KEYS):
            vals = []
            for zm in (first, 'random', 'ancestor'):
                for k in CRIT_K:
                    v = crit_value(e3, kind, g, zm, k) if zm else None
                    vals.append(spd(v) if v is not None else '--')
            head = f'\\multirow{{{len(KEYS)}}}{{*}}{{{tname}}}' if i == 0 else ''
            L.append(f'{head} & {TAB[g]} & ' + ' & '.join(vals) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular}', '\\end{table*}']
    write('tab-criteria.tex', L)


def tab_surgery(sf, clo):
    """Exp-5 ablations: mean query time of the index without the technique / with it (median of three
    runs per method order, geometric mean of the two orders; load_surgery, load_closure).
    tab-surgery.tex (Table 7): surgery-first on rule 2 for every configuration. The up-closure
    ablation (load_closure) is kept in paper-data.json but not reported.
    tab-surgery-full.tex (appendix): surgery-first on rules 2 and 3 for every configuration.
    No bold for values below 1 here: ratios within 0.93-1.05 are run-to-run noise (the |Z| <= 2
    up-closure cells, where both variants run the same code, span this range).
    Table 7 shades every speedup with \\heat (figures/heat_macros.tex, from the figure session):
    no tint at 1 or below, the deepest tint at the largest speedup of the table, on a log scale."""
    cols = [(zm, k) for zm in ZMODES for k in CRIT_K]
    gains = {(g, zm, k): sf[('rule2', g, zm, k)]['gain'] for g in KEYS for zm, k in cols}
    hi = max(gains.values())
    L = [header('results/surgery-order{,-r2,-r3} (git 88b068a, seed 1)'),
         '\\input{figures/heat_macros}',
         '\\begin{table}[t]\\centering\\small',
         '\\caption{Speedup from surgery-first exploration on rule 2.}',
         '\\label{tab:surgery}', '\\setlength{\\tabcolsep}{3pt}',
         colwidth_tabular('l' + 'r' * len(cols)), '\\toprule',
         '& \\multicolumn{3}{c}{Random} & \\multicolumn{3}{c}{Ancestor}\\\\',
         '\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}',
         'Graph & ' + ' & '.join(f'$k={k}$' for _, k in cols) + '\\\\', '\\midrule']
    for g in KEYS:
        L.append(f'{TAB[g]} & ' + ' & '.join(f'\\heat{{{gains[(g, zm, k)]:.6g}}}{{{sig(gains[(g, zm, k)])}}}{{1}}{{{hi:.6g}}}'
                                            for zm, k in cols) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-surgery.tex', L)

    # Appendix: rules 2 and 3 in the layout of Table 6, one block per rule, shaded on one scale (the largest speedup
    # of both rules).
    hi = max(sf[(kind, g, zm, k)]['gain'] for kind in ('rule2', 'rule3') for g in KEYS for zm, k in cols)
    L = [header('results/surgery-order{,-r2,-r3} (git 88b068a, seed 1)'),
         '\\input{figures/heat_macros}',
         '\\begin{table}[tp]\\centering\\small',
         '\\caption{Speedup from surgery-first exploration on rules 2 and 3.}',
         '\\label{tab:surgery-full}', '\\setlength{\\tabcolsep}{3pt}',
         colwidth_tabular('l' + 'r' * len(cols)), '\\toprule',
         '& \\multicolumn{3}{c}{Random} & \\multicolumn{3}{c}{Ancestor}\\\\',
         '\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}',
         'Graph & ' + ' & '.join(f'$k={k}$' for _, k in cols) + '\\\\']
    for kind, tname in (('rule2', 'Rule 2'), ('rule3', 'Rule 3')):
        L.append('\\midrule')
        L.append(f'\\multicolumn{{{len(cols) + 1}}}{{@{{}}l}}{{\\emph{{{tname}}}}}\\\\')
        for g in KEYS:
            L.append(f'{TAB[g]} & ' + ' & '.join(
                f'\\heat{{{sf[(kind, g, zm, k)]["gain"]:.6g}}}{{{sig(sf[(kind, g, zm, k)]["gain"])}}}{{1}}{{{hi:.6g}}}'
                for zm, k in cols) + '\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-surgery-full.tex', L)


REPS = ('', '-r2', '-r3')   # three repetitions of each Exp-5 ablation run (directory suffixes)


def load_closure():
    """Up-closure ablation (Exp-5): mean query time of the core-tree index without up-closures
    ("ctskel-noclo") divided by that with them ("ctskel"). Both variants share one index and answer the
    same requests in one process (git 88b068a, seed 1: the requests of E2 seed 1). The variant timed
    second gains from the caches warmed by the first: at |Z| <= 2, where no closure is built and both run
    the same code, the second was up to 1.5x faster. Every graph therefore runs in both orders,
    closure-ablation{REP}/clo-G (ctskel first) and closure-ablation-rev{REP}/clo-G (ctskel-noclo first),
    three times each (REPS). Single queries occasionally stall for ~17 ms, which doubles a mean of a few
    microseconds, so the ratio of an order is the median of its three runs, and the reported ratio is the
    geometric mean of the two orders. Checks: no mismatches; the counters of both variants equal those
    of ablation-final default-G and noclosure-G, which answer the same requests."""
    A = RES / 'ablation-final'
    out = {}
    for g in KEYS:
        de = {(r['zmode'], r['zsize']): meth(r)['ctskel'] for r in rows(A / f'default-{g}.stdout')}
        nc = {(r['zmode'], r['zsize']): meth(r)['ctskel'] for r in rows(A / f'noclosure-{g}.stdout')}
        ratios = {}   # (order, zmode, k) -> [ratio per repetition]
        closures = {}
        for order, base, first in (('with_first', 'closure-ablation', 'ctskel'), ('without_first', 'closure-ablation-rev', 'ctskel-noclo')):
            for rep in REPS:
                d = base + rep
                for r in rows(RES / d / f'clo-{g}.stdout'):
                    names = [m['method'] for m in r['methods']]
                    check(names[0] == first and set(names) == {'ctskel', 'ctskel-noclo'}, f'{d}/clo-{g}: unexpected method order {names}')
                    ms = meth(r); key = (r['zmode'], r['zsize'])
                    c, n = ms['ctskel'], ms['ctskel-noclo']
                    check(c['mismatches_vs_first'] == 0 and n['mismatches_vs_first'] == 0, f'{d}/clo-{g} {key}: mismatches')
                    for f in COUNTERS + ('separated',):
                        check(c[f] == de[key][f], f'{d}/clo-{g} {key}: ctskel {f} differs from ablation-final default')
                        check(n[f] == nc[key][f], f'{d}/clo-{g} {key}: ctskel-noclo {f} differs from ablation-final noclosure')
                    ratios.setdefault((order,) + key, []).append(n['mean_us'] / c['mean_us'])
                    closures[key] = c['closures'] / QUERIES_PER_RUN
        for zm in ZMODES:
            for k in ZSIZES:
                rec = {'graph': g, 'zmode': zm, 'k': k}
                for order in ('with_first', 'without_first'):
                    v = ratios[(order, zm, k)]
                    check(len(v) == len(REPS), f'{g} {zm}-{k} {order}: {len(v)} repetitions')
                    rec[f'ratios_{order}'] = v
                    rec[f'ratio_{order}'] = statistics.median(v)
                rec['ratio'] = math.sqrt(rec['ratio_with_first'] * rec['ratio_without_first'])
                rec['closures_per_query'] = closures[(zm, k)]
                rec['source'] = rel(RES / 'closure-ablation{,-rev}{,-r2,-r3}' / f'clo-{g}.stdout')
                out[(g, zm, k)] = rec
    return out


def load_index_stats(e1):
    """Shortcut-arc statistics of TrekIndex (appendix, Exp-3): index-stats/stats.jsonl, one build
    per graph with the counters of Algorithm 1 (git 10b50d8). Checked against the builds of Table 5."""
    p = RES / 'index-stats' / 'stats.jsonl'
    check(p.exists() and (p.parent / 'DONE').exists(), f'missing source {rel(p)} (or DONE)')
    out = {}
    for l in p.read_text().splitlines():
        o = json.loads(l); g = pathlib.Path(o['graph']).stem
        if g in DROPPED:
            continue
        check(g in KEYS, f'{rel(p)}: unknown graph {g}')
        c = e1[g]
        check((o['width'], o['core_n'], o['label_entries'], o['core_shortcuts']) ==
              (c['width'], c['core_n'], c['label_entries'], c['shortcuts']), f'{g}: index-stats differ from Table 5')
        out[g] = o
    check(set(out) == set(KEYS), f'{rel(p)}: graphs {sorted(out)}')
    USED.add(rel(p))
    return out


def load_surgery():
    """Surgery-first ablation (Exp-5): mean query time of the core-tree index without surgery-first
    ("ctskel") divided by that with it ("ctskel-sf"), rules 2 and 3. As for the up-closures, both
    variants share one index and answer the same requests in one process, so every graph runs in both
    orders, surgery-order{REP}/sf-ab-G (off first) and sf-ba-G (on first), three times each (REPS),
    git 88b068a, seed 1; the gain is the geometric mean over the two orders of the median of three runs."""
    out = {}
    for g in KEYS:
        runs = {}
        for rep in REPS:
            D = RES / ('surgery-order' + rep)
            for order, job, first in (('off_first', 'sf-ab', 'ctskel'), ('on_first', 'sf-ba', 'ctskel-sf')):
                for r in rows(D / f'{job}-{g}.stdout'):
                    ms = meth(r)
                    check(list(ms)[0] == first, f'{D.name}/{job}-{g}: unexpected method order {list(ms)}')
                    check(all(m['mismatches_vs_first'] == 0 for m in ms.values()), f'{D.name}/{job}-{g}: mismatches')
                    runs.setdefault((r['kind'], r['zmode'], r['k']), {}).setdefault(order, []).append(
                        ms['ctskel']['mean_us'] / ms['ctskel-sf']['mean_us'])
        for (kind, zm, k), v in runs.items():
            check(set(v) == {'off_first', 'on_first'} and all(len(x) == len(REPS) for x in v.values()),
                  f'{g} {kind} {zm}-{k}: missing order or repetition')
            a, b = statistics.median(v['off_first']), statistics.median(v['on_first'])
            out[(kind, g, zm, k)] = {'kind': kind, 'graph': g, 'zmode': zm, 'k': k,
                                     'gains_off_first': v['off_first'], 'gains_on_first': v['on_first'],
                                     'gain_off_first': a, 'gain_on_first': b, 'gain': math.sqrt(a * b),
                                     'source': rel(RES / 'surgery-order{,-r2,-r3}' / f'sf-{{ab,ba}}-{g}.stdout')}
    for kind in ('rule2', 'rule3'):
        for g in KEYS:
            for zm in ZMODES:
                for k in CRIT_K:
                    check((kind, g, zm, k) in out, f'surgery-order: {kind} {g} {zm}-{k} missing')
    return out


def tab_apps(apps):
    """Exp-7 table: screening (CauseNet, ASER), the replay with observed covariates, and the PC workload, with the
    speedup over the fastest of the six baselines of the row. CEG screening and the replay with all parents are
    left out; model testing is replaced by the PC workload. The fastest method of each
    row is bold on the band of TrekIndex's column in the appendix tables (\\trekcell)."""
    cols = [b for b, _ in BASELINES]
    L = [header(f'results/apps-all (screen-*, med-causenet-observed), {PC_RUN} (pc-*); Geiger, d-collision: {CLASSIC_RUN}'),
         '\\input{figures/heat_macros}',
         '\\begin{table}[t]\\centering\\scriptsize',
         '\\caption{Total query time of the application workloads.}',
         '\\label{tab:apps}', '\\setlength{\\tabcolsep}{1.6pt}',
         colwidth_tabular('l' + 'r' * (len(cols) + 3)), '\\toprule',
         'Graph & Queries & ' + ' & '.join(FULL_ABBR[c] for c in cols) + ' & \\Method{} & Speedup\\\\']
    by = {(a['workload'], a['graph']): a for a in apps}
    labels = {'screen': 'Covariate screening', 'replay-observed': 'Literature replay', 'pc': 'Conditional-independence tests of PC'}
    for w, graphs in APP_TABLE + (('pc', tuple(KEYS)),):
        L.append('\\midrule')
        L.append(f'\\multicolumn{{{len(cols) + 4}}}{{@{{}}l}}{{\\emph{{{labels[w]}}}}}\\\\')
        for g in graphs:
            a = by.get((w, g))
            if a is None or a['missing']:
                L.append(f'{TAB[g]} & ' + ' & '.join([RED] * (len(cols) + 3)) + '\\\\'); continue
            t = a['baseline_query_s']
            done = {c: t[c] for c in cols if c in t and c not in a.get('pending', ())}
            best = min(list(done) + ['index'], key=lambda c: a['index_query_s'] if c == 'index' else done[c])
            mark = lambda c, s: bestcell(s) if c == best else s
            vals = [RED if c in a.get('pending', ()) else mark(c, tsec(t[c])) if c in t else '--' for c in cols]
            L.append(f'{TAB[g]} & {num(a["queries"])} & ' + ' & '.join(vals)
                     + f' & {mark("index", tsec(a["index_query_s"]))} & {spd(a["query_speedup"])}\\\\')
    L += ['\\bottomrule', '\\end{tabular*}', '\\end{table}']
    write('tab-apps.tex', L)


# ------------------------------------------------------------------- macros
def crit_speedup(recs):
    """Index speedup over the fastest of the six baselines for the summed time of recs (one cell, or the
    cells of one group)."""
    check(all(r[b] is not None for r in recs for b, _ in BASELINES), f'criteria: a baseline is missing ({recs[0]["source"]})')
    t, b = min((sum(r[b] for r in recs), b) for b, _ in BASELINES)
    return (t / sum(r['ctskel'] for r in recs), b)


# Appendix tables: the mean query time of every method behind Figures 5 and 6. A
# baseline's cell is shaded by its slowdown over TrekIndex in the same row (\\heat, figures/heat_macros.tex, log
# scale from 1x to the largest slowdown of the table); TrekIndex's times are bold.
def timecell(us, idx_us, hi):
    """A baseline's time in ms, on a tint (8-60%) that grows with its slowdown over TrekIndex on a log scale up
    to hi (the 95th percentile of the table's slowdowns). The tint is computed here, since pgfmath cannot take
    slowdowns above 16383."""
    t = 8 + 52 * max(0.0, min(1.0, math.log(us / idx_us) / math.log(hi)))
    return f'\\heatp{{{t:.1f}}}{{{sig(us / 1e3)}}}'


def trekcell(us):
    return f'\\trekcell{{{sig(us)}}}'   # TrekIndex's time in microseconds, on the band of its column


def bestcell(s):
    return f'\\trekcell{{{s}}}'   # the fastest time of a row in Tables 5 and 7, styled like TrekIndex's column above


def p95(xs):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(0.95 * len(xs)))]


FULL_ABBR = {'bidirectional': 'BB', 'bidirectional-anc': 'Anc-BB', 'cifly': 'CIfly', 'rp': 'RP',
             'geiger': 'Geiger', 'dcollision': 'd-collision'}


def load_answers():
    """The share of d-separated answers per configuration of Figure 5 (e2-merged, seed 1; every method gives the
    same answers) and the share of queries for which each criterion holds per criterion, graph, kind and size of
    set (criteria-all; every method gives the same answers), shown in the appendix time tables."""
    sep = {}
    for g in KEYS:
        p = RES / 'e2-merged' / f'bench-s1-{g}.stdout'
        for r in rows(p):
            ms = meth(r)
            check(len({m['separated'] for m in ms.values()}) == 1, f'{rel(p)}: methods disagree on the answers')
            sep[(g, r['zmode'], r['zsize'])] = ms['bidirectional']['separated'] / r['queries']
    hold = {}
    for kind, key, ctext, tname, zms, _ in CRITERIA:
        for g in KEYS:
            p = RES / 'criteria-all' / f'{CRIT_ALL_JOB[kind]}-{g}.stdout'
            for r in rows(p):
                if r['kind'] != kind:
                    continue
                ms = meth(r)
                check(len({m['holds'] for m in ms.values()}) == 1, f'{rel(p)}: methods disagree on the answers')
                hold[(kind, g, r['zmode'], r['k'])] = ms['ctskel']['holds'] / r['queries']
    return sep, hold


def longtable(src, spec, ncols, caption, label, head, rows, font='\\footnotesize', sep='3pt'):
    """A longtable of the appendix (one-column pages): it breaks across pages like text, so no
    float leaves the rest of a page empty, and it repeats its header on every page."""
    return [src, '\\begingroup' + font + f'\\setlength{{\\tabcolsep}}{{{sep}}}\\setlength{{\\LTpre}}{{6pt}}\\setlength{{\\LTpost}}{{6pt}}',
            f'\\begin{{longtable}}{{{spec}}}', f'\\caption{{{caption}}}\\label{{{label}}}\\\\', '\\toprule', *head, '\\midrule',
            '\\endfirsthead', f'\\multicolumn{{{ncols}}}{{@{{}}l}}{{\\emph{{Table~\\thetable{{}} (continued)}}}}\\\\', '\\toprule', *head,
            '\\midrule', '\\endhead', '\\bottomrule', '\\endlastfoot', *rows, '\\end{longtable}', '\\endgroup']


def tab_dsep_full(e2, sep):
    """Table: mean d-separation query time of every method and share of d-separated answers per cell (Figure 5)."""
    cells = e2['cells']; base = [b for b, _ in BASELINES]; nb = len(base)
    hi = p95([cells[(g, zm, k)][b]['mean_us'] / cells[(g, zm, k)]['ctskel']['mean_us']
              for g in KEYS for zm in ZMODES for k in ZSIZES for b in base])
    half = 'r' * (nb + 2)
    head = [f'& & \\multicolumn{{{nb + 2}}}{{c}}{{Random $Z$}} & \\multicolumn{{{nb + 2}}}{{c}}{{Ancestor $Z$}}\\\\',
            f'\\cmidrule(lr){{3-{4 + nb}}}\\cmidrule(lr){{{5 + nb}-{6 + 2 * nb}}}',
            f'& & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\% & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\%\\\\',
            'Graph & $|Z|$ & ' + ' & '.join(([FULL_ABBR[b] for b in base] + ['\\Method{}', 'Sep.']) * 2) + '\\\\']
    rows = []
    for gi, g in enumerate(KEYS):
        if gi:
            rows.append('\\midrule')
        for i, k in enumerate(ZSIZES):
            vals = []
            for zm in ZMODES:
                c = cells[(g, zm, k)]; ix = c['ctskel']['mean_us']
                vals += [timecell(c[b]['mean_us'], ix, hi) for b in base] + [trekcell(ix), f'{100 * sep[(g, zm, k)]:.0f}']
            rows.append(f'{TAB[g] if i == 0 else ""} & {k} & ' + ' & '.join(vals) + '\\\\')
    write('tab-dsep-full.tex', longtable(header('results/e2-merged (bench-s1-*, rp-s1-*)') + '\n\\input{figures/heat_macros}',
          '@{}lr' + half + '@{\\hspace{10pt}}' + half + '@{}', 2 * (nb + 2) + 2,
          'Mean d-separation query time of every method and share of d-separated answers, per configuration; a darker cell of a baseline is slower relative to \\Method{}.',
          'tab:dsep-full', head, rows))


def tab_crit_full(cb, hold):
    """Tables: mean query time of every method and share of queries for which the criterion holds, for random and
    ancestor sets of every size (one longtable), and for the parent sets of the back-door and the child sets of the
    front-door criterion (the data of Figure 6). +R: mean of the two method orders."""
    base = [b for b, _ in BASELINES]; nb = len(base)
    names = [FULL_ABBR[m] for m in base] + ['\\Method{}']
    src = header(f'results/criteria-all, {CLASSIC_RUN} (as Figure 6)') + '\n\\input{figures/heat_macros}'
    hi = p95([r[m] / r['ctskel'] for key, r in cb.items() if key[1] in KEYS for m in base if r.get(m) is not None])

    def row(key):
        r = cb[key]
        return (['--' if r.get(m) is None else timecell(r[m], r['ctskel'], hi) for m in base] + [trekcell(r['ctskel'])]
                + [f'{100 * hold[key]:.0f}'])
    half = 'r' * (nb + 2)
    head = [f'& & & \\multicolumn{{{nb + 2}}}{{c}}{{Random}} & \\multicolumn{{{nb + 2}}}{{c}}{{Ancestor}}\\\\',
            f'\\cmidrule(lr){{4-{5 + nb}}}\\cmidrule(lr){{{6 + nb}-{7 + 2 * nb}}}',
            f'& & & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\% & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\%\\\\',
            'Criterion & Graph & $k$ & ' + ' & '.join((names + ['Holds']) * 2) + '\\\\']
    rows = []
    for ci, (kind, key, ctext, tname, zms, _) in enumerate(CRITERIA):
        if ci:
            rows.append('\\midrule')
        for gi, g in enumerate(KEYS):
            if gi:
                rows.append('\\addlinespace[2pt]')
            for j, k in enumerate(CRIT_K):
                lead = (tname if gi == 0 and j == 0 else '') + ' & ' + (TAB[g] if j == 0 else '')
                rows.append(f'{lead} & {k} & ' + ' & '.join(row((kind, g, 'random', k)) + row((kind, g, 'ancestor', k))) + '\\\\')
    write('tab-crit-full.tex', longtable(src, '@{}llr' + half + '@{\\hspace{8pt}}' + half + '@{}', 2 * (nb + 2) + 3,
          'Mean query time of every method and share of queries for which the criterion holds, for random and ancestor sets of 1, 4, and 16 vertices; a darker cell of a baseline is slower relative to \\Method{}.',
          'tab:crit-full', head, rows, sep='2.4pt'))
    # Parent and child sets side by side, as random and ancestor sets above: the back-door criterion (parents)
    # on the left, the front-door criterion (children) on the right.
    pc = [(kind, tname, zm) for kind, key, ctext, tname, zms, _ in CRITERIA for zm in ('parents', 'children') if zm in zms]
    check([z for _, _, z in pc] == ['parents', 'children'], 'one criterion with parent sets, then one with child sets')
    head = [f'& & \\multicolumn{{{nb + 2}}}{{c}}{{{pc[0][1]}, parent sets}} & \\multicolumn{{{nb + 2}}}{{c}}{{{pc[1][1]}, child sets}}\\\\',
            f'\\cmidrule(lr){{3-{4 + nb}}}\\cmidrule(lr){{{5 + nb}-{6 + 2 * nb}}}',
            f'& & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\% & \\multicolumn{{{nb}}}{{c}}{{Baselines (ms)}} & \\textmu s & \\%\\\\',
            'Graph & $k$ & ' + ' & '.join((names + ['Holds']) * 2) + '\\\\']
    rows = []
    for gi, g in enumerate(KEYS):
        if gi:
            rows.append('\\addlinespace[2pt]')
        for j, k in enumerate(CRIT_K):
            rows.append(f'{TAB[g] if j == 0 else ""} & {k} & ' + ' & '.join(row((pc[0][0], g, 'parents', k)) + row((pc[1][0], g, 'children', k))) + '\\\\')
    write('tab-crit-full-pc.tex', longtable(src, '@{}lr' + half + '@{\\hspace{8pt}}' + half + '@{}', 2 * (nb + 2) + 2,
          'Mean query time of every method and share of queries for which the criterion holds, for the back-door criterion with parent sets and the front-door criterion with child sets of 1, 4, and 16 vertices; a darker cell of a baseline is slower relative to \\Method{}.',
          'tab:crit-full-pc', head, rows, sep='2.4pt'))


def dat_appendix(sep, hold, ab, sf):
    """Data of the appendix figures (figures/app_plots.tex), the same numbers as the tables above:
    gen/app-sep.dat, the share of d-separated answers (%) per graph and configuration of Figure 5;
    gen/app-idx-<K>.dat, the mean query time (microseconds) of BB, trek traversal, trek search and TrekIndex and
    the edge tests of trek search per query (Exp-4, every configuration);
    gen/app-sf-<rule>.dat, the speedup from surgery-first exploration on rules 2 and 3 (Exp-6);
    gen/app-hold-<C>.dat, the share of queries for which criterion C holds (%), one row per graph and one column
    per kind and size of set (r1 = random sets of one vertex, a16 = ancestor sets of 16, p4 = parent sets of 4, ...)."""
    cols = ' '.join(MAC[g] for g in KEYS)
    L = [header('results/e2-merged (bench-s1-*)'), '% Share of d-separated answers (%).', f'zmode k {cols}']
    for zm in ZMODES:
        for k in ZSIZES:
            L.append(f'{zm} {k} ' + ' '.join(f'{100 * sep[(g, zm, k)]:.2f}' for g in KEYS))
    write('app-sep.dat', L)
    for g in KEYS:
        L = [header(f'results/{ABLA_RUN} (abl-{g}; code-ablB, seed 1)'),
             f'% {TXT[g]}: mean query time in microseconds, and edge tests of trek search per query.',
             'zmode k bidir trav edge index tests']
        for zm in ZMODES:
            for k in ZSIZES:
                c = ab[(g, zm, k)]
                L.append(f'{zm} {k} {c["bi_us"]:.3f} {c["trav_us"]:.3f} {c["edge_us"]:.3f} {c["idx_us"]:.3f} {c["tests"]:.3f}')
        write(f'app-idx-{MAC[g]}.dat', L)
    for kind, key in (('rule2', 'RuleTwo'), ('rule3', 'RuleThree')):
        L = [header('results/surgery-order{,-r2,-r3} (git 88b068a, seed 1)'),
             f'% {kind}: speedup from surgery-first exploration.', f'zmode k {cols}']
        for zm in ZMODES:
            for k in CRIT_K:
                L.append(f'{zm} {k} ' + ' '.join(f'{sf[(kind, g, zm, k)]["gain"]:.4f}' for g in KEYS))
        write(f'app-sf-{key}.dat', L)
    for kind, key, ctext, tname, zms, _ in CRITERIA:
        L = [header('results/criteria-all'), f'% {ctext}: share of queries for which the criterion holds (%).',
             'graph ' + ' '.join(f'{zm[0]}{k}' for zm in zms for k in CRIT_K)]
        for g in KEYS:
            L.append(f'{MAC[g]} ' + ' '.join(f'{100 * hold[(kind, g, zm, k)]:.2f}' for zm in zms for k in CRIT_K))
        write(f'app-hold-{key}.dat', L)


def macros(e1, e2, e3, apps, abl, et, sf, clo, cb, ist):
    M = []   # (name, latex value, raw value, note)

    def add(name, val, raw=None, note=''):
        check(re.fullmatch(r'[A-Za-z]+', name), f'bad macro name {name}')
        M.append((name, val, raw, note))

    # E1
    gmax = max(KEYS, key=lambda g: e1[g]['build_s'])
    add('BuildMaxSec', sig(e1[gmax]['build_s']), e1[gmax]['build_s'], 'max core-tree build time over the graphs, seconds')
    add('BuildMaxMin', f'{e1[gmax]["build_s"] / 60:.1f}', e1[gmax]['build_s'] / 60, 'same, minutes')
    add('BuildMaxGraph', TXT[gmax], gmax)
    gv = max(KEYS, key=lambda g: e1[g]['n'])
    add('MaxVerticesMillion', f'{e1[gv]["n"] / 1e6:.1f}', e1[gv]['n'], f'largest |V| over the graphs ({TXT[gv]}), millions')
    cn = e1['causenet-full-dag']
    add('CauseNetPlainEntries', words(cn['plain']['label_entries']), cn['plain']['label_entries'], '2-hop on the full trek DAG')
    add('CauseNetCoreEntries', words(cn['label_entries']), cn['label_entries'], 'core-tree index')
    add('CauseNetEntryRatio', sig(cn['plain']['label_entries'] / cn['label_entries']), cn['plain']['label_entries'] / cn['label_entries'])
    add('CauseNetPlainBuildSec', sig(cn['plain']['build_s']), cn['plain']['build_s'], 'mean of three builds')
    add('CauseNetCoreBuildSec', sig(cn['build_s']), cn['build_s'], 'mean of three builds')
    cf = e1['dcol-cnfdag']
    add('CNFDAGPlainEntries', words(cf['plain']['label_entries']), cf['plain']['label_entries'], '2-hop on the full trek DAG')
    add('CNFDAGCoreEntries', words(cf['label_entries']), cf['label_entries'], 'core-tree index')
    add('CNFDAGEntryRatio', sig(cf['plain']['label_entries'] / cf['label_entries']), cf['plain']['label_entries'] / cf['label_entries'])
    # index size (labels + peripheral data) of the plain 2-hop index / TrekIndex on the graphs whose
    # periphery holds more than 95% of the vertices (the text names this threshold)
    big = [g for g in KEYS if e1[g]['peeled_frac'] > 0.95]
    check(big == ['dcol-cnfdag', 'causenet-full-dag'], f'Exp-1/Exp-4 text: graphs with >95% periphery are now {big}')
    rs = sorted((e1[g]['plain']['label_bytes'] / e1[g]['index_bytes'], g) for g in big)
    add('PeriphSizeRatioMin', sig(rs[0][0]), rs[0][0], f'plain index bytes / TrekIndex bytes, {TXT[rs[0][1]]}')
    add('PeriphSizeRatioMax', sig(rs[-1][0]), rs[-1][0], TXT[rs[-1][1]])

    # E2
    cells = list(e2['cells'].values())
    add('DsepCells', num(len(cells)), len(cells))
    lo = min(cells, key=lambda c: c['speedup_vs_fastest']); hi = max(cells, key=lambda c: c['speedup_vs_fastest'])
    add('DsepMinAll', sig(lo['speedup_vs_fastest']), lo['speedup_vs_fastest'], 'vs fastest baseline')
    add('DsepMinAllGraph', TXT[lo['graph']], lo['graph'])
    add('DsepMinAllCell', f'{lo["zmode"]}, $|Z|={lo["k"]}$', f'{lo["zmode"]}-{lo["k"]}')
    big = [c for c in cells if e1[c['graph']]['n'] >= 10**6]
    lb = min(big, key=lambda c: c['speedup_vs_fastest'])
    add('DsepMinMillion', sig(lb['speedup_vs_fastest']), lb['speedup_vs_fastest'], f'graphs with >= 1M vertices; cell {lb["zmode"]}-{lb["k"]}')
    add('DsepMinMillionGraph', TXT[lb['graph']], lb['graph'])
    add('DsepMaxAll', sig(hi['speedup_vs_fastest']), hi['speedup_vs_fastest'], f'cell {hi["zmode"]}-{hi["k"]}')
    add('DsepMaxAllGraph', TXT[hi['graph']], hi['graph'])
    hc = max((c for c in cells if KIND[c['graph']] == 'causal'), key=lambda c: c['speedup_vs_fastest'])
    add('DsepMaxCausal', sig(hc['speedup_vs_fastest']), hc['speedup_vs_fastest'], f'five causal graphs; cell {hc["zmode"]}-{hc["k"]}')
    add('DsepMaxCausalGraph', TXT[hc['graph']], hc['graph'])
    for g in KEYS:
        gc = [c for c in cells if c['graph'] == g]
        a = min(gc, key=lambda c: c['speedup_vs_fastest']); b = max(gc, key=lambda c: c['speedup_vs_fastest'])
        add(f'DsepMin{MAC[g]}', sig(a['speedup_vs_fastest']), a['speedup_vs_fastest'], f'{a["zmode"]}-{a["k"]} vs {a["fastest_baseline"]}')
        add(f'DsepMax{MAC[g]}', sig(b['speedup_vs_fastest']), b['speedup_vs_fastest'], f'{b["zmode"]}-{b["k"]} vs {b["fastest_baseline"]}')

    # E2 tail latency: index p99 (mean of the per-seed p99) and best baseline p99 / index p99
    hp = max(cells, key=lambda c: c['ctskel']['p99_us'])
    add('PNineNineMaxIndex', tus(hp['ctskel']['p99_us']), hp['ctskel']['p99_us'], 'largest index p99 over the 96 cells')
    add('PNineNineMaxIndexWhich', f'{TXT[hp["graph"]]}, {hp["zmode"]}-{hp["k"]}')
    lr = min(cells, key=lambda c: p99_ratio(c)[0]); r, b = p99_ratio(lr)
    add('PNineNineMinRatio', sig(r), r, 'min over cells of (smallest baseline p99 / index p99)')
    add('PNineNineMinRatioWhich', f'{TXT[lr["graph"]]}, {lr["zmode"]}-{lr["k"]}, vs {BASE_TEXT[b]}')

    # E3: the criteria against the fastest of the four baselines (criteria-all). Exp-2 reports the
    # mean query time over the set sizes CRIT_K for each criterion, graph, and kind of set (a group);
    # the CritCell* macros keep the per-size cells, which the text also reports.
    speedup = crit_speedup
    groups, cfgs = [], []
    for kind, key, ctext, _, zms, _ in CRITERIA:
        gs = []
        for g in KEYS:
            for zm in zms:
                recs = [cb[(kind, g, zm, k)] for k in CRIT_K]
                for k, r in zip(CRIT_K, recs):
                    sp, b = speedup([r])
                    cfgs.append((sp, ctext, g, zm, k, b))
                sp, b = speedup(recs)
                gs.append((sp, ctext, g, zm, b))
        groups += gs
        a, z = min(gs), max(gs)
        add(f'CritMin{key}', sig(a[0]), a[0], f'{TXT[a[2]]}, {a[3]} sets, mean over sizes, vs {BASE_TEXT[a[4]]}')
        add(f'CritMax{key}', sig(z[0]), z[0], f'{TXT[z[2]]}, {z[3]} sets, mean over sizes, vs {BASE_TEXT[z[4]]}')
    add('CritGroups', num(len(groups)), len(groups), 'criterion x graph x kind of set')
    a, z = min(groups), max(groups)
    add('CritMinAll', sig(a[0]), a[0], f'mean over sizes, vs {BASE_TEXT[a[4]]} (fastest of the six baselines)')
    add('CritMinAllWhich', f'{a[1]}, {TXT[a[2]]}, {a[3]} sets')
    add('CritMaxAll', sig(z[0]), z[0], f'{z[1]}, {TXT[z[2]]}, {z[3]} sets, mean over sizes, vs {BASE_TEXT[z[4]]}')
    add('CritGroupsSlower', num(sum(1 for c in groups if c[0] < 1)), sum(1 for c in groups if c[0] < 1))
    cz = max(c for c in groups if KIND[c[2]] == 'causal')
    add('CritMaxCausal', sig(cz[0]), cz[0], f'largest group speedup on the real causal graphs ({cz[1]}, {TXT[cz[2]]}, {cz[3]} sets), vs {BASE_TEXT[cz[4]]}')
    add('CritMaxAllGraph', TXT[z[2]], None)
    # Figure 6 shows random sets only: the Exp-2 text discusses these groups (one per criterion and graph)
    rg = [c for c in groups if c[3] == 'random']
    add('CritRandGroups', num(len(rg)), len(rg), 'criterion x graph, random sets (Figure 6)')
    a2, z2 = min(rg), max(rg)
    add('CritRandMin', sig(a2[0]), a2[0], f'random sets, mean over sizes, vs {BASE_TEXT[a2[4]]}')
    add('CritRandMinWhich', f'{a2[1]}, {TXT[a2[2]]}')
    add('CritRandMax', sig(z2[0]), z2[0], f'random sets, mean over sizes, vs {BASE_TEXT[z2[4]]}')
    add('CritRandMaxWhich', f'{z2[1]}, {TXT[z2[2]]}')
    add('CritRandSlower', num(sum(1 for c in rg if c[0] < 1)), sum(1 for c in rg if c[0] < 1), 'random groups where the index is slower')
    # finer analysis of Figure 6: by graph, and the effect of the labels on the criteria with directed-path tests
    reach_kinds = {'backdoor', 'rule3', 'frontdoor', 'seqbackdoor', 'mediation'}   # criteria with Reach tests (Table 2)
    big = [g for g in KEYS if e1[g]['peeled_frac'] > 0.95]
    pm = min(c[0] for c in rg if c[2] in big)
    add('CritRandPeriphMin', sig(pm), pm, f'random groups on the graphs with >95% periphery ({", ".join(TXT[g] for g in big)})')
    # graph of the smallest random-group speedup of each criterion (Exp-2 text)
    lowg = {min((c for c in rg if c[1] == ctext), key=lambda c: c[0])[2] for kind, key, ctext, _, zms, _ in CRITERIA}
    add('CritRandMinGraphs', ' and '.join(TXT[g] for g in KEYS if g in lowg), None, 'graphs that give the smallest random-group speedup of some criterion')
    add('CritRandMinGraphCount', num(len(lowg)), len(lowg))
    ctext_of = {kind: ctext for kind, key, ctext, *_ in CRITERIA}
    rr = [c for c in rg if any(ctext_of[k] == c[1] for k in reach_kinds)]
    add('CritRandReachGroups', num(len(rr)), len(rr), 'random groups of the criteria with directed-path tests')
    from collections import Counter
    for b_, cnt in Counter(c[4] for c in rg).items():
        add(f'CritRandFastest{b_.replace("-", "").title()}', num(cnt), cnt, f'random groups whose fastest baseline is {BASE_TEXT[b_]}')
    add('CritCells', num(len(cfgs)), len(cfgs), 'per set size')
    slow = sorted(c for c in cfgs if c[0] < 1)
    add('CritCellsSlower', num(len(slow)), len(slow), 'per-size cells where the index is slower than the fastest of the six baselines')
    add('CritCellsSlowerGraphs', ' and '.join(TXT[g] for g in KEYS if any(c[2] == g for c in slow)), None)
    c = min(cfgs)
    add('CritCellMin', sig(c[0]), c[0], f'vs {BASE_TEXT[c[5]]}')
    add('CritCellMinWhich', f'{c[1]}, {TXT[c[2]]}, {c[3]}-{c[4]}')
    add('CritCellSlowdownMax', sig(1 / c[0]), 1 / c[0], 'largest slowdown of the index, = 1 / CritCellMin')

    # surgery-first
    for var, name in (('-off', 'Off'), ('', 'On')):
        cs = [(e3['cells'][(kind + var, g, zm, k)]['speedup'], ctext, g, zm, k)
              for kind, ctext in (('rule2', 'rule 2'), ('rule3', 'rule 3')) for g in KEYS for zm in ZMODES for k in CRIT_K]
        a = min(cs)
        add(f'SFMin{name}', sig(a[0]), a[0])
        add(f'SFMin{name}Which', f'{a[1]}, {TXT[a[2]]}, {a[3]}-{a[4]}')

    # surgery-first gain per cell: index time off / on, geometric mean of both method orders
    gains = [(sf[(kind, g, zm, k)]['gain'], ctext, g, zm, k)
             for kind, ctext in (('rule2', 'rule 2'), ('rule3', 'rule 3')) for g in KEYS for zm in ZMODES for k in CRIT_K]
    a, b = min(gains), max(gains)
    add('SFGainMin', sig(a[0]), a[0], 'index time off / on, over the rule 2/3 cells')
    add('SFGainMinWhich', f'{a[1]}, {TXT[a[2]]}, {a[3]}-{a[4]}')
    add('SFGainMax', sig(b[0]), b[0])
    add('SFGainMaxWhich', f'{b[1]}, {TXT[b[2]]}, {b[3]}-{b[4]}')
    add('SFCells', num(len(gains)), len(gains))
    # Exp-6 text on rule 2, the rule of Table 7: every configuration gains, the gain grows with the size of the set,
    # and ancestor sets gain more than random sets of the same size
    g2 = {(g, zm, k): sf[('rule2', g, zm, k)]['gain'] for g in KEYS for zm in ZMODES for k in CRIT_K}
    check(min(g2.values()) > 1, 'Exp-6 text: surgery-first speeds up every rule-2 configuration')
    check(all(g2[(g, zm, a)] < g2[(g, zm, b)] for g in KEYS for zm in ZMODES for a, b in zip(CRIT_K, CRIT_K[1:])),
          'Exp-6 text: the rule-2 gain grows with the size of the set')
    check(all(g2[(g, 'ancestor', k)] > g2[(g, 'random', k)] for g in KEYS for k in CRIT_K),
          'Exp-6 text: ancestor sets gain more than random sets of the same size on rule 2')
    add('SFTwoMax', sig(max(g2.values())), max(g2.values()), 'surgery-first gain on rule 2 (Table 7), largest')

    # E1 facts used in the text
    add('IndexMinSize', nbytes(min(e1[g]['index_bytes'] for g in KEYS)), min(e1[g]['index_bytes'] for g in KEYS), 'labels plus tree part')
    add('IndexMaxSize', nbytes(max(e1[g]['index_bytes'] for g in KEYS)), max(e1[g]['index_bytes'] for g in KEYS))
    add('PeakRSSMax', nbytes(max(e1[g]['peak_rss_bytes'] for g in KEYS)), max(e1[g]['peak_rss_bytes'] for g in KEYS))

    # Exp-6: edge tests by class (edge-tests, seed 1)
    def pct(x):
        return f'{100 * x:.0f}'
    for g in ('causenet-full-dag', 'dcol-cnfdag'):
        lo = min((et[(g, zm, k)]['share']['proj'], zm, k) for zm in ZMODES for k in ZSIZES)
        add(f'ETProjMin{MAC[g]}', pct(lo[0]), lo[0], f'% of edge tests decided by the projections alone, minimum over all cells ({lo[1]}-{lo[2]})')
    for k, name in ((1, 'One'), (32, 'ThirtyTwo')):
        v = et[('syn-ba-2000k', 'ancestor', k)]['share']['refuted']
        add(f'ETRefutedBAtwoAnc{name}', pct(v), v, f'% of edge tests refuted by the labels, BA 2M, ancestor-{k}')
    rs = max((et[(g, 'random', k)]['share']['searched'], g, k) for g in KEYS for k in ZSIZES)
    add('ETSearchRandMax', pct(rs[0]), rs[0], f'% of edge tests that search, max over graphs, {TXT[rs[1]]}, random-{rs[2]}')
    spq = max((r['searches_per_query'], key) for key, r in et.items() if key[0] in KEYS)
    add('ETSearchPQMax', f'{spq[0]:.2f}', spq[0], f'searches per query, {TXT[spq[1][0]]}, {spq[1][1]}-{spq[1][2]}')
    add('ETSearchPQMaxWhich', f'{TXT[spq[1][0]]}, {spq[1][1]}, $|Z|={spq[1][2]}$')   # as DsepMinAllCell

    # E6: core-tree vs plain 2-hop query time (ratio = plain / core-tree)
    ac = list(abl.values())
    a = min(ac, key=lambda c: c['ratio']); b = max(ac, key=lambda c: c['ratio'])
    add('AblCTMinRatio', sig(a['ratio']), a['ratio'], 'plain 2-hop mean / core-tree mean')
    add('AblCTMinRatioWhich', f'{TXT[a["graph"]]}, {a["zmode"]}-{a["k"]}')
    add('AblCTMaxRatio', sig(b['ratio']), b['ratio'])
    add('AblCTMaxRatioWhich', f'{TXT[b["graph"]]}, {b["zmode"]}-{b["k"]}')
    for g in KEYS:
        gc = [c for c in ac if c['graph'] == g]
        if not gc:
            continue
        a = min(gc, key=lambda c: c['ratio']); b = max(gc, key=lambda c: c['ratio'])
        add(f'AblCTMin{MAC[g]}', sig(a['ratio']), a['ratio'], f'{a["zmode"]}-{a["k"]}')
        add(f'AblCTMax{MAC[g]}', sig(b['ratio']), b['ratio'], f'{b["zmode"]}-{b["k"]}')

    # Exp-3, appendix: shortcut arcs on the graphs with w = 3 (index-stats)
    w3 = [g for g in KEYS if ist[g]['width'] == 3]
    sc = {g: ist[g]['shortcuts_created'] for g in w3}
    lo = min(w3, key=lambda g: sc[g]); hi = max(w3, key=lambda g: sc[g])
    add('ShortcutMin', count(sc[lo]), sc[lo], f'shortcut arcs created, {TXT[lo]}')
    add('ShortcutMax', count(sc[hi]), sc[hi], f'shortcut arcs created, {TXT[hi]}')
    per = max(sc[g] / ist[g]['peeled'] for g in w3)
    add('ShortcutPerPeeledMax', sig(per), per, 'shortcut arcs created per peeled vertex, max over w = 3 graphs')
    pm = max(w3, key=lambda g: ist[g]['max_parallel_all'])
    add('ParallelMax', num(ist[pm]['max_parallel_all']), ist[pm]['max_parallel_all'], 'most parallel shortcut arcs between two states')
    add('ParallelMaxGraph', TXT[pm])
    add('WSGraphs', num(len(w3)), len(w3), 'graphs with w = 3')

    # E5
    by = {a['name']: a for a in apps}

    def app(name, macro, field, fmt):
        a = by.get(name)
        if a is None or a['missing'] or a.get(field) is None:
            add(macro, RED, None, f'missing: {name}')
        else:
            add(macro, fmt(a[field]), a[field], name)

    app('screen-causenet-full-dag', 'ScreenCauseNet', 'query_speedup', sig)
    app('screen-aser-concept-dag', 'ScreenASER', 'query_speedup', sig)
    app('med-causenet-observed', 'ReplayCauseNet', 'query_speedup', sig)
    app('med-causenet-all', 'ReplayCauseNetAll', 'query_speedup', sig)
    scr = [by.get(n) for n in ('screen-causenet-full-dag', 'screen-aser-concept-dag')]
    par = {(a['observed_frac'], a['cap']) for a in scr if a and not a['missing']}
    if len(par) == 1 and None not in next(iter(par)):
        of, cap = par.pop()
        add('ScreenObservedPct', f'{of * 100:.0f}', of, 'screening: probability that a vertex is observed (%)')
        add('ScreenCap', num(cap), cap, 'screening: most observed parents in Z')
    else:
        add('ScreenObservedPct', RED, None, 'screening parameters missing or differ'); add('ScreenCap', RED, None)
    # range of the speedups of Table 7 (one per row, over the fastest baseline of the row), for the abstract
    rowsp = []
    for w, graphs in APP_TABLE:
        for g in graphs:
            a = next((x for x in apps if x['workload'] == w and x['graph'] == g and not x['missing']), None)
            if a is not None:
                rowsp.append((a['query_speedup'], w, g))
    rowsp += [(a['query_speedup'], 'pc', a['graph']) for a in apps if a['workload'] == 'pc' and not a['missing']]
    if rowsp:
        lo, hi = min(rowsp), max(rowsp)
        add('AppMin', sig(lo[0]), lo[0], f'smallest Table 7 speedup ({lo[1]}, {TXT[lo[2]]})')
        add('AppMax', sig(hi[0]), hi[0], f'largest Table 7 speedup ({hi[1]}, {TXT[hi[2]]})')
    pcs = [a for a in apps if a['workload'] == 'pc']
    pcok = [a for a in pcs if not a['missing']]
    if pcs and len(pcok) == len(pcs):
        lo = min(pcok, key=lambda a: a['query_speedup']); hi = max(pcok, key=lambda a: a['query_speedup'])
        add('PCMin', sig(lo['query_speedup']), lo['query_speedup'], TXT[lo['graph']]); add('PCMinGraph', TXT[lo['graph']])
        add('PCMax', sig(hi['query_speedup']), hi['query_speedup'], TXT[hi['graph']]); add('PCMaxGraph', TXT[hi['graph']])
        add('PCGraphs', num(len(pcok)), len(pcok), 'graphs with an ancestral region of the requested size')
        for k, f in (('PCRegions', 'pc_regions'), ('PCRegionMin', 'pc_region_min'), ('PCRegionMax', 'pc_region_max'),
                     ('PCLevel', 'pc_level'), ('PCSample', 'queries')):
            v = {a[f] for a in pcok}; add(k, num(v.pop()) if len(v) == 1 else RED, None, f)
        tmin = min(a['pc_tests'] for a in pcok); tmax = max(a['pc_tests'] for a in pcok)
        add('PCTestsMinMillion', f'{tmin / 1e6:.1f}', tmin, 'tests of the full PC workload, smallest graph')
        add('PCTestsMaxMillion', f'{tmax / 1e6:.1f}', tmax, 'tests of the full PC workload, largest graph')
        for a in pcok:
            add(f'PC{MAC[a["graph"]]}', sig(a['query_speedup']), a['query_speedup'], f'{a["name"]} vs {a["fastest_baseline"]}')
        fbs = {a['fastest_baseline'] for a in pcok}
        fbn = FULL_ABBR[fbs.pop()] if len(fbs) == 1 else None
        add('PCFastestBaseline', fbn if fbn else RED, None, 'the fastest baseline on every PC graph (its abbreviation in the Competitors list)')
        # traversals without an ancestral restriction (Bidirectional, CIfly, Geiger, d-collision) relative to TrekIndex
        plain = [min(a['baseline_query_s'][b] for b in ('bidirectional', 'cifly', *CLASSIC)) / a['index_query_s'] for a in pcok]
        add('PCPlainMin', sig(min(plain)), min(plain), 'min over graphs of min(BB, CIfly, Geiger, d-collision) / TrekIndex')
        add('PCPlainMax', sig(max(plain)), max(plain), 'max over graphs of min(BB, CIfly, Geiger, d-collision) / TrekIndex')
    else:
        for k in ('PCMin', 'PCMinGraph', 'PCMax', 'PCMaxGraph', 'PCGraphs', 'PCRegions', 'PCRegionMin', 'PCRegionMax',
                  'PCLevel', 'PCSample', 'PCTestsMinMillion', 'PCTestsMaxMillion', 'PCFastestBaseline', 'PCPlainMin', 'PCPlainMax'):
            add(k, RED, None, 'missing: pc-v2')
    mk = [a for a in apps if a['workload'] == 'markov' and not a['missing']]
    mq = {a['queries'] for a in mk}
    add('MarkovQueries', num(mq.pop()) if len(mq) == 1 else RED, None, 'queries per model-testing job')
    mb = {(a['buckets'], a['queries_per_bucket']) for a in mk}
    add('MarkovBuckets', num(next(iter(mb))[0]) if len(mb) == 1 else RED, None, 'buckets by |Pa(x)| and sampling mode per job')
    add('MarkovQueriesPerBucket', num(next(iter(mb))[1]) if len(mb) == 1 else RED, None)
    for a in mk:
        add(f'Markov{MAC[a["graph"]]}', sig(a['query_speedup']), a['query_speedup'], a['name'])
    lo = min(mk, key=lambda a: a['query_speedup']); hi = max(mk, key=lambda a: a['query_speedup'])
    add('MarkovMin', sig(lo['query_speedup']), lo['query_speedup'], TXT[lo['graph']])
    add('MarkovMinGraph', TXT[lo['graph']])
    add('MarkovMax', sig(hi['query_speedup']), hi['query_speedup'], TXT[hi['graph']])
    add('MarkovMaxGraph', TXT[hi['graph']])
    # covariate sets of the replay with all parents (statistics written with the requests)
    ws = REPLAY_STATS or RES / 'causenet-mediation.stats.json'
    w = json.loads(ws.read_text())['W_sizes']['all']; USED.add(rel(ws))
    add('ReplayAllWMedian', num(w['p50']), w['p50'], 'median |W|, replay with all parents as covariates')
    add('ReplayAllWMax', num(w['max']), w['max'])
    # the text names the fastest baseline of each workload
    fb = {a['name']: a['fastest_baseline'] for a in apps if not a['missing']}
    check(all(fb[n] == 'bidirectional' for n in fb if n.startswith('screen-')), 'Exp-7 text: screening is no longer fastest with Bidirectional')
    check(fb.get('med-causenet-observed', 'bidirectional-anc') == 'bidirectional-anc', 'Exp-7 text: replay (observed) is no longer fastest with Ancestral bidirectional')
    check(all(fb[a['name']] == 'cifly' for a in mk), 'Exp-7 text: model testing is no longer fastest with CIfly')

    # checks
    add('Mismatches', num(e2['mismatches'] + e3['mismatches']), e2['mismatches'] + e3['mismatches'])
    add('MethodRuns', num(e2['runs'] + e3['runs']), e2['runs'] + e3['runs'],
        'E2: ancestral bidirectional, CIfly, index vs bidirectional; E3: index vs traversal; 1,000 queries each')
    names = [m[0] for m in M]
    check(len(names) == len(set(names)), 'duplicate macro')
    return M


def appendix_macros(e2, cb, sf, ans):
    """Macros of the analysis in the appendix (Additional Results for Section 7). Every claim of
    that text is checked here: the shares of answers come from tab_answers, the margins from the mean times that
    Tables 10-13 show (the fastest of the six baselines over TrekIndex), and the surgery-first gains from Table 8."""
    M = []

    def add(name, val, raw=None, note=''):
        check(re.fullmatch(r'[A-Za-z]+', name), f'bad macro name {name}')
        M.append((name, val, raw, note))

    sep, hold = ans
    cells = list(e2['cells'].values())
    nb = sum(1 for c in cells if c['fastest_baseline'] == 'bidirectional')
    add('AppDsepBidirFastest', num(nb), nb, 'd-separation configurations whose fastest baseline is Bidirectional')
    # shares of d-separated answers: random Z falls from |Z| = 1 to 32 on every graph, ancestor Z rises
    lo_k, hi_k = ZSIZES[0], ZSIZES[-1]
    check(all(sep[(g, 'random', hi_k)] < sep[(g, 'random', lo_k)] for g in KEYS), 'appendix text: random-Z share of separated answers no longer falls on every graph')
    check(all(sep[(g, 'ancestor', hi_k)] > sep[(g, 'ancestor', lo_k)] for g in KEYS), 'appendix text: ancestor-Z share of separated answers no longer rises on every graph')
    gr = max(KEYS, key=lambda g: sep[(g, 'random', lo_k)] - sep[(g, 'random', hi_k)])
    ga = max(KEYS, key=lambda g: sep[(g, 'ancestor', hi_k)] - sep[(g, 'ancestor', lo_k)])
    pct = lambda v: f'{100 * v:.0f}'
    add('AppSepRandGraph', TXT[gr]); add('AppSepRandFirst', pct(sep[(gr, 'random', lo_k)])); add('AppSepRandLast', pct(sep[(gr, 'random', hi_k)]))
    add('AppSepAncGraph', TXT[ga]); add('AppSepAncFirst', pct(sep[(ga, 'ancestor', lo_k)])); add('AppSepAncLast', pct(sep[(ga, 'ancestor', hi_k)]))
    # margins of the criteria from the displayed mean times
    ms = [b for b, _ in BASELINES]
    mg = {key: min(r[b] for b in ms if r.get(b) is not None) / r['ctskel'] for key, r in cb.items() if key[1] in KEYS}
    check(len(mg) == sum(len(zms) for *_, zms, _ in CRITERIA) * len(KEYS) * len(CRIT_K), 'appendix text: criteria cells missing')
    add('AppCritCells', num(len(mg)), len(mg))
    words = {1: 'One', 4: 'Four', 16: 'Sixteen'}
    prev = None
    for k in CRIT_K:
        v, key = min((v, key) for key, v in mg.items() if key[3] == k)
        add(f'AppCritMin{words[k]}', sig(v), v, f'smallest margin for sets of {k} vertices')
        add(f'AppCritMin{words[k]}Which', f'{dict((c[0], c[2]) for c in CRITERIA)[key[0]]}, {TXT[key[1]]}, {key[2]}')
        check(prev is None or v < prev, 'appendix text: the smallest margin no longer shrinks as the sets grow')
        prev = v
    check(min(mg.values()) > 1, 'appendix text: TrekIndex is no longer faster than every baseline in every configuration')
    # shares of queries for which a criterion holds
    hv, hk = min((v, k) for k, v in hold.items()); Hv, Hk = max((v, k) for k, v in hold.items())
    name = dict((c[0], c[2]) for c in CRITERIA)
    add('AppHoldMin', pct(hv), hv, f'{name[hk[0]]}, {TXT[hk[1]]}, {hk[2]}-{hk[3]}')
    add('AppHoldMinWhich', f'{name[hk[0]]}, {TXT[hk[1]]}, {hk[2]}-{hk[3]}')
    add('AppHoldMax', pct(Hv), Hv); add('AppHoldMaxWhich', f'{name[Hk[0]]}, {TXT[Hk[1]]}, {Hk[2]}-{Hk[3]}')
    check(all(hold[(kind, g, 'random', CRIT_K[-1])] < hold[(kind, g, 'random', CRIT_K[0])] for kind, *_ in CRITERIA for g in KEYS),
          'appendix text: with random sets, the share of queries that hold no longer falls for every criterion and graph')
    # surgery-first on rule 3 (Table 8)
    r3 = {(g, zm, k): sf[('rule3', g, zm, k)]['gain'] for g in KEYS for zm in ZMODES for k in CRIT_K}
    rr = [v for (g, zm, k), v in r3.items() if zm == 'random']; ra = [v for (g, zm, k), v in r3.items() if zm == 'ancestor']
    add('AppSFThreeRandMin', sig(min(rr)), min(rr)); add('AppSFThreeRandMax', sig(max(rr)), max(rr))
    add('AppSFThreeAncMax', sig(max(ra)), max(ra))
    check(all(r3[(g, 'random', CRIT_K[i])] < r3[(g, 'random', CRIT_K[i + 1])] for g in KEYS for i in range(len(CRIT_K) - 1)),
          'appendix text: on rule 3 with random sets, the gain no longer grows with the size of the set on every graph')
    slow = sorted((v, g, k) for (g, zm, k), v in r3.items() if zm == 'ancestor' and v < 1)
    add('AppSFThreeSlowerCount', num(len(slow)), len(slow), 'rule 3, ancestor cells below 1')
    add('AppSFThreeSlowerMin', sig(slow[0][0]) if slow else RED, slow[0][0] if slow else None)
    add('AppSFThreeSlowerWhich', '; '.join(f'{TXT[g]}, ancestor-{k}' for _, g, k in slow) if slow else RED)
    # rule 2 gains more with ancestor sets, rule 3 with random sets, for every graph and size (Figure fig:app:sf)
    r2 = {(g, zm, k): sf[('rule2', g, zm, k)]['gain'] for g in KEYS for zm in ZMODES for k in CRIT_K}
    check(all(r2[(g, 'ancestor', k)] > r2[(g, 'random', k)] for g in KEYS for k in CRIT_K),
          'appendix text: on rule 2, ancestor sets no longer gain more than random sets for every graph and size')
    check(all(r3[(g, 'random', k)] > r3[(g, 'ancestor', k)] for g in KEYS for k in CRIT_K),
          'appendix text: on rule 3, random sets no longer gain more than ancestor sets for every graph and size')
    return M


def appendix_analysis_macros(e2, cb, ans, ab):
    """Further analysis of the appendix figures (fastest baselines, answer shares, edge tests); every claim of the
    text is checked here."""
    M = []

    def add(name, val, raw=None, note=''):
        check(re.fullmatch(r'[A-Za-z]+', name), f'bad macro name {name}')
        M.append((name, val, raw, note))
    sep, hold = ans
    name = dict((c[0], c[2]) for c in CRITERIA)
    BA = ('syn-ba-500k', 'syn-ba-2000k')
    # d-separation: which baseline is the fastest (Figure 5)
    fb = {key: c['fastest_baseline'] for key, c in e2['cells'].items()}
    check(set(fb.values()) <= {'bidirectional', 'bidirectional-anc', 'cifly'}, f'appendix text: fastest d-separation baselines {set(fb.values())}')
    anc = [key for key, b in fb.items() if b == 'bidirectional-anc']; cif = [key for key, b in fb.items() if b == 'cifly']
    check(all(g in BA for g, _, _ in anc), 'appendix text: Anc-BB is no longer the fastest baseline only on the BA graphs')
    add('AppDsepAncFastest', num(len(anc)), len(anc)); add('AppDsepAncMaxZ', num(max(k for *_, k in anc)), max(k for *_, k in anc))
    check(len(cif) == 1, 'appendix text: CIfly is no longer the fastest baseline in exactly one configuration')
    g, zm, k = cif[0]; add('AppDsepCiflyWhich', f'{TXT[g]}, {zm}, $|Z|={k}$')
    lo_k, hi_k = ZSIZES[0], ZSIZES[-1]
    check(all(e2['cells'][(g, zm, hi_k)]['speedup_vs_fastest'] < e2['cells'][(g, zm, lo_k)]['speedup_vs_fastest'] for g in KEYS for zm in ZMODES),
          'appendix text: the d-separation speedup no longer shrinks from |Z| = 1 to 32 for both kinds of sets on every graph')
    # criteria: fastest baselines, ancestor vs random margins, parent and child sets (Figures 10-16)
    base = [b for b, _ in BASELINES]
    best = {key: min((r[b], b) for b in base if r.get(b) is not None) for key, r in cb.items() if key[1] in KEYS}
    from collections import Counter
    cnt = Counter(b for _, b in best.values())
    check(set(cnt) <= {'bidirectional', 'bidirectional-anc', 'cifly'}, f'appendix text: fastest criteria baselines {set(cnt)}')
    add('AppCritBidirFastest', num(cnt['bidirectional']), cnt['bidirectional'])
    add('AppCritAncFastest', num(cnt['bidirectional-anc']), cnt['bidirectional-anc'])
    small = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine'}
    add('AppCritCiflyFastest', small.get(cnt['cifly'], num(cnt['cifly'])), cnt['cifly'])
    mg = {key: v / cb[key]['ctskel'] for key, (v, _) in best.items()}
    pairs = [(kind, g, k) for kind, *_ in CRITERIA for g in KEYS for k in CRIT_K]
    smaller = sum(mg[(kind, g, 'ancestor', k)] < mg[(kind, g, 'random', k)] for kind, g, k in pairs)
    add('AppCritAncSmaller', num(smaller), smaller); add('AppCritPairs', num(len(pairs)), len(pairs))
    check(smaller > len(pairs) / 2, 'appendix text: ancestor sets no longer give smaller margins than random sets in most combinations')
    for zm, mac in (('parents', 'Parent'), ('children', 'Child')):
        vs = [v for key, v in mg.items() if key[2] == zm]
        check(len(vs) == len(KEYS) * len(CRIT_K), f'appendix text: {zm} configurations')
        add(f'AppCrit{mac}Min', sig(min(vs)), min(vs)); add(f'AppCrit{mac}Max', sig(max(vs)), max(vs))
    # shares of queries that hold: parent sets of the back-door criterion, child sets of the front-door criterion
    check(all(hold[('backdoor', g, 'parents', CRIT_K[-1])] > hold[('backdoor', g, 'parents', CRIT_K[0])] for g in KEYS),
          'appendix text: with parent sets, the back-door criterion no longer holds more often for larger sets on every graph')
    ch = [v for (kind, g, zm, k), v in hold.items() if zm == 'children']
    add('AppHoldChildMin', f'{100 * min(ch):.0f}', min(ch)); add('AppHoldChildMax', f'{100 * max(ch):.0f}', max(ch))
    # effect of the index in every configuration (Figures fig:app:idx, fig:app:tests)
    cells = [ab[(g, zm, k)] for g in KEYS for zm in ZMODES for k in ZSIZES]
    tr = [c['trav_us'] / c['bi_us'] for c in cells]; se = [c['edge_us'] / c['bi_us'] for c in cells]
    ie = [c['edge_us'] / c['idx_us'] for c in cells]
    add('AppAblaTravMin', sig(min(tr)), min(tr)); add('AppAblaTravMax', sig(max(tr)), max(tr))
    add('AppAblaEdgeMin', sig(min(se)), min(se)); add('AppAblaEdgeMax', sig(max(se)), max(se))
    check(min(se) > 1, 'appendix text: trek search is no longer slower than BB in every configuration')
    add('AppAblaIdxEdgeMin', sig(min(ie)), min(ie)); add('AppAblaIdxEdgeMax', sig(max(ie)), max(ie))
    check(all(ab[(g, zm, ZSIZES[i])]['tests'] < ab[(g, zm, ZSIZES[i + 1])]['tests'] for g in KEYS for zm in ZMODES for i in range(len(ZSIZES) - 1)),
          'appendix text: the edge tests per query no longer grow with |Z| on every graph')
    t1 = max(ab[(g, zm, lo_k)]['tests'] for g in KEYS for zm in ZMODES)
    t32 = [ab[(g, zm, hi_k)]['tests'] for g in KEYS for zm in ZMODES]
    add('AppTestsFirstMax', f'{t1:.2f}', t1); add('AppTestsLastMin', f'{min(t32):.2f}', min(t32)); add('AppTestsLastMax', f'{max(t32):.1f}', max(t32))
    sl = [ab[(g, zm, hi_k)]['edge_us'] / ab[(g, zm, lo_k)]['edge_us'] for g in KEYS for zm in ZMODES]
    check(min(sl) > 1, 'appendix text: trek search is no longer slower at |Z| = 32 than at |Z| = 1 on every graph')
    add('AppEdgeSlowMin', sig(min(sl)), min(sl)); add('AppEdgeSlowMax', sig(max(sl)), max(sl))
    return M


def app_plain_macros(e1):
    """Appendix text on witness contraction (tab:app:plain): the construction of TrekIndex against that of 2-hop labels
    on the full trek DAG; every claim of the text is checked here."""
    M = []

    def add(name, val, raw=None, note=''):
        check(re.fullmatch(r'[A-Za-z]+', name), f'bad macro name {name}')
        M.append((name, val, raw, note))
    BA = ('syn-ba-500k', 'syn-ba-2000k')
    per = [g for g in KEYS if g not in BA]
    check(all(e1[g]['peeled_frac'] == 0 for g in BA) and all(e1[g]['peeled_frac'] > 0 for g in per),
          'appendix text: the graphs without a periphery are no longer exactly the BA graphs')
    # construction on the graphs with a periphery: build time and index size both shrink
    for what, mac, f in (('build time', 'Build', lambda g: e1[g]['plain']['build_s'] / e1[g]['build_s']),
                         ('index size', 'Size', lambda g: e1[g]['plain']['label_bytes'] / e1[g]['index_bytes'])):
        rs = sorted((f(g), g) for g in per)
        check(rs[0][0] > 1, f'appendix text: witness contraction no longer reduces the {what} on every graph with a periphery')
        add(f'AppPlain{mac}Min', sig(rs[0][0]), rs[0][0], f'plain / TrekIndex {what}'); add(f'AppPlain{mac}MinGraph', TXT[rs[0][1]], rs[0][1])
        add(f'AppPlain{mac}Max', sig(rs[-1][0]), rs[-1][0]); add(f'AppPlain{mac}MaxGraph', TXT[rs[-1][1]], rs[-1][1])
    # BA graphs: nothing peeled, same labels; build time and index size within a few percent
    check(all(e1[g]['plain']['label_entries'] == e1[g]['label_entries'] for g in BA), 'appendix text: the BA labels now differ')
    dev = max(abs(x - 1) for g in BA for x in (e1[g]['build_s'] / e1[g]['plain']['build_s'],
                                                e1[g]['index_bytes'] / e1[g]['plain']['label_bytes']))
    add('AppPlainBADevPct', f'{math.ceil(100 * dev)}', dev, 'max relative difference of build time and index size, BA graphs')
    return M


def write_numbers(M):
    L = [header('results/{round2,round2-extra,round2d,round3,e2-merged}'),
         '% Macros named ...Sec / ...Min carry no unit. Speedups carry no \\times sign.']
    for name, val, raw, note in M:
        L.append(f'\\newcommand{{\\{name}}}{{{val}}}' + (f'  % {note}' if note else ''))
    write('numbers.tex', L)


# --------------------------------------------------------------------- main
def jsonable(e1, e2, e3, apps, M, missing, abl, et, sf, clo):
    def ck(t):
        return '|'.join(str(x) for x in t)
    return {
        'generated_by': SCRIPT,
        'conventions': 'per cell: mean over seeds of per-seed means, mean over seeds of per-seed p99; speedup = baseline mean / '
                       'index mean; times in microseconds unless the key says _s or _bytes',
        'graphs': [{'key': g, 'table': TAB[g], 'text': TXT[g], 'macro': MAC[g], 'kind': KIND[g]} for g in KEYS],
        'e1_construction': e1,
        'e2_dsep': {'cells': {ck(k): v for k, v in e2['cells'].items()}, 'checked_runs': e2['runs'], 'mismatches': e2['mismatches'],
                    'rp_runs': e2['rp_runs'], 'rp_runs_separated_count_agrees': e2['rp_count_agree']},
        'e3_criteria': {'cells': {ck(k): v for k, v in sorted(e3['cells'].items(), key=lambda kv: ck(kv[0]))},
                        'note': 'key kind|graph|zmode|k; rule2/rule3 = surgery-first on (ctskel-sf), rule2-off/rule3-off = off (ctskel)',
                        'checked_runs': e3['runs'], 'mismatches': e3['mismatches']},
        'e5_applications': apps,
        'e6_ablation_coretree_vs_plain': {'note': 'ratio = plain 2-hop (zskel) mean / core-tree (ctskel) mean, means over seeds 1-3',
                                          'cells': {ck(k): v for k, v in abl.items()}},
        'edge_tests': {'note': 'seed 1; H_Z edge tests of the core-tree index by class; definitions in load_edgetests()',
                       'cells': {ck(k): v for k, v in et.items()}},
        'exp5_surgery_first': {'note': 'gain = index time off / on, geometric mean of both method orders; load_surgery()',
                               'cells': {ck(k): v for k, v in sf.items()}},
        'exp5_up_closures': {'note': 'ratio = index time without / with up-closures, geometric mean of both method orders; load_closure()',
                             'cells': {ck(k): v for k, v in clo.items()}},
        'macros': [{'name': n, 'latex': v, 'value': r, 'note': no} for n, v, r, no in M],
        'missing_sources': missing,
        'sources': sorted(USED),
    }


def main():
    try:
        e1 = load_e1(); e2 = load_e2(); e3 = load_e3()
        for g in KEYS:
            check(e2['nm'][g] == (e1[g]['n'], e1[g]['m']), f'{g}: n/m differ between round2 and e2-merged')
        missing = []
        apps = load_apps(missing)
        abl = load_ablation(); et = load_edgetests(); sf = load_surgery(); clo = load_closure(); ist = load_index_stats(e1)
        OUT.mkdir(parents=True, exist_ok=True)
        tab_datasets(e1); tab_build(e1); tab_dsep(e2); tab_p99(e2); dat_dsep(e2); tab_criteria(e3); tab_surgery(sf, clo)
        dat_crit(e3); tab_apps(apps); dat_ablation(abl); dat_edgetests(et); tab_decomp(e1, abl); dat_decomp(e1, abl)
        tab_app_plain(e1)
        cb = load_crit_bars(e3); dat_crit_bars(cb)
        ans = load_answers(); tab_dsep_full(e2, ans[0]); tab_crit_full(cb, ans[1])
        ab = load_index_ablation(); tab_index_ablation(ab); dat_appendix(ans[0], ans[1], ab, sf)
        M = (macros(e1, e2, e3, apps, abl, et, sf, clo, cb, ist) + appendix_macros(e2, cb, sf, ans) + index_ablation_macros(ab)
             + appendix_analysis_macros(e2, cb, ans, ab) + app_plain_macros(e1))
        write_numbers(M)
        (OUT / 'paper-data.json').write_text(json.dumps(jsonable(e1, e2, e3, apps, M, missing, abl, et, sf, clo), indent=1) + '\n')
    except SourceError as e:
        sys.exit(f'make_paper_artifacts: {e}')
    for m in missing:
        print('MISSING (red ?? emitted):', m, file=sys.stderr)
    print('wrote', ', '.join(sorted(p.name for p in OUT.iterdir())))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Generate the tables, plot data and text macros of the paper from the results.')
    ap.add_argument('--results', type=pathlib.Path, default=RES, help='directory with one sub-directory per run (default: %(default)s)')
    ap.add_argument('--out', type=pathlib.Path, default=OUT, help='output directory (default: %(default)s)')
    ap.add_argument('--replay-stats', type=pathlib.Path, help='causenet-mediation.stats.json (default: <results>/causenet-mediation.stats.json)')
    a = ap.parse_args()
    RES, OUT, REPLAY_STATS = a.results.resolve(), a.out, a.replay_stats
    main()
