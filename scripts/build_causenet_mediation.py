#!/usr/bin/env python3
"""Map the published mediation hypotheses of Toonsi et al. (Bioinformatics
2026) onto CauseNet concepts and emit back-door identification checks for
`zskeleton replay` on the CauseNet DAG (tools/causenet_dag.cpp).

Name matching: lower-case; drop ICD qualifiers (", unspecified", " NOS",
"unspecified "), parenthesised text and punctuation; exact match against
CauseNet concepts. Drugs fall back to their first word ("vincristine
liposome" -> "vincristine"). Hypotheses whose three names do not all match
are skipped and counted.

Covariates W for each hypothesis (A indication, D drug, O outcome):
  observed: parents of A, D or O in CauseNet that are themselves names of a
            DD-CKG disease or drug (measurable in a cohort); main workload.
  all:      all parents of A, D or O; stress workload (CauseNet hubs have
            tens of thousands of parents).
Checks per W, as in mediation identification: A->O | W, A->D | W,
D->O | W u {A}.
"""
import argparse, array, collections, csv, json, pathlib, re

def norm(s):
    s = s.lower()
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r",?\s*unspecified\b", " ", s)
    s = re.sub(r"\bunspecified\s+", " ", s)
    s = re.sub(r"\bnos\b", " ", s)
    s = re.sub(r"[^a-z0-9 ]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True, type=pathlib.Path)
    ap.add_argument("--graph", required=True, type=pathlib.Path)
    ap.add_argument("--study", required=True, type=pathlib.Path, help="Mediation-Analysis-using-Causal-Knowledge-Graph checkout")
    ap.add_argument("--output", required=True, type=pathlib.Path)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    concept = {}
    with open(a.names, encoding="utf-8") as f:
        for line in f:
            i, name, _, _ = line.rstrip("\n").split("\t")
            concept.setdefault(norm(name), int(i))
    # Study vocabulary (what a cohort can observe) and hypotheses with names.
    rows = []
    for fn in ["results/disease_progression_significant_results.tsv", "results/comorbidity_significant_results.tsv"]:
        rows += list(csv.reader(open(a.study / fn), delimiter="\t"))[1:]
    triples = sorted({(r[1], r[3], r[5]) for r in rows})   # indication, drug, outcome names
    vocab = set()
    for ind, drug, out in triples:
        vocab |= {norm(ind), norm(out), norm(drug), norm(drug).split(" ")[0]}
    observed = {concept[v] for v in vocab if v in concept}
    def drug_id(d):
        k = norm(d)
        if k in concept: return concept[k]
        k = k.split(" ")[0] if k else k
        return concept.get(k)
    mapped = []
    for ind, drug, out in triples:
        x, m, y = concept.get(norm(ind)), drug_id(drug), concept.get(norm(out))
        if None not in (x, m, y) and len({x, m, y}) == 3:
            mapped.append((x, m, y))
    need = {v for t in mapped for v in t}
    pa = collections.defaultdict(list)
    with open(a.graph, "rb") as f:
        assert f.read(8) == b"FDGBIN01"
        n, mcount = array.array("Q"), None
        hdr = array.array("Q"); hdr.fromfile(f, 2); mcount = hdr[1]
        chunk = 1 << 22
        left = mcount
        while left:
            k = min(chunk, left); buf = array.array("Q"); buf.fromfile(f, 2 * k); left -= k
            for i in range(0, 2 * k, 2):
                if buf[i + 1] in need: pa[buf[i + 1]].append(buf[i])
    reqs = {"observed": [], "all": []}
    sizes = {"observed": [], "all": []}
    def B(lst, xs, ys, zs):
        lst.append("B %d %s %d %s %d %s" % (len(xs), " ".join(map(str, xs)), len(ys), " ".join(map(str, ys)),
                                              len(zs), " ".join(map(str, sorted(zs)))))
    for x, m, y in mapped:
        allp = (set(pa[x]) | set(pa[m]) | set(pa[y])) - {x, m, y}
        for name, w in (("observed", allp & observed), ("all", allp)):
            sizes[name].append(len(w))
            B(reqs[name], [x], [y], w); B(reqs[name], [x], [m], w); B(reqs[name], [m], [y], w | {x})
    for name, lst in reqs.items():
        (a.output / f"causenet-{name}.requests").write_text("\n".join(lst) + "\n")
    pct = lambda v, p: sorted(v)[min(len(v) - 1, int(p * len(v)))] if v else 0
    stats = {"study_unique_triples": len(triples), "mapped_triples": len(mapped),
             "observable_concepts": len(observed),
             "W_sizes": {k: {"p50": pct(v, .5), "p90": pct(v, .9), "max": max(v) if v else 0} for k, v in sizes.items()},
             "requests": {k: len(v) for k, v in reqs.items()}}
    (a.output / "causenet-mediation.stats.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps(stats))

if __name__ == "__main__":
    main()
