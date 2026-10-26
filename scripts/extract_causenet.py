#!/usr/bin/env python3
"""Stream causenet-full.jsonl.bz2 (Heindorf et al., CIKM 2020; CC BY 4.0)
into a TSV of causal relations: cause <TAB> effect <TAB> support <TAB> source types.
support = number of provenance records (sentences, list entries, infoboxes).
The archive is decompressed on the fly; nothing large is written."""
import bz2, json, sys

def main(src, dst):
    n = bad = loops = 0
    with bz2.open(src, "rt", encoding="utf-8") as f, open(dst, "w", encoding="utf-8") as out:
        for line in f:
            try:
                r = json.loads(line)
                c = r["causal_relation"]["cause"]["concept"]
                e = r["causal_relation"]["effect"]["concept"]
            except Exception:
                bad += 1
                continue
            if c == e:
                loops += 1
                continue
            srcs = r.get("sources", [])
            types = sorted({s.get("type", "?") for s in srcs})
            c, e = c.replace("\t", " "), e.replace("\t", " ")
            out.write(f"{c}\t{e}\t{len(srcs)}\t{','.join(types)}\n")
            n += 1
            if n % 1000000 == 0:
                print(f"{n} relations", file=sys.stderr, flush=True)
    print(json.dumps({"relations": n, "malformed": bad, "self_loops": loops}))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
