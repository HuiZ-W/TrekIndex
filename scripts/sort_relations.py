#!/usr/bin/env python3
"""Sort a relations TSV (from scripts/extract_ikraph_causal.sh) into the order used for the iKraph DAG.

tools/causenet_dag.cpp numbers vertices in order of first appearance, so the vertex IDs of the iKraph DAG
depend on the line order of its input. The paper's file was built from the relations sorted by GNU sort
under the en_US.UTF-8 locale, whose collation ignores tabs and semicolons at the first level; for these
lines (digits, tabs, semicolons) that order is: by the line with tabs and semicolons removed, ties broken
by the bytes of the line. This script reproduces it on any platform; the result is byte-identical to the
paper's ikraph-causal-dag.bin (see the SHA-256 in the dataset README).

usage: sort_relations.py IN.tsv OUT.tsv
"""
import sys


def main(src, dst):
    with open(src, 'rb') as f:
        lines = f.read().split(b'\n')
    if lines and lines[-1] == b'':
        lines.pop()
    lines.sort(key=lambda l: (l.translate(None, b'\t;'), l))
    with open(dst, 'wb') as f:
        f.write(b'\n'.join(lines) + b'\n')


if __name__ == '__main__':
    main(*sys.argv[1:3])
