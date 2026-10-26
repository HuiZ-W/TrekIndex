#!/usr/bin/env python3
"""Independent brute-force check of the trek-DAG reduction to the query graph H_Z (paper, Section 4).

Oracle: Lauritzen ancestral moralization. Candidate: connectivity of H_Z whose
edges are Z-avoiding treks, searched explicitly. Pure Python, no project code.
"""
# Brute-force check: d-sep(x,y|Z) == x,y disconnected in H_Z, where H_Z has nodes {x,y} ∪ Z and
# an edge (a,b) iff the trek DAG minus Z has a path start(a) -> end(b).
# Trek DAG: U_v -> U_p (go up to parent), U_v -> D_v (switch at top), D_v -> D_c (go down).
# start(x)=U_x, end(y)=D_y; for z in Z: start(z)=U_z without switching at z, end(z)=D_z reached from a parent.
# Interior states of Z vertices are forbidden.
import itertools, random

def ancestors(n, pa, S):
    seen=set(S); st=list(S)
    while st:
        v=st.pop()
        for p in pa[v]:
            if p not in seen: seen.add(p); st.append(p)
    return seen

def dsep_moral(n, pa, x, y, Z):
    A=ancestors(n, pa, {x,y}|Z)
    adj={v:set() for v in A}
    for v in A:
        ps=[p for p in pa[v]]
        for p in ps: adj[v].add(p); adj[p].add(v)
        for p,q in itertools.combinations(ps,2): adj[p].add(q); adj[q].add(p)
    seen={x}; st=[x]
    while st:
        v=st.pop()
        if v==y: return False
        for w in adj[v]:
            if w not in Z and w not in seen: seen.add(w); st.append(w)
    return True

def trek_reach(n, pa, ch, Z, a, b):
    # start at U_a; if a in Z must take a parent edge first (no switch at a, no stay)
    # target: D_b; if b in Z it must be entered from a parent (D_b only reached via D_p->D_b, p not in Z)
    start=[('U',a)]
    seen=set(start); st=list(start)
    while st:
        s,v=st.pop()
        nxt=[]
        if s=='U':
            for p in pa[v]: nxt.append(('U',p))
            if v not in Z: nxt.append(('D',v))
        else:
            for c in ch[v]: nxt.append(('D',c))
        for (t,w) in nxt:
            if (t,w)==('D',b) and not (b in Z and t=='D' and w==b and s=='U' and v==b):
                return True
            if w in Z: continue   # Z vertices cannot be interior
            if (t,w) not in seen: seen.add((t,w)); st.append((t,w))
    return False

def dsep_skeleton(n, pa, ch, x, y, Z):
    nodes=[x,y]+sorted(Z)
    comp={x}; st=[x]
    while st:
        a=st.pop()
        for b in nodes:
            if b in comp: continue
            if trek_reach(n,pa,ch,Z,a,b) or trek_reach(n,pa,ch,Z,b,a):
                comp.add(b); st.append(b)
    return y not in comp

def random_dag(n, p, rng):
    order=list(range(n)); rng.shuffle(order)
    pa=[[] for _ in range(n)]; ch=[[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<p:
                u,v=order[i],order[j]; pa[v].append(u); ch[u].append(v)
    return pa,ch

def main():
    rng = random.Random(7); checks = 0
    for _ in range(600):                       # exhaustive (x, y, Z), n <= 7
        n = rng.randint(2, 7); pa, ch = random_dag(n, rng.choice([0.2, 0.35, 0.5, 0.7]), rng)
        for x in range(n):
            for y in range(x + 1, n):
                rest = [v for v in range(n) if v not in (x, y)]
                for r in range(len(rest) + 1):
                    for Z in itertools.combinations(rest, r):
                        Z = set(Z); checks += 1
                        assert dsep_moral(n, pa, x, y, Z) == dsep_skeleton(n, pa, ch, x, y, Z), (n, pa, x, y, Z)
    rng = random.Random(11)
    for _ in range(4000):                      # sampled, n = 8..14
        n = rng.randint(8, 14); pa, ch = random_dag(n, rng.choice([0.1, 0.2, 0.3, 0.5]), rng)
        for _ in range(30):
            x, y = rng.sample(range(n), 2); rest = [v for v in range(n) if v not in (x, y)]
            Z = set(rng.sample(rest, rng.randint(0, len(rest)))); checks += 1
            assert dsep_moral(n, pa, x, y, Z) == dsep_skeleton(n, pa, ch, x, y, Z), (n, pa, x, y, Z)
    print(f"Z-skeleton reduction: {checks} checks, 0 mismatches")

if __name__ == "__main__":
    main()
