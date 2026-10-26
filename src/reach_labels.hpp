#pragma once
#include "core.hpp"
#include <numeric>
namespace rs {
// Pruned 2-hop reachability labels on the input DAG G (Cohen et al. 2003; pruned construction as in
// Yano et al. 2013 and distribution labeling, Jin and Wang 2013), built once, single-threaded. They
// answer u ~> v in G (u = v included) and do not depend on any query. The "+R" baselines use them
// only for the directed-path (Reach) tests of a criterion: exactly where plain reachability decides
// a test, and as a refutation filter (no path at all) before their own search otherwise.
struct ReachLabels {
  std::vector<std::vector<uint32_t>> lin, lout;   // hub ranks, increasing
  double build_ms = 0;
  explicit ReachLabels(const Graph &g) : lin(g.n), lout(g.n) {
    auto t = Clock::now();
    const int n = g.n;
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    auto key = [&](int v) { return (uint64_t(g.pa[v].size()) + 1) * (uint64_t(g.ch[v].size()) + 1); };
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return key(a) > key(b); });
    std::vector<uint32_t> seen(n, 0), hub(n, 0);   // per-BFS marks; hub marks of the current hub's label
    std::vector<int> q;
    uint32_t ep = 0;
    for (int r = 0; r < n; ++r) {
      const int h = order[r];
      for (int dir = 0; dir < 2; ++dir) {   // 0: forward (fills lin), 1: backward (fills lout)
        ++ep;
        for (uint32_t x : (dir == 0 ? lout[h] : lin[h])) hub[x] = ep;
        auto covered = [&](int u) {   // h ~> u (dir 0) or u ~> h (dir 1) already answered by the labels
          for (uint32_t x : (dir == 0 ? lin[u] : lout[u])) if (hub[x] == ep) return true;
          return false;
        };
        q.assign(1, h); seen[h] = ep;
        for (size_t i = 0; i < q.size(); ++i) {
          const int u = q[i];
          if (covered(u)) continue;
          (dir == 0 ? lin[u] : lout[u]).push_back(uint32_t(r));
          for (int w : (dir == 0 ? g.ch[u] : g.pa[u])) if (seen[w] != ep) { seen[w] = ep; q.push_back(w); }
        }
      }
    }
    build_ms = elapsed(t);
  }
  static bool meet(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b) {
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) { if (a[i] == b[j]) return true; if (a[i] < b[j]) ++i; else ++j; }
    return false;
  }
  bool reach(int u, int v) const { return meet(lout[u], lin[v]); }
  size_t entries() const { size_t e = 0; for (auto &l : lin) e += l.size(); for (auto &l : lout) e += l.size(); return e; }
  size_t bytes() const { return 4 * entries() + 2 * lin.size() * sizeof(std::vector<uint32_t>); }
};

} // namespace rs
