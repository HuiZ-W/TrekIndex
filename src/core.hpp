#pragma once
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace rs {
using Label = std::vector<uint32_t>;
using Clock = std::chrono::steady_clock;
inline double elapsed(Clock::time_point t) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}
// Time per phase of the baselines, for the cost breakdown (single-threaded accumulators that the
// benchmark loops reset before each timed pass): sep_ns = the d-separation searches of a criterion,
// anc_ns = the ancestor marking inside them (ancestral bidirectional, RP-separation). The rest of a
// criterion's time is its directed-path (Reach) tests and bookkeeping.
struct PhaseTimes { double sep_ns = 0, anc_ns = 0; void reset() { sep_ns = anc_ns = 0; } };
inline PhaseTimes &phase_times() { static PhaseTimes p; return p; }
struct PhaseTimer {
  double &acc; Clock::time_point t0;
  explicit PhaseTimer(double &a) : acc(a), t0(Clock::now()) {}
  ~PhaseTimer() { acc += std::chrono::duration<double, std::nano>(Clock::now() - t0).count(); }
};
inline bool intersects(const Label &a, const Label &b) {
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j])
      return true;
    if (a[i] < b[j])
      ++i;
    else
      ++j;
  }
  return false;
}
inline Label unite(const Label &a, const Label &b) {
  Label c;
  c.reserve(a.size() + b.size());
  std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));
  return c;
}
struct Graph {
  int n;
  size_t m = 0;
  std::vector<std::vector<int>> pa, ch;
  std::vector<int> order, rank;
  Graph(int n_, const std::vector<std::pair<int, int>> &edges)
      : n(n_), pa(n_), ch(n_), rank(n_) {
    if (n_ < 0 || n_ > std::numeric_limits<int>::max()/2)
      throw std::invalid_argument("graph size");
    auto es = edges;
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    for (auto [u, v] : es) {
      if (u < 0 || u >= n || v < 0 || v >= n || u == v)
        throw std::invalid_argument("edge");
      pa[v].push_back(u);
      ch[u].push_back(v);
    }
    m = es.size();
    std::vector<int> d(n);
    std::deque<int> q;
    for (int v = 0; v < n; ++v) {
      d[v] = pa[v].size();
      if (!d[v])
        q.push_back(v);
    }
    while (!q.empty()) {
      int v = q.front();
      q.pop_front();
      rank[v] = order.size();
      order.push_back(v);
      for (int c : ch[v])
        if (!--d[c])
          q.push_back(c);
    }
    if ((int)order.size() != n)
      throw std::invalid_argument("cycle");
  }
  std::vector<uint8_t> context(const std::vector<int> &zs) const {
    std::vector<uint8_t> z(n);
    for (int v : zs) {
      if (v < 0 || v >= n)
        throw std::invalid_argument("conditioning vertex");
      z[v] = 1;
    }
    return z;
  }
  void validate(int x, int y, const std::vector<uint8_t> &z) const {
    if (x < 0 || y < 0 || x >= n || y >= n || x == y || z[x] || z[y])
      throw std::invalid_argument("endpoints must be distinct and outside Z");
  }
  std::vector<uint8_t> ancestors(const std::vector<uint8_t> &z) const {
    auto a = z;
    for (auto it = order.rbegin(); it != order.rend(); ++it)
      if (a[*it])
        for (int p : pa[*it])
          a[p] = 1;
    return a;
  }
};
// Two arrival directions, target-aware and allocation-reusing. Both directions
// at y accept.
struct BayesBall {
  const Graph *g;
  std::vector<uint8_t> z, a;
  std::vector<uint32_t> seen;
  uint32_t epoch = 0;
  std::vector<int> todo;
  BayesBall(const Graph &g_, const std::vector<int> &zs)
      : g(&g_), z(g_.context(zs)), a(g_.ancestors(z)), seen(2 * g_.n) {}
  size_t bytes() const {
    return sizeof(*this) + z.capacity() + a.capacity() +
           4 * (seen.capacity() + todo.capacity());
  }
  bool separated(int x, int y) {
    g->validate(x, y, z);
    if (++epoch == 0) {
      std::fill(seen.begin(), seen.end(), 0);
      epoch = 1;
    }
    todo.clear();
    auto push = [&](int s) {
      if (seen[s] != epoch) {
        seen[s] = epoch;
        todo.push_back(s);
      }
    };
    push(2 * x);
    while (!todo.empty()) {
      int s = todo.back();
      todo.pop_back();
      int v = s / 2, d = s % 2;
      if (v == y)
        return false;
      if (d == 0 && !z[v]) {
        for (int p : g->pa[v])
          push(2 * p);
        for (int c : g->ch[v])
          push(2 * c + 1);
      }
      if (d == 1) {
        if (!z[v])
          for (int c : g->ch[v])
            push(2 * c + 1);
        if (a[v])
          for (int p : g->pa[v])
            push(2 * p);
      }
    }
    return true;
  }
};
// Direct rule-table semantics: collider transitions are enabled at conditioned
// nodes; visiting a conditioned descendant and bouncing back activates its
// ancestors implicitly.
struct CIflyLike {
  const Graph *g;
  std::vector<uint8_t> z;
  std::vector<uint32_t> seen;
  uint32_t epoch = 0;
  std::vector<int> todo;
  CIflyLike(const Graph &g_, const std::vector<int> &zs)
      : g(&g_), z(g_.context(zs)), seen(2 * g_.n) {}
  size_t bytes() const {
    return sizeof(*this) + z.capacity() +
           4 * (seen.capacity() + todo.capacity());
  }
  bool separated(int x, int y) {
    g->validate(x, y, z);
    if (++epoch == 0) {
      std::fill(seen.begin(), seen.end(), 0);
      epoch = 1;
    }
    todo.clear();
    auto push = [&](int s) {
      if (seen[s] != epoch) {
        seen[s] = epoch;
        todo.push_back(s);
      }
    };
    push(2 * x);
    while (!todo.empty()) {
      int s = todo.back();
      todo.pop_back();
      int v = s / 2, d = s % 2;
      if (v == y)
        return false;
      if (!z[v]) {
        if (d == 0)
          for (int p : g->pa[v])
            push(2 * p);
        for (int c : g->ch[v])
          push(2 * c + 1);
      } else if (d == 1)
        for (int p : g->pa[v])
          push(2 * p);
    }
    return true;
  }
};
struct DSU {
  std::vector<int> p, sz;
  explicit DSU(int n) : p(n), sz(n, 1) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) {
    while (x != p[x]) {
      p[x] = p[p[x]];
      x = p[x];
    }
    return x;
  }
  void join(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y)
      return;
    if (sz[x] < sz[y])
      std::swap(x, y);
    p[y] = x;
    sz[x] += sz[y];
  }
};
struct UpdateStats {
  size_t ancestor_delta = 0, hyperedges = 0, region = 0, labels_evaluated = 0,
         labels_changed = 0;
  bool fallback = false;
};
struct MutableRoot {
  const Graph *g;
  std::vector<uint8_t> z, a;
  std::vector<int> support, cid;
  std::vector<std::vector<int>> hedge, inc;
  std::map<int, std::vector<int>> members;
  std::vector<Label> labels;
  std::vector<uint32_t> queued;
  uint32_t epoch = 0;
  MutableRoot(const Graph &gg, const std::vector<int> &zs)
      : g(&gg), z(gg.context(zs)) {
    rebuild();
  }
  Label compute(int v) const {
    if (z[v])
      return {};
    if (a[v])
      return {uint32_t(2 * cid[v])};
    Label l;
    for (int p : g->pa[v])
      l = unite(l, labels[p]);
    if (l.empty())
      l.push_back(2 * v + 1);
    return l;
  }
  std::vector<int> edge(int v) const {
    std::vector<int> e;
    if (!a[v])
      return e;
    if (!z[v])
      e.push_back(v);
    for (int p : g->pa[v])
      if (!z[p])
        e.push_back(p);
    return e;
  }
  void rebuild() {
    int n = g->n;
    a = g->ancestors(z);
    support.assign(n, 0);
    cid.assign(n, -1);
    hedge.assign(n, {});
    inc.assign(n, {});
    members.clear();
    labels.assign(n, {});
    queued.assign(n, 0);
    epoch = 0;
    DSU d(n);
    for (int v = 0; v < n; ++v) {
      support[v] = z[v];
      for (int c : g->ch[v])
        support[v] += a[c];
      hedge[v] = edge(v);
      for (int w : hedge[v]) {
        inc[w].push_back(v);
        d.join(hedge[v][0], w);
      }
    }
    std::map<int, int> minimum;
    for (int v = 0; v < n; ++v)
      if (a[v] && !z[v]) {
        int r = d.find(v);
        if (!minimum.count(r))
          minimum[r] = v;
      }
    for (int v = 0; v < n; ++v)
      if (a[v] && !z[v]) {
        cid[v] = minimum[d.find(v)];
        members[cid[v]].push_back(v);
      }
    for (int v : g->order)
      labels[v] = compute(v);
  }
  bool separated(int x, int y) const {
    g->validate(x, y, z);
    return !intersects(labels[x], labels[y]);
  }
  UpdateStats toggle(int v, double fallback_fraction = 0.5) {
    if (v < 0 || v >= g->n)
      throw std::invalid_argument("update vertex");
    UpdateStats st;
    int delta = z[v] ? -1 : 1;
    z[v] = !z[v];
    support[v] += delta;
    std::vector<int> todo{v}, da;
    while (!todo.empty()) {
      int w = todo.back();
      todo.pop_back();
      bool now = support[w] > 0;
      if (now == bool(a[w]))
        continue;
      int sign = now ? 1 : -1;
      a[w] = now;
      da.push_back(w);
      for (int p : g->pa[w]) {
        support[p] += sign;
        todo.push_back(p);
      }
    }
    st.ancestor_delta = da.size();
    std::set<int> changed(da.begin(), da.end());
    changed.insert(v);
    for (int c : g->ch[v])
      if (a[c] || !hedge[c].empty())
        changed.insert(c);
    st.hyperedges = changed.size();
    std::set<int> oldcomps;
    std::unordered_set<int> region;
    auto touch = [&](int w) {
      region.insert(w);
      if (cid[w] >= 0)
        oldcomps.insert(cid[w]);
    };
    for (int w : da)
      touch(w);
    touch(v);
    for (int h : changed) {
      for (int w : hedge[h])
        touch(w);
      for (int w : edge(h))
        touch(w);
    }
    for (int c : oldcomps)
      for (int w : members.at(c))
        region.insert(w);
    st.region = region.size();
    if (region.size() > fallback_fraction * g->n) {
      st.fallback = true;
      rebuild();
      st.labels_evaluated = g->n;
      return st;
    }
    // Delete incidence first. Forward splits require connectivity
    // recomputation, not DSU rollback.
    for (int h : changed) {
      for (int w : hedge[h]) {
        auto &es = inc[w];
        es.erase(std::find(es.begin(), es.end(), h));
      }
      hedge[h] = edge(h);
      for (int w : hedge[h])
        inc[w].push_back(h);
    }
    for (int c : oldcomps)
      members.erase(c);
    for (int w : region)
      cid[w] = -1;
    std::unordered_set<int> visited_h;
    for (int start : region)
      if (a[start] && !z[start] && cid[start] < 0) {
        std::vector<int> comp{start};
        cid[start] = start;
        int lo = start;
        for (size_t i = 0; i < comp.size(); ++i) {
          int w = comp[i];
          lo = std::min(lo, w);
          for (int h : inc[w])
            if (visited_h.insert(h).second)
              for (int u : hedge[h])
                if (cid[u] < 0) {
                  if (!region.count(u))
                    throw std::logic_error("repair escaped region");
                  cid[u] = start;
                  comp.push_back(u);
                }
        }
        for (int w : comp)
          cid[w] = lo;
        members.emplace(lo, std::move(comp));
      }
    if (++epoch == 0) {
      std::fill(queued.begin(), queued.end(), 0);
      epoch = 1;
    }
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                        std::greater<std::pair<int, int>>>
        q;
    auto enqueue = [&](int w) {
      if (queued[w] != epoch) {
        queued[w] = epoch;
        q.emplace(g->rank[w], w);
      }
    };
    for (int w : region)
      enqueue(w);
    while (!q.empty()) {
      int w = q.top().second;
      q.pop();
      ++st.labels_evaluated;
      auto l = compute(w);
      if (l != labels[w]) {
        ++st.labels_changed;
        labels[w] = std::move(l);
        for (int c : g->ch[w])
          enqueue(c);
      }
    }
    return st;
  }
  size_t bytes() const {
    size_t b = sizeof(*this) + z.capacity() + a.capacity() +
               4 * (support.capacity() + cid.capacity() + queued.capacity());
    for (auto *vv : {&hedge, &inc}) {
      b += vv->capacity() * sizeof(std::vector<int>);
      for (auto &v : *vv)
        b += 4 * v.capacity();
    }
    b += labels.capacity() * sizeof(Label);
    for (auto &v : labels)
      b += 4 * v.capacity();
    for (auto &kv : members)
      b += sizeof(kv) + 4 * kv.second.capacity() + 3 * sizeof(void *);
    return b;
  }
};
} // namespace rs
