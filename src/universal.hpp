#pragma once
#include "core.hpp"
#include <array>
#include <limits>

namespace rs {
// This compiler never receives a conditioning set. Terminals are 0 and 1.
// A node tests a vertex's membership in the query's Z.
struct CompileLimit : std::runtime_error {
  using std::runtime_error::runtime_error;
};
struct BDDNode { uint32_t var, low, high; };
struct BDDKey {
  uint32_t a, b, c;
  bool operator==(const BDDKey &o) const {
    return a == o.a && b == o.b && c == o.c;
  }
};
struct BDDHash {
  size_t operator()(const BDDKey &k) const {
    size_t h = k.a;
    h = h * 1000003ULL ^ k.b;
    return h * 1000003ULL ^ k.c;
  }
};
struct BDDBuilder {
  std::vector<BDDNode> nodes{{UINT32_MAX, 0, 0}, {UINT32_MAX, 1, 1}};
  std::unordered_map<BDDKey, uint32_t, BDDHash> unique, memo;
  size_t max_nodes, max_calls;
  size_t calls = 0, peak_memo = 0, allocated = 2, peak_nodes = 2, collections = 0;
  size_t memo_hits = 0, memo_misses = 0, gc_memo_retained = 0;
  double gc_ms = 0;
  struct ApplyFrame { uint32_t a,b,v=0,low=0,high_a=0,high_b=0; int stage=0; };
  struct RestrictFrame { uint32_t root,low=0; int stage=0; };
  std::vector<ApplyFrame> apply_stack;
  std::vector<RestrictFrame> restrict_stack;
  BDDBuilder(size_t n, size_t c) : max_nodes(n), max_calls(c) {
    if (n < 2 || n >= UINT32_MAX || !c)
      throw std::invalid_argument("invalid BDD resource limits");
  }
  void tick() {
    if (++calls > max_calls) throw CompileLimit("apply_or_restrict_call_limit");
  }
  uint32_t node(uint32_t v, uint32_t l, uint32_t h) {
    if (l == h) return l;
    BDDKey key{v, l, h};
    auto it = unique.find(key);
    if (it != unique.end()) return it->second;
    if (nodes.size() >= max_nodes) throw CompileLimit(max_nodes==UINT32_MAX-1 ? "bdd_node_id_capacity" : "bdd_node_limit");
    uint32_t id = uint32_t(nodes.size());
    nodes.push_back({v, l, h});
    ++allocated; peak_nodes = std::max(peak_nodes, nodes.size());
    unique.emplace(key, id);
    return id;
  }
  // Exact tracing collection. IDs are remapped only at safe phase boundaries;
  // the caller must supply every root that subsequent computation will use.
  void collect(std::vector<uint32_t> &roots, bool keep_memo = false) {
    auto began = Clock::now();
    std::vector<uint8_t> live(nodes.size());
    live[0] = live[1] = 1;
    for (auto r : roots) live[r] = 1;
    for (size_t id = nodes.size(); id-- > 2;)
      if (live[id]) live[nodes[id].low] = live[nodes[id].high] = 1;
    std::vector<uint32_t> remap(nodes.size());
    std::vector<BDDNode> fresh;
    fresh.reserve(std::count(live.begin(), live.end(), uint8_t(1)));
    fresh.push_back(nodes[0]); fresh.push_back(nodes[1]); remap[1] = 1;
    for (size_t id = 2; id < nodes.size(); ++id)
      if (live[id]) {
        auto v = nodes[id]; remap[id] = uint32_t(fresh.size());
        fresh.push_back({v.var, remap[v.low], remap[v.high]});
      }
    for (auto &r : roots) r = remap[r];
    nodes.swap(fresh);
    unique.clear();
    if (keep_memo) {
      // A cached operation is valid after compaction only when both operands
      // and its result survive. Cache entries do not keep garbage alive.
      decltype(memo) retained;
      for (const auto &entry : memo) {
        const auto &key = entry.first;
        if (live[key.b] && live[key.c] && live[entry.second])
          retained.emplace(BDDKey{key.a, remap[key.b], remap[key.c]},
                           remap[entry.second]);
      }
      gc_memo_retained += retained.size();
      memo.swap(retained);
    } else memo.clear();
    for (uint32_t id = 2; id < nodes.size(); ++id) {
      auto v = nodes[id]; unique.emplace(BDDKey{v.var,v.low,v.high},id);
    }
    ++collections;
    gc_ms += elapsed(began);
  }
  // op=0 AND, op=1 OR. All arguments are copied before recursive allocation.
  uint32_t apply(uint32_t op, uint32_t a, uint32_t b) {
    // Explicit reusable stacks avoid graph-dependent C++ recursion limits.
    apply_stack.clear();apply_stack.push_back({a,b});uint32_t result=0;
    while(!apply_stack.empty()) {
      auto &f=apply_stack.back();
      if(f.stage==0) {
        tick();if(f.a>f.b)std::swap(f.a,f.b);
        bool terminal=true;
        if(f.a==f.b)result=f.a;
        else if(op==0&&f.a==0)result=0;
        else if(op==0&&f.a==1)result=f.b;
        else if(op==1&&f.a==0)result=f.b;
        else if(op==1&&(f.a==1||f.b==1))result=1;
        else terminal=false;
        if(terminal){apply_stack.pop_back();continue;}
        auto found=memo.find({op,f.a,f.b});
        if(found!=memo.end()){++memo_hits;result=found->second;apply_stack.pop_back();continue;}
        ++memo_misses;
        auto na=nodes[f.a],nb=nodes[f.b];f.v=std::min(na.var,nb.var);
        uint32_t la=na.var==f.v?na.low:f.a,lb=nb.var==f.v?nb.low:f.b;
        f.high_a=na.var==f.v?na.high:f.a;f.high_b=nb.var==f.v?nb.high:f.b;
        f.stage=1;apply_stack.push_back({la,lb});
      } else if(f.stage==1) {
        f.low=result;f.stage=2;auto ha=f.high_a,hb=f.high_b;apply_stack.push_back({ha,hb});
      } else {
        result=node(f.v,f.low,result);
        if(memo.size()>=max_nodes)memo.clear();
        memo.emplace(BDDKey{op,f.a,f.b},result);peak_memo=std::max(peak_memo,memo.size());
        apply_stack.pop_back();
      }
    }
    return result;
  }
  // Restrict both endpoint variables to zero; queries reject endpoints in Z.
  uint32_t endpoints_zero(uint32_t root,uint32_t x,uint32_t y,
                         std::unordered_map<uint32_t,uint32_t>&cache) {
    restrict_stack.clear();restrict_stack.push_back({root});uint32_t result=0;
    while(!restrict_stack.empty()) {
      auto &f=restrict_stack.back();
      if(f.stage==0) {
        tick();
        if(f.root<2||nodes[f.root].var>std::max(x,y)) {
          result=f.root;restrict_stack.pop_back();continue;
        }
        auto found=cache.find(f.root);
        if(found!=cache.end()){result=found->second;restrict_stack.pop_back();continue;}
        auto n=nodes[f.root];f.stage=(n.var==x||n.var==y)?3:1;
        restrict_stack.push_back({n.low});
      } else if(f.stage==1) {
        f.low=result;f.stage=2;auto high=nodes[f.root].high;restrict_stack.push_back({high});
      } else {
        if(f.stage==2)result=node(nodes[f.root].var,f.low,result);
        cache.emplace(f.root,result);restrict_stack.pop_back();
      }
    }
    return result;
  }
};

// Query-local scratch. Stamping includes input ingestion in O(|Z|), with no
// O(n) bitmap clear except once at construction or on a 32-bit epoch wrap.
struct QueryContext {
  std::vector<uint32_t> mark;
  uint32_t epoch = 0;
  explicit QueryContext(int n) : mark(n) {}
  bool contains(int v) const { return mark[v] == epoch; }
  void set(const std::vector<int> &z, int x, int y) {
    if (x < 0 || y < 0 || x >= int(mark.size()) || y >= int(mark.size()) || x == y)
      throw std::invalid_argument("invalid endpoints");
    if (++epoch == 0) {
      std::fill(mark.begin(), mark.end(), 0);
      epoch = 1;
    }
    for (int v : z) {
      if (v < 0 || v >= int(mark.size()))
        throw std::invalid_argument("conditioning vertex");
      mark[v] = epoch;
    }
    if (contains(x) || contains(y))
      throw std::invalid_argument("endpoints must be outside Z");
  }
};

struct UniversalBDD {
  int n;
  std::vector<BDDNode> nodes;
  std::vector<uint32_t> roots; // All ordered endpoint pairs, n*n root IDs.
  size_t created_nodes = 0, apply_calls = 0, peak_memo = 0, peak_nodes = 0, collections = 0;
  std::string variable_order;
  UniversalBDD(const Graph &g, const std::string &order = "topo",
               size_t max_nodes = UINT32_MAX-1, size_t max_calls = std::numeric_limits<size_t>::max(),
               int max_vertices = 256, const std::string &guards = "ancestor")
      : n(g.n), variable_order(order) {
    // Resource guard also bounds recursion depth and the quadratic state matrix.
    if (g.n > max_vertices || g.n > 512)
      throw CompileLimit("vertex_limit");
    if (order != "topo" && order != "reverse" && order != "id")
      throw std::invalid_argument("variable order must be topo, reverse, or id");
    if (guards != "ancestor" && guards != "local")
      throw std::invalid_argument("guards must be ancestor or local");
    std::vector<int> vertices = g.order;
    if (order == "reverse") std::reverse(vertices.begin(), vertices.end());
    if (order == "id") std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<uint32_t> rank(n);
    for (int i = 0; i < n; ++i) rank[vertices[i]] = i;
    BDDBuilder b(max_nodes, max_calls);
    int states = 2 * n;
    std::vector<uint32_t> reach(size_t(states) * states);
    auto at = [&](int x, int y) -> uint32_t & { return reach[size_t(x)*states+y]; };
    auto edge = [&](int x, int y, uint32_t guard) {
      at(x,y) = b.apply(1, at(x,y), guard);
    };
    std::vector<uint32_t> activation(n);
    for (auto it = g.order.rbegin(); it != g.order.rend(); ++it) {
      int v = *it;
      activation[v] = b.node(rank[v], 0, 1);
      if (guards == "ancestor")
        for (int c : g.ch[v]) activation[v] = b.apply(1, activation[v], activation[c]);
    }
    for (int v = 0; v < n; ++v) {
      auto no = b.node(rank[v], 1, 0);
      at(2*v,2*v) = at(2*v+1,2*v+1) = 1;
      // Local guards use z_v and activate ancestors by down/up walks.
      // Ancestor guards use OR_{w in De*(v)} z_w, compiled once symbolically.
      for (int p : g.pa[v]) {
        edge(2*v, 2*p, no);
        edge(2*v+1, 2*p, activation[v]);
      }
      for (int c : g.ch[v]) {
        edge(2*v, 2*c+1, no);
        edge(2*v+1, 2*c+1, no);
      }
    }
    // Boolean-semiring Floyd-Warshall. The state graph can contain cycles.
    // Evaluation commutes with AND/OR, so this closes every Z simultaneously.
    std::vector<int> elimination;
    for (int v : g.order) { elimination.push_back(2*v); elimination.push_back(2*v+1); }
    for (int k : elimination) {
      for (int i = 0; i < states; ++i) {
        if (i == k || !at(i,k)) continue;
        auto ik = at(i,k);
        for (int j = 0; j < states; ++j) {
          if (j == k || !at(k,j) || at(i,j) == 1) continue;
          at(i,j) = b.apply(1, at(i,j), b.apply(0, ik, at(k,j)));
        }
      }
      // The closure matrix is the complete live root set at this boundary.
      if (b.nodes.size() >= std::min(size_t(10000), max_nodes/4)) b.collect(reach);
    }
    roots.resize(size_t(n)*n);
    for (int x = 0; x < n; ++x)
      for (int y = 0; y < n; ++y) {
        if (x == y) continue;
        auto r = b.apply(1, at(2*x,2*y), at(2*x,2*y+1));
        std::unordered_map<uint32_t, uint32_t> cache;
        roots[size_t(x)*n+y] = b.endpoints_zero(r, rank[x], rank[y], cache);
      }
    created_nodes = b.allocated;
    peak_nodes = b.peak_nodes;
    collections = b.collections;
    apply_calls = b.calls;
    peak_memo = b.peak_memo;
    // Keep only nodes reachable from final roots; remap rank tests to vertex IDs.
    // Construction IDs are topological: every child has a smaller ID.
    std::vector<uint8_t> live(b.nodes.size());
    live[0] = live[1] = 1;
    for (auto r : roots) live[r] = 1;
    for (size_t id = b.nodes.size(); id-- > 2;)
      if (live[id]) {
        live[b.nodes[id].low] = live[b.nodes[id].high] = 1;
      }
    std::vector<uint32_t> remap(b.nodes.size());
    nodes.reserve(std::count(live.begin(), live.end(), uint8_t(1)));
    nodes.push_back(b.nodes[0]); nodes.push_back(b.nodes[1]); remap[1] = 1;
    for (size_t id = 2; id < b.nodes.size(); ++id)
      if (live[id]) {
        auto v = b.nodes[id];
        remap[id] = uint32_t(nodes.size());
        nodes.push_back({uint32_t(vertices[v.var]), remap[v.low], remap[v.high]});
      }
    for (auto &r : roots) r = remap[r];
  }
  bool separated(int x, int y, const std::vector<int> &z, QueryContext &scratch) const {
    if (int(scratch.mark.size()) != n) throw std::invalid_argument("scratch size");
    scratch.set(z, x, y);
    uint32_t r = roots[size_t(x)*n+y];
    while (r > 1) {
      auto v = nodes[r];
      r = scratch.contains(v.var) ? v.high : v.low;
    }
    return r == 0;
  }
  size_t bytes() const {
    return sizeof(*this) + nodes.capacity()*sizeof(BDDNode) + roots.capacity()*4 +
           variable_order.capacity();
  }
};

// Fair arbitrary-Z traversal: one reusable workspace, no per-query ancestor
// computation or O(n) context initialization. Same local rules as CIflyLike.
struct StreamingCIfly {
  const Graph *g;
  QueryContext context;
  std::vector<uint32_t> seen;
  std::vector<int> todo;
  explicit StreamingCIfly(const Graph &gg) : g(&gg), context(gg.n), seen(2*gg.n) {}
  bool separated(int x, int y, const std::vector<int> &z) { return run<false>(x, y, z); }
  // Separation in G with the out-edges of x removed (back-door graph G_{\underline x}).
  bool separated_cut_out(int x, int y, const std::vector<int> &z) { return run<true>(x, y, z); }
  template<bool Cut> bool run(int x, int y, const std::vector<int> &z) {
    auto previous = context.epoch;
    context.set(z, x, y);
    if (previous == UINT32_MAX) std::fill(seen.begin(), seen.end(), 0);
    todo.clear();
    auto push = [&](int s) {
      if (seen[s] != context.epoch) {
        seen[s] = context.epoch; todo.push_back(s);
      }
    };
    push(2*x);
    while (!todo.empty()) {
      int s = todo.back(); todo.pop_back();
      int v = s/2, d = s%2;
      if (v == y) return false;
      if (!context.contains(v)) {
        if (!d) for (int p : g->pa[v]) { if (Cut && p == x) continue; push(2*p); }
        if (!(Cut && v == x)) for (int c : g->ch[v]) push(2*c+1);
      } else if (d) for (int p : g->pa[v]) { if (Cut && p == x) continue; push(2*p); }
    }
    return true;
  }
  size_t bytes() const {
    return sizeof(*this) + 4*(context.mark.capacity()+seen.capacity()+todo.capacity());
  }
};
} // namespace rs
