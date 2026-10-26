// Builds the TrekIndex (core-tree index) of the running example of the paper
// (Figure 1) with the default settings and prints what the index stores: the
// width, the core, the hub order, the 2-hop labels of the core states, and the
// chain and ancestor masks of every peeled vertex. It then answers Q1 and
// prints which layer decided its edge tests. Used for Figure 3(c).
//
//   make build/dump_running_example && ./build/dump_running_example
#include "coretree.hpp"
#include <iostream>

using namespace rs;

static const char *kKind = "UPDE";   // kUp=0, kLeave=1 (P), kDown=2, kEnter=3 (E)

int main() {
  const std::vector<std::string> name = {"Age", "Drug", "Disease", "Protein", "Gene",
                                         "Mutation", "Recovery", "SideEffect", "Biomarker",
                                         "HospitalStay"};
  auto id = [&](const std::string &s) {
    for (size_t i = 0; i < name.size(); ++i)
      if (name[i] == s) return int(i);
    throw std::invalid_argument(s);
  };
  const std::vector<std::pair<std::string, std::string>> named = {
      {"Age", "Drug"},       {"Age", "Disease"},      {"Drug", "Protein"},
      {"Protein", "Disease"}, {"Gene", "Protein"},     {"Drug", "Recovery"},
      {"Gene", "Recovery"},  {"Mutation", "Gene"},    {"Drug", "SideEffect"},
      {"Recovery", "Biomarker"}, {"Recovery", "HospitalStay"}};
  std::vector<std::pair<int, int>> edges;
  for (auto &[u, v] : named) edges.emplace_back(id(u), id(v));
  Graph g(int(name.size()), edges);

  CoreTreeZSkeleton ct(g);   // defaults of the experiments: degree order, automatic width
  const auto &L = ct.inner.labels;
  auto core_state = [&](uint32_t s) {   // core state id -> "U_Drug"
    return std::string(1, kKind[s & 3]) + "_" + name[ct.core_vertex[s / 4]];
  };

  std::cout << "width " << ct.width << ", peeled " << ct.peeled << ", shortcuts " << ct.sc.count
            << "\ncore:";
  for (int v : ct.core_vertex) std::cout << ' ' << name[v];
  std::cout << "\n\nhub order (rank: state):\n";
  for (uint32_t k = 0; k < L.states; ++k) std::cout << "  " << k << ": " << core_state(L.hub_state[k]) << '\n';

  std::cout << "\nlabels (hubs as states, in rank order):\n";
  for (uint32_t s = 0; s < L.states; ++s) {
    std::cout << "  " << core_state(s) << "  Lout = {";
    for (uint64_t i = L.out_off[s]; i < L.out_off[s + 1]; ++i)
      std::cout << (i > L.out_off[s] ? ", " : "") << core_state(L.hub_state[L.out_hub[i]]);
    std::cout << "}  Lin = {";
    for (uint64_t i = L.in_off[s]; i < L.in_off[s + 1]; ++i)
      std::cout << (i > L.in_off[s] ? ", " : "") << core_state(L.hub_state[L.in_hub[i]]);
    std::cout << "}\n";
  }

  auto mask = [&](uint16_t m, uint32_t kd) {   // 4 bits for start/end kind kd -> "{U, P}"
    std::string out = "{";
    for (int b = 0; b < 4; ++b)
      if ((m >> (4 * kd + b)) & 1) out += std::string(out.size() > 1 ? ", " : "") + kKind[b];
    return out + "}";
  };
  std::cout << "\npeeled vertices (chain, then per own state the chain states reached with F empty):\n";
  for (int v = 0; v < g.n; ++v) {
    if (ct.is_core(v)) continue;
    const int q = ct.pidx(v);
    std::cout << "  " << name[v] << ": chain {";
    for (uint32_t i = ct.ch_off[q]; i < ct.ch_off[q + 1]; ++i)
      std::cout << (i > ct.ch_off[q] ? ", " : "") << name[ct.ch_v[i]];
    std::cout << "}\n";
    for (uint32_t i = ct.ch_off[q]; i < ct.ch_off[q + 1]; ++i)
      for (uint32_t kd = 0; kd < 4; ++kd)
        std::cout << "    " << kKind[kd] << '_' << name[v] << ": fwd " << name[ct.ch_v[i]] << ' '
                  << mask(ct.ch_fwd[i], kd) << ", bwd " << name[ct.ch_v[i]] << ' '
                  << mask(ct.ch_bwd[i], kd) << '\n';
  }

  // Q1: is Drug d-separated from Gene given {Biomarker}?
  auto s = ct.make_scratch();
  ZSkeletonStats st;
  const bool sep = ct.separated(id("Drug"), id("Gene"), {id("Biomarker")}, s, &st);
  std::cout << "\nQ1 separated = " << sep << "; edge checks " << st.edge_checks << ", refuted "
            << st.neg_reach << ", certified " << st.cert_free << ", searches " << st.searches
            << ", fringe-resolved " << st.fringe_resolved << ", core routes " << st.core_routes
            << ", mask builds " << st.mask_builds << ", projections by masks " << st.proj_fast
            << ", by search " << st.proj_dfs << '\n';
  for (bool empty_z : {true}) {
    (void)empty_z;
    ZSkeletonStats st0;
    const bool sep0 = ct.separated(id("Drug"), id("Gene"), {}, s, &st0);
    std::cout << "Q1 with Z empty: separated = " << sep0 << "; edge checks " << st0.edge_checks
              << ", refuted " << st0.neg_reach << ", certified " << st0.cert_free << ", searches "
              << st0.searches << '\n';
  }
  // Q2: separation condition of the back-door criterion, Sep(Drug, Disease, {Age}; {}, {Drug}).
  {
    ZSkeletonStats st2;
    Surgery cut;
    cut.cut_out = {id("Drug")};
    const bool sep2 = ct.separated_general({id("Drug")}, {id("Disease")}, {id("Age")}, cut, s, &st2);
    std::cout << "Q2 separation condition holds = " << sep2 << "; edge checks " << st2.edge_checks
              << ", refuted " << st2.neg_reach << ", certified " << st2.cert_free << ", searches "
              << st2.searches << " (found " << st2.search_found << ", states " << st2.search_states << ")\n";
  }
  return 0;
}
