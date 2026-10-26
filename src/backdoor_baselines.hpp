#pragma once
#include "paper_algorithms.hpp"
#include "polytree.hpp"
#include "reach_labels.hpp"
#include "universal.hpp"
namespace rs {
// Pearl's back-door criterion for a single treatment x and outcome y with a
// single-pair d-separation baseline: (i) no z in Z is a proper descendant of x
// (graph search from x); (ii) x ⊥ y | Z in G_{\underline x}, answered by the
// baseline's own algorithm on the graph without the out-edges of x
// (Sep::separated_cut_out). Used for the application experiments, so that every
// baseline answers exactly the criterion the index answers.
template<class Sep> struct BackdoorBaseline {
  const Graph *g;
  Sep sep;
  std::vector<uint32_t> seen,in_z;uint32_t epoch=0;
  std::vector<int> stack;
  template<class... A> explicit BackdoorBaseline(const Graph &gg,A&&... a)
      : g(&gg),sep(gg,std::forward<A>(a)...),seen(gg.n),in_z(gg.n) {}
  // +R: condition (i) from the reachability labels of G (exact: plain reachability); the separation
  // test (ii) is the baseline's own, unchanged.
  const ReachLabels *R=nullptr;
  void use_labels(const ReachLabels &r){R=&r;}
  bool valid(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs) {
    if(xs.size()!=1||ys.size()!=1)throw std::invalid_argument("back-door baselines take one treatment and one outcome");
    const int x=xs[0],y=ys[0];
    if(++epoch==0){std::fill(seen.begin(),seen.end(),0);std::fill(in_z.begin(),in_z.end(),0);epoch=1;}
    if(R){for(int z:zs)if(z!=x&&R->reach(x,z))return false;PhaseTimer pt(phase_times().sep_ns);return sep.separated_cut_out(x,y,zs);}
    for(int z:zs)in_z[z]=epoch;
    stack.clear();
    for(int c:g->ch[x])if(seen[c]!=epoch){seen[c]=epoch;stack.push_back(c);}
    while(!stack.empty()) {
      const int v=stack.back();stack.pop_back();
      if(in_z[v]==epoch)return false;
      for(int c:g->ch[v])if(seen[c]!=epoch){seen[c]=epoch;stack.push_back(c);}
    }
    PhaseTimer pt(phase_times().sep_ns);
    return sep.separated_cut_out(x,y,zs);
  }
};
} // namespace rs
