#pragma once
#include "zskeleton.hpp"

namespace rs {
// Causal criteria answered by one graph-only Z-skeleton index, with graph
// surgery applied at query time (paper, Sections 2.2 and 4).
// Engine is ZSkeleton (pure index) or AdaptiveZSkeleton; both expose
// separated_general(xs,ys,zs,Surgery,scratch,stats). Ancestor/descendant
// side conditions always use the index labels (exact).
struct FrontdoorResult {
  bool valid=false;
  bool unintercepted=false;     // condition (i) violated: a directed X -> Y path avoids M
  bool backdoor_x_m=false;      // condition (ii) violated: open back-door path X <- ... M
  bool backdoor_m_y=false;      // condition (iii) violated: back-door path M <- ... Y not blocked by X
};
struct BackdoorResult {
  bool valid=false;
  bool descendant_in_z=false;   // condition (i) violated; witness x -> ... -> z
  bool open_backdoor_path=false;// condition (ii) violated
  int x=-1,z=-1;
};
// Index is the structure answering the side conditions (ZSkeleton or
// CoreTreeZSkeleton); Engine answers separated_general with the same Scratch.
template<class Engine,class Index=ZSkeleton> struct CausalQueries {
  using Scratch=typename std::remove_const<Index>::type::Scratch;
  const Index *index;
  Engine *engine;
  CausalQueries(const Index &idx,Engine &eng):index(&idx),engine(&eng) {}

  // Pearl's back-door criterion for (X, Y) and candidate set Z:
  // (i) no vertex of Z is a descendant of X; (ii) X ⊥ Y | Z in G_{\underline X}.
  // Sufficient for adjustment, not complete (the complete criterion of
  // Shpitser–VanderWeele–Robins / Perković et al. needs edge-level surgery on
  // proper causal paths); `valid` means "satisfies the back-door criterion".
  BackdoorResult backdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,
                    Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys,&zs});
    if(xs.empty()||ys.empty())throw std::invalid_argument("back-door needs nonempty X and Y");
    BackdoorResult r;
    if(descendant_in(xs,zs,r,s))return r;
    r.open_backdoor_path=!engine->separated_general(xs,ys,zs,Surgery{{},xs},s,st);
    r.valid=!r.open_backdoor_path;
    return r;
  }
  // Rule 1: P(y|do(x),z,w) = P(y|do(x),w)  if (Y ⊥ Z | X,W) in G_{\overline X}.
  bool rule1(const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &xs,const std::vector<int> &ws,
             Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys,&zs,&ws});
    if(ys.empty()||zs.empty())return true;
    return engine->separated_general(ys,zs,unite(xs,ws),Surgery{xs,{}},s,st);
  }
  // Rule 2: P(y|do(x),do(z),w) = P(y|do(x),z,w)  if (Y ⊥ Z | X,W) in G_{\overline X \underline Z}.
  bool rule2(const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &xs,const std::vector<int> &ws,
             Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys,&zs,&ws});
    if(ys.empty()||zs.empty())return true;
    return engine->separated_general(ys,zs,unite(xs,ws),Surgery{xs,zs},s,st);
  }
  // Rule 3: P(y|do(x),do(z),w) = P(y|do(x),w)  if (Y ⊥ Z | X,W) in
  // G_{\overline X \overline{Z(W)}}, where Z(W) = Z minus the ancestors of W
  // in G_{\overline X}.
  bool rule3(const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &xs,const std::vector<int> &ws,
             Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys,&zs,&ws});
    if(ys.empty()||zs.empty())return true;
    auto cut=xs;
    std::vector<uint8_t> anc;index->reaches_avoiding_many(zs,ws,xs,anc,s,st);
    for(size_t i=0;i<zs.size();++i)if(!anc[i])cut.push_back(zs[i]);
    return engine->separated_general(ys,zs,unite(xs,ws),Surgery{cut,{}},s,st);
  }

  // Pearl's front-door criterion for (X, Y) and M: (i) M intercepts every directed path
  // from X to Y; (ii) X ⊥ M in G_{\underline X} (no unblocked back-door path X to M);
  // (iii) M ⊥ Y | X in G_{\underline M} (X blocks every back-door path M to Y).
  FrontdoorResult frontdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ms,
                            Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys,&ms});
    if(xs.empty()||ys.empty())throw std::invalid_argument("front-door needs nonempty X and Y");
    FrontdoorResult r;
    for(int x:xs)if(index->reaches_avoiding(x,ys,ms,s,st)){r.unintercepted=true;return r;}
    if(!ms.empty()) {
      if(!engine->separated_general(xs,ms,{},Surgery{{},xs},s,st)){r.backdoor_x_m=true;return r;}
      if(!engine->separated_general(ms,ys,xs,Surgery{{},ms},s,st)){r.backdoor_m_y=true;return r;}
    }
    r.valid=true;
    return r;
  }
  // Sequential back-door criterion (Pearl and Robins 1995) for treatments xs in temporal
  // order with covariates ws[k] observed before xs[k]: for every k, ws[k] contains no
  // descendant of xs[k..n-1], and Y ⊥ x_k | x_1..x_{k-1}, W_1..W_k in
  // G_{\underline{x_k}, \overline{x_{k+1} .. x_n}}.
  bool sequential_backdoor(const std::vector<int> &xs,const std::vector<std::vector<int>> &ws,const std::vector<int> &ys,
                           Scratch &s,ZSkeletonStats *st=nullptr) const {
    if(ws.size()!=xs.size())throw std::invalid_argument("one covariate set per treatment");
    std::vector<int> all;for(auto &w:ws)all.insert(all.end(),w.begin(),w.end());
    disjoint({&xs,&ys,&all});
    if(xs.empty()||ys.empty())throw std::invalid_argument("sequential back-door needs treatments and outcomes");
    std::vector<int> cond;
    for(size_t k=0;k<xs.size();++k) {
      for(size_t j=k;j<xs.size();++j)for(int w:ws[k])if(index->is_descendant(xs[j],w,s))return false;
      cond.insert(cond.end(),ws[k].begin(),ws[k].end());
      std::vector<int> later(xs.begin()+long(k)+1,xs.end());
      if(!engine->separated_general(ys,{xs[k]},cond,Surgery{later,{xs[k]}},s,st))return false;
      cond.push_back(xs[k]);
    }
    return true;
  }
  // Natural direct and indirect effects of X on Y through M identified by the mediation
  // formula with covariates W: Pearl's conditions A-1..A-4 (2001, 2014) with the W-specific
  // effects identified by (sequential) back-door adjustment:
  //   W contains no descendant of X or M;  M ⊥ Y | X, W in G_{\underline M};
  //   X ⊥ M | W in G_{\underline X};      X ⊥ Y | W in G_{\underline X \overline M}.
  bool mediation(const std::vector<int> &xs,const std::vector<int> &ms,const std::vector<int> &ys,const std::vector<int> &ws,
                 Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ms,&ys,&ws});
    if(xs.empty()||ms.empty()||ys.empty())throw std::invalid_argument("mediation needs nonempty X, M and Y");
    for(int w:ws){for(int x:xs)if(index->is_descendant(x,w,s))return false;for(int m:ms)if(index->is_descendant(m,w,s))return false;}
    if(!engine->separated_general(ms,ys,unite(xs,ws),Surgery{{},ms},s,st))return false;
    if(!engine->separated_general(xs,ms,ws,Surgery{{},xs},s,st))return false;
    return engine->separated_general(xs,ys,ws,Surgery{ms,xs},s,st);
  }
  // s-recoverability of P(y | x) under selection bias (Bareinboim, Tian and Pearl 2014,
  // Thm. 1): recoverable without external data iff Y ⊥ S | X once the selection node S,
  // with parents ps and no children, is added to G. A path reaches S through p -> S with
  // p a non-collider, so this is Y ⊥ (ps \ X) | X in G; S itself is never materialized.
  bool selection_recoverable(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ps,
                             Scratch &s,ZSkeletonStats *st=nullptr) const {
    disjoint({&xs,&ys});
    if(ys.empty())throw std::invalid_argument("selection needs nonempty Y");
    std::unordered_set<int> in_x(xs.begin(),xs.end()),in_y(ys.begin(),ys.end()),seen;std::vector<int> ends;
    for(int p:ps) {
      if(p<0||p>=index->g->n)throw std::invalid_argument("vertex out of range");
      if(in_y.count(p))return false;   // y -> S
      if(!in_x.count(p)&&seen.insert(p).second)ends.push_back(p);
    }
    return ends.empty()||engine->separated_general(ys,ends,xs,Surgery{},s,st);
  }

 private:
  // Condition (i): some z in Z is a proper descendant of X. One label merge
  // per (x, z) pair (D_x reaches D_z in T); no traversal of G.
  bool descendant_in(const std::vector<int> &xs,const std::vector<int> &zs,BackdoorResult &r,Scratch &s) const {
    for(int x:xs)for(int z:zs)if(index->is_descendant(x,z,s)){r.descendant_in_z=true;r.x=x;r.z=z;return true;}
    return false;
  }
  static std::vector<int> unite(const std::vector<int> &a,const std::vector<int> &b) {
    auto r=a;r.insert(r.end(),b.begin(),b.end());return r;
  }
  void disjoint(std::initializer_list<const std::vector<int>*> sets) const {
    std::unordered_set<int> seen;
    for(auto *set:sets)for(int v:*set) {
      if(v<0||v>=index->g->n)throw std::invalid_argument("vertex out of range");
      if(!seen.insert(v).second)throw std::invalid_argument("causal query sets must be pairwise disjoint without duplicates");
    }
  }
};
} // namespace rs
