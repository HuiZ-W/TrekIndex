#pragma once
#include "zskeleton.hpp"
#include <unordered_map>
#include <unordered_set>
#include <queue>

namespace rs {
// Core–tree Z-skeleton index = TrekIndex (paper, Section 5: witness contraction of the
// periphery, 2-hop labels on the core), width w (1–4 in practice).
//
// Min-degree elimination on the undirected skeleton of G with a width cap:
// while some vertex has at most `width` alive neighbours, eliminate it. Its
// portals are those neighbours (they get a fill edge when there are two). The
// core K is the set of never-eliminated vertices.
//
// Eliminating v contracts the trek DAG: every T-path through v's four states
// between two different portals becomes a *shortcut arc* between portal
// states, with a witness (v's states used, plus the witnesses of the arcs it
// composes). Detours that return to the same portal are dropped (Lemma A:
// replaceable by the fork U_c -> D_c, or forbidden). The 2-hop labels are
// built on T[K] plus the shortcut arcs between core states. At query time a
// shortcut is blocked when a witness state is forbidden; the certificate adds
// the tails of blocked arcs to its mask set, the search skips them.
//
// A fringe (peeled) state is handled by projection: a DFS over the local arcs
// of peeled vertices towards their portals, recording the reachable states of
// every portal up to the core. Two states connect at a common projected state,
// or through the core between projected core states (core machinery).
struct CoreTreeZSkeleton {
  struct Arc { uint32_t from,to; int32_t wit; };           // full-graph state ids; wit: witness node or -1
  struct Wit { int v; uint8_t mask; int32_t left,right; }; // v's states used, and the composed arcs' witnesses
  const Graph *g;
  int width;
  std::vector<int> core_id;         // >= 0: core index; < 0: peeled, index -core_id-1
  std::vector<int> core_vertex;     // core index -> vertex
  std::vector<Arc> arcs;std::vector<Wit> wits;   // shortcut arcs only (direct T arcs are derived from the graph)
  // Per peeled vertex p: its chain = transitive closure of its portals, in
  // BFS order. Flags: bit 0 = direct portal of p, bit 1 = that portal is a
  // parent of p (else a child; fill-edge portals carry no direct arcs).
  // ch_fwd/ch_bwd: F-free reachable states of the chain vertex, 4 bits per
  // start/end kind (H2H-style ancestor arrays, adapted; see project()).
  std::vector<uint32_t> ch_off;std::vector<int> ch_v;std::vector<uint8_t> ch_flag;std::vector<uint16_t> ch_fwd,ch_bwd;
  std::vector<uint8_t> chain_wit;                 // per p: some chain vertex (p included) has local shortcut arcs
  std::vector<uint32_t> lsc_off,lsc_arc;          // per p: local shortcut arcs (ids), incident to p's states
  Shortcuts sc;std::vector<uint32_t> cw_off,cw_arc;   // per p: core shortcut arcs it witnesses
  bool use_anc_masks=true;
  mutable ZSkeletonStats *st_=nullptr;   // projection counters (set per query by run_h)
  size_t peeled=0,max_depth=0,fill_edges=0;
  double peel_ms=0;
  Graph core;        // declared after everything build() fills, before the labels
  ZSkeleton inner;
  // Node-level up-closure budget: closure_budget + closure_per_node * h states
  // (0 disables it); used only when h >= closure_min_h and the closure leaves
  // at most closure_max_down down-phase start states.
  size_t closure_budget=8,closure_per_node=4,closure_min_h=6,closure_max_down=4;
  // Output-sensitive probing (kept as an ablation switch; measured slower than
  // the closure + label merges). 0 = off.
  size_t probe_threshold=0;
  // Start the H_Z exploration from the endpoint side whose live states graph surgery
  // restricts, when only one side is (rule 2 cuts the out-edges of Z, rule 3 the in-edges
  // of Z(W)); H_Z is undirected, so either side gives the same answer. Off (--surgery-first 0):
  // always from X, the behaviour of git 3e4cc65.
  bool surgery_first=true;
  // Try only dominant projected states on core routes (off only for the ablation).
  bool prune_dominated=true;
  // An endpoint whose parents are all conditioned (or have their out-edges cut) has no open
  // path that leaves or enters it through a parent: the parent would be a conditioned
  // non-collider. It is treated like an in-cut endpoint (start D, end U), which sharpens the
  // T test and, with surgery_first, starts the exploration from its side.
  bool parents_blocked_live=false;   // under evaluation (2026-10-14)
  // Bytes of the fringe structures (everything except the core labels).
  size_t fringe_bytes()const{return arcs.size()*sizeof(Arc)+wits.size()*sizeof(Wit)+4*(core_id.size()+ch_off.size()+ch_v.size()+lsc_off.size()+lsc_arc.size()+cw_off.size()+cw_arc.size())+ch_flag.size()+2*(ch_fwd.size()+ch_bwd.size())+chain_wit.size();}
  std::string fringe_breakdown()const {
    auto mb=[](size_t b){return std::to_string(b>>20);};
    return "{\"arcs\":"+mb(arcs.size()*sizeof(Arc))+",\"wits\":"+mb(wits.size()*sizeof(Wit))+",\"chain\":"+mb(4*(ch_off.size()+ch_v.size())+ch_flag.size()+2*(ch_fwd.size()+ch_bwd.size()))
      +",\"lsc\":"+mb(4*(lsc_off.size()+lsc_arc.size()))+",\"cw\":"+mb(4*(cw_off.size()+cw_arc.size()))+",\"core_id\":"+mb(4*core_id.size())
      +",\"n_arcs\":"+std::to_string(arcs.size())+",\"n_chain\":"+std::to_string(ch_v.size())+"}";
  }
  bool is_core(int v)const{return core_id[v]>=0;}
  int pidx(int v)const{return -core_id[v]-1;}
  int n()const{return g->n;}

  // w = 0 selects the width automatically (see choose_width).
  CoreTreeZSkeleton(const Graph &graph,bool progress=false,const std::string &order="degree",bool compress=false,int w=0,int bp=0)
    :g(&graph),width(w>0?w:choose_width(graph)),core_id(graph.n,-1),core(build(graph,progress)),inner(core,progress,order,compress,&sc,bp) {}
  // Automatic width: wider peeling shrinks the core labels but every eliminated
  // chain adds shortcut arcs at its portals, which the core search and the
  // closure must scan. Measured on eight graphs, the shortcut arcs pay off only
  // when the graph is tree-like enough that the fringe dominates: width 3 when
  // pendant peeling alone removes at least 80 % of the vertices, else width 1.
  static int choose_width(const Graph &gr) {
    const int n=gr.n;std::vector<uint32_t> deg(n);std::vector<int> st;
    for(int v=0;v<n;++v){deg[v]=uint32_t(gr.pa[v].size()+gr.ch[v].size());if(deg[v]<=1)st.push_back(v);}
    std::vector<uint8_t> alive(n,1);size_t peeled=0;
    while(!st.empty()){int v=st.back();st.pop_back();if(!alive[v])continue;alive[v]=0;++peeled;
      for(int p:gr.pa[v])if(alive[p]&&--deg[p]<=1)st.push_back(p);
      for(int c:gr.ch[v])if(alive[c]&&--deg[c]<=1)st.push_back(c);}
    return n>0&&peeled>=size_t(0.8*n)?3:1;
  }

  struct Scratch {
    ZSkeleton::Scratch cs;                       // core scratch (marks on core ids)
    std::vector<uint32_t> zmark,inmark,outmark;  // full-graph marks, epochs below
    uint32_t zep=0,cutep=0;
    std::vector<uint32_t> arc_block;uint32_t arc_epoch=0;   // blocked shortcut arcs (global arc ids)
    std::vector<uint32_t> dedup;uint32_t dep=0;
    std::vector<int> cz,ci,co;                   // core ids of the query's Z / cut vertices (reused buffers)
    std::vector<std::vector<uint32_t>> fwd,bwd;  // cached projections per H node: reachable states (full ids)
    std::vector<uint8_t> fwd_ready,bwd_ready;
    std::vector<uint32_t> tmp_fwd,tmp_bwd,dfs;
    std::vector<int> nodes;std::vector<uint8_t> conditioned,side,reached,upblocked;
    std::vector<int> stack_nodes;
    // Up-closure per core H node (see up_closure): 0 = not computed, 1 = complete, 2 = over budget.
    std::vector<uint8_t> clo_state;std::vector<std::vector<uint32_t>> clo_up,clo_down;
    std::vector<uint32_t> upvis,upstack;uint32_t up_epoch=0;
    bool clo_hit=false,no_live=false;            // set by the last edge test (statistics only)
    // Output-sensitive exploration: per node, the mask-set indices of its
    // reachable end states (probes); per discovered node, the OR of the row
    // bits of its projected start states.
    std::vector<std::vector<uint32_t>> probe_idx;std::vector<uint64_t> rowacc,rowtmp;bool probing=false;
    // Cached certificate bits per node (row of its single projected core start
    // state, column of its single projected core end state), as in ZSkeleton.
    std::vector<uint64_t> rowc,colc;std::vector<uint8_t> rowc_ready,colc_ready;
    Scratch(int core_n,int n,size_t narcs):cs(core_n),zmark(n),inmark(n),outmark(n),arc_block(narcs),dedup(n),upvis(4*size_t(core_n)) {
      cs.arc_block=arc_block.data();
    }
    bool in_z(int v)const{return zmark[v]==zep;}
    bool in_cut(int v)const{return inmark[v]==cutep;}
    bool out_cut(int v)const{return outmark[v]==cutep;}
  };
  Scratch make_scratch()const{return Scratch(core.n,g->n,arcs.size());}

  bool separated(int x,int y,const std::vector<int> &zs,Scratch &s,ZSkeletonStats *st=nullptr)const {
    return separated_general({x},{y},zs,Surgery{},s,st);
  }
  bool separated_sets(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,Scratch &s,ZSkeletonStats *st=nullptr)const {
    return separated_general(xs,ys,zs,Surgery{},s,st);
  }
  bool separated_general(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut,Scratch &s,ZSkeletonStats *st=nullptr)const {
    if(xs.empty()||ys.empty())throw std::invalid_argument("empty endpoint set");
    mark(zs,cut,s);
    if(++s.dep==0){std::fill(s.dedup.begin(),s.dedup.end(),0);s.dep=1;}
    for(const auto *set:{&xs,&ys})for(int v:*set) {
      if(v<0||v>=g->n||s.in_z(v)||s.dedup[v]==s.dep)throw std::invalid_argument("endpoint sets must be disjoint and outside Z");
      s.dedup[v]=s.dep;
    }
    if(st)++st->queries;
    s.nodes.clear();s.conditioned.clear();s.side.clear();
    for(int v:xs){s.nodes.push_back(v);s.conditioned.push_back(0);s.side.push_back(1);}
    for(int v:ys){s.nodes.push_back(v);s.conditioned.push_back(0);s.side.push_back(2);}
    for(int v:zs){if(s.dedup[v]==s.dep)continue;s.dedup[v]=s.dep;   // duplicate Z entries are legal input
      s.nodes.push_back(v);s.conditioned.push_back(1);s.side.push_back(0);}
    return run_h(s,st);
  }
  // v is a proper descendant of u in G.
  bool is_descendant(int u,int v,Scratch &s)const {
    if(u<0||v<0||u>=g->n||v>=g->n)throw std::invalid_argument("vertex");
    if(u==v)return false;
    mark({},Surgery{},s);
    return connect_states(s,state(u,kDown),state(v,kDown),nullptr);
  }
  // Some w in targets has a directed path u -> ... -> w whose interior avoids `avoid`.
  bool reaches_avoiding(int u,const std::vector<int> &targets,const std::vector<int> &avoid,Scratch &s,ZSkeletonStats *st=nullptr)const {
    if(targets.empty())return false;
    mark(avoid,Surgery{},s);
    for(int w:targets)if(connect_states(s,state(u,kDown),state(w,kEnter),st))return true;
    return false;
  }
  // Batched form for rule 3: out[i] = reaches_avoiding(us[i], targets, avoid).
  // One marking pass, one mask build and one backward projection per target,
  // shared by all sources.
  void reaches_avoiding_many(const std::vector<int> &us,const std::vector<int> &targets,const std::vector<int> &avoid,
                             std::vector<uint8_t> &out,Scratch &s,ZSkeletonStats *st=nullptr)const {
    out.assign(us.size(),0);
    if(targets.empty()||us.empty())return;
    mark(avoid,Surgery{},s);
    st_=st;
    // Target projections are computed lazily: the scan over targets stops at
    // the first hit, so most targets are never projected.
    if(s.bwd.size()<targets.size())s.bwd.resize(targets.size());
    s.bwd_ready.assign(targets.size(),0);
    for(size_t i=0;i<us.size();++i) {
      const uint32_t a=state(us[i],kDown);
      project(s,a,true,s.tmp_fwd);
      for(size_t j=0;j<targets.size()&&!out[i];++j) {
        const uint32_t t=state(targets[j],kEnter);
        if(!s.bwd_ready[j]){s.bwd_ready[j]=1;project(s,t,false,s.bwd[j]);}
        out[i]=connect_proj(s,a,t,s.tmp_fwd,s.bwd[j],st);
      }
    }
  }

 private:
  // ---- construction -------------------------------------------------------
  static inline uint64_t ekey(uint32_t a,uint32_t b){if(a>b)std::swap(a,b);return (uint64_t(a)<<32)|b;}
  Graph build(const Graph &gr,bool progress) {
    auto start=Clock::now();
    const int n=gr.n;
    std::vector<uint32_t> deg(n);
    for(int v=0;v<n;++v)deg[v]=uint32_t(gr.pa[v].size()+gr.ch[v].size());
    std::vector<uint8_t> alive(n,1);
    std::unordered_map<int,std::vector<int>> fill;                 // fill-edge neighbours (skeleton only)
    std::unordered_set<uint64_t> fillset;
    std::unordered_map<int,std::vector<uint32_t>> inc;             // alive vertex -> incident shortcut arc ids
    std::vector<int> order;order.reserve(n);
    std::priority_queue<std::pair<uint32_t,int>,std::vector<std::pair<uint32_t,int>>,std::greater<>> pq;
    for(int v=0;v<n;++v)pq.push({deg[v],v});
    std::vector<int> nb;std::vector<Arc> ins,outs;
    struct Portal{int v;uint8_t kind;};                            // kind 0 parent, 1 child, 2 fill (no direct arcs)
    std::vector<std::vector<Portal>> portals(n);std::vector<std::vector<uint32_t>> lsc_lists(n);
    while(!pq.empty()) {
      const uint32_t d=pq.top().first;const int v=pq.top().second;pq.pop();
      if(!alive[v]||deg[v]!=d)continue;
      if(d>uint32_t(width))break;
      nb.clear();auto &pv=portals[v];
      for(int p:gr.pa[v])if(alive[p]){nb.push_back(p);pv.push_back({p,0});}
      for(int c:gr.ch[v])if(alive[c]){nb.push_back(c);pv.push_back({c,1});}
      if(auto it=fill.find(v);it!=fill.end())for(int u:it->second)if(alive[u]){nb.push_back(u);pv.push_back({u,2});}
      std::sort(nb.begin(),nb.end());nb.erase(std::unique(nb.begin(),nb.end()),nb.end());
      // Arcs incident to v in the contracted trek graph: direct T arcs to alive neighbours (transient) plus stored shortcuts.
      ins.clear();outs.clear();
      for(int p:gr.pa[v])if(alive[p]){ins.push_back({state(p,kDown),state(v,kEnter),-1});outs.push_back({state(v,kLeave),state(p,kUp),-1});}
      for(int c:gr.ch[v])if(alive[c]){ins.push_back({state(c,kLeave),state(v,kUp),-1});outs.push_back({state(v,kDown),state(c,kEnter),-1});}
      if(auto it=inc.find(v);it!=inc.end()){
        for(uint32_t id:it->second){const Arc &a=arcs[id];if(vertex_of(a.to)==v)ins.push_back(a);else outs.push_back(a);lsc_lists[v].push_back(id);}
        inc.erase(it);
      }
      // Compose: in-arc alpha->sigma, internal sigma ~> sigma', out-arc sigma'->beta, different portals.
      for(const Arc &i:ins)for(const Arc &o:outs) {
        if(vertex_of(i.from)==vertex_of(o.to))continue;           // self-detour: dropped (Lemma A)
        uint32_t s1=i.to&3,s2=o.from&3;uint8_t mask=0;
        if(s1==s2)mask=uint8_t(1u<<s1);
        else if(s1==kUp&&(s2==kLeave||s2==kDown))mask=uint8_t(1u<<kUp|1u<<s2);
        else if(s1==kEnter&&s2==kDown)mask=uint8_t(1u<<kEnter|1u<<kDown);
        else continue;
        wits.push_back({v,mask,i.wit,o.wit});
        Arc na{i.from,o.to,int32_t(wits.size()-1)};
        uint32_t id=add_arc(na);
        inc[vertex_of(na.from)].push_back(id);inc[vertex_of(na.to)].push_back(id);
      }
      // Remove v's shortcut arcs from the neighbours' incidence lists.
      for(int u:nb){auto it=inc.find(u);if(it==inc.end())continue;auto &L=it->second;
        L.erase(std::remove_if(L.begin(),L.end(),[&](uint32_t id){return vertex_of(arcs[id].from)==v||vertex_of(arcs[id].to)==v;}),L.end());if(L.empty())inc.erase(it);}
      alive[v]=0;order.push_back(v);
      for(int u:nb)deg[u]--;
      // The portals become a clique in the skeleton (fill edges), as in any elimination ordering.
      for(size_t x=0;x<nb.size();++x)for(size_t y=x+1;y<nb.size();++y)
        if(fillset.insert(ekey(nb[x],nb[y])).second){fill[nb[x]].push_back(nb[y]);fill[nb[y]].push_back(nb[x]);deg[nb[x]]++;deg[nb[y]]++;fill_edges++;}
      for(int u:nb)pq.push({deg[u],u});
    }
    peeled=order.size();
    // Depth (for reporting): parent = earliest-eliminated portal, or the core.
    std::vector<uint32_t> pos(n,UINT32_MAX),depth(n,0);for(uint32_t i=0;i<order.size();++i)pos[order[i]]=i;
    for(auto it=order.rbegin();it!=order.rend();++it){int v=*it;uint32_t best=UINT32_MAX;int bp=-1;
      for(const Portal &pt:portals[v]){uint32_t pu=alive[pt.v]?UINT32_MAX-1:pos[pt.v];if(pu<best){best=pu;bp=pt.v;}}
      depth[v]=bp<0||alive[bp]?1:depth[bp]+1;max_depth=std::max<size_t>(max_depth,depth[v]);}
    // Core and peeled indices.
    int k=0,np=0;for(int v=0;v<n;++v){if(alive[v]){core_id[v]=k++;core_vertex.push_back(v);}else core_id[v]=-(++np);}
    std::vector<std::pair<int,int>> es;
    for(int v=0;v<n;++v)if(alive[v])for(int c:gr.ch[v])if(alive[c])es.emplace_back(core_id[v],core_id[c]);
    // Shortcut arcs between core states, and the peeled vertices witnessing them.
    sc.out.assign(4*size_t(k),{});sc.in.assign(4*size_t(k),{});sc.count=0;
    std::vector<std::vector<uint32_t>> cw(np);
    for(auto &[u,L]:inc)for(uint32_t id:L){const Arc &a=arcs[id];if(vertex_of(a.from)!=u)continue;   // each arc once, from its tail
      uint32_t cf=state(core_id[vertex_of(a.from)],a.from&3),ct=state(core_id[vertex_of(a.to)],a.to&3);
      sc.out[cf].push_back({ct,id});sc.in[ct].push_back({cf,id});sc.count++;
      witnesses(a.wit,[&](int wv){cw[pidx(wv)].push_back(id);});}
    cw_off.assign(np+1,0);for(int q=0;q<np;++q){auto &L=cw[q];std::sort(L.begin(),L.end());L.erase(std::unique(L.begin(),L.end()),L.end());cw_off[q+1]=cw_off[q]+uint32_t(L.size());}
    cw_arc.resize(cw_off[np]);for(int q=0;q<np;++q)std::copy(cw[q].begin(),cw[q].end(),cw_arc.begin()+cw_off[q]);
    // Local shortcut arcs and chains (portal closure, BFS order).
    lsc_off.assign(np+1,0);ch_off.assign(np+1,0);
    std::vector<int> chain;std::vector<uint8_t> flags;
    for(int v=0;v<n;++v){if(alive[v])continue;const int q=pidx(v);
      lsc_off[q+1]=lsc_off[q]+uint32_t(lsc_lists[v].size());lsc_arc.insert(lsc_arc.end(),lsc_lists[v].begin(),lsc_lists[v].end());
      chain.clear();flags.clear();
      auto add=[&](int u,uint8_t f){for(size_t i=0;i<chain.size();++i)if(chain[i]==u){flags[i]|=f;return;}chain.push_back(u);flags.push_back(f);};
      for(const Portal &pt:portals[v])add(pt.v,pt.kind==2?0:uint8_t(1|(pt.kind==0?2:0)));
      for(size_t i=0;i<chain.size();++i)if(!alive[chain[i]])for(const Portal &pt:portals[chain[i]])add(pt.v,0);
      ch_off[q+1]=ch_off[q]+uint32_t(chain.size());ch_v.insert(ch_v.end(),chain.begin(),chain.end());ch_flag.insert(ch_flag.end(),flags.begin(),flags.end());
    }
    ch_fwd.assign(ch_v.size(),0);ch_bwd.assign(ch_v.size(),0);
    chain_wit.assign(np,0);
    for(int v=0;v<n;++v){if(alive[v])continue;const int q=pidx(v);bool w=lsc_off[q+1]>lsc_off[q];
      for(uint32_t i=ch_off[q];i<ch_off[q+1]&&!w;++i){int u=ch_v[i];if(!alive[u]){int r=pidx(u);w=lsc_off[r+1]>lsc_off[r];}}
      chain_wit[q]=w;}
    // Ancestor masks: projections with nothing forbidden, per start/end kind (the DFS path, masks off).
    {
      Scratch tmp(k,n,arcs.size());tmp.zep=1;tmp.cutep=1;tmp.arc_epoch=1;   // fresh epochs: no marks
      const bool saved=use_anc_masks;use_anc_masks=false;
      std::vector<uint32_t> out;
      for(int v=0;v<n;++v){if(alive[v])continue;const int q=pidx(v);
        auto at=[&](int u){for(uint32_t i=ch_off[q];i<ch_off[q+1];++i)if(ch_v[i]==u)return i;throw std::logic_error("projection left the chain");};
        for(uint32_t kd=0;kd<4;++kd){
          project(tmp,state(v,kd),true,out);for(uint32_t x:out){int u=vertex_of(x);if(u!=v)ch_fwd[at(u)]|=uint16_t((1u<<(x&3))<<(4*kd));}
          project(tmp,state(v,kd),false,out);for(uint32_t x:out){int u=vertex_of(x);if(u!=v)ch_bwd[at(u)]|=uint16_t((1u<<(x&3))<<(4*kd));}
        }
      }
      use_anc_masks=saved;
    }
    peel_ms=elapsed(start);
    if(progress)std::cerr<<"{\"event\":\"peeled\",\"width\":"<<width<<",\"peeled\":"<<peeled<<",\"core_n\":"<<k<<",\"core_m\":"<<es.size()<<",\"shortcuts\":"<<sc.count<<",\"arcs\":"<<arcs.size()<<",\"fill\":"<<fill_edges<<",\"max_depth\":"<<max_depth<<",\"ms\":"<<peel_ms<<"}\n";
    return Graph(k,es);
  }
  uint32_t add_arc(const Arc &a){arcs.push_back(a);return uint32_t(arcs.size()-1);}
  template<class F> void witnesses(int32_t w,F f)const {
    while(w>=0){const Wit &x=wits[w];f(x.v);if(x.right>=0)witnesses(x.right,f);w=x.left;}
  }

  // ---- query marks ----------------------------------------------------------
  // Forbidden kinds of a vertex: all four if conditioned, P/E if its in-edges are cut, U/D if its out-edges are cut.
  uint8_t forbidden_mask(const Scratch &s,int v)const {
    if(s.in_z(v))return 15;
    uint8_t m=0;
    if(s.in_cut(v))m|=uint8_t(1u<<kLeave|1u<<kEnter);
    if(s.out_cut(v))m|=uint8_t(1u<<kUp|1u<<kDown);
    return m;
  }
  bool forbidden(const Scratch &s,uint32_t sigma)const{return forbidden_mask(s,vertex_of(sigma))>>(sigma&3)&1;}
  // A shortcut arc is blocked iff some witness state is forbidden.
  bool arc_blocked(const Scratch &s,int32_t w)const {
    while(w>=0){const Wit &x=wits[w];if(forbidden_mask(s,x.v)&x.mask)return true;if(x.right>=0&&arc_blocked(s,x.right))return true;w=x.left;}
    return false;
  }
  bool arc_open(const Scratch &s,uint32_t id)const{return s.arc_block[id]!=s.arc_epoch;}

  void mark(const std::vector<int> &zs,const Surgery &cut,Scratch &s)const {
    if(++s.zep==0){std::fill(s.zmark.begin(),s.zmark.end(),0);s.zep=1;}
    for(int v:zs){if(v<0||v>=g->n)throw std::invalid_argument("conditioning vertex");s.zmark[v]=s.zep;}
    if(++s.cutep==0){std::fill(s.inmark.begin(),s.inmark.end(),0);std::fill(s.outmark.begin(),s.outmark.end(),0);s.cutep=1;}
    for(int v:cut.cut_in){if(v<0||v>=g->n)throw std::invalid_argument("surgery vertex");s.inmark[v]=s.cutep;}
    for(int v:cut.cut_out){if(v<0||v>=g->n)throw std::invalid_argument("surgery vertex");s.outmark[v]=s.cutep;}
    // Core marks: only core vertices; the fringe is handled by projection.
    auto &cz=s.cz,&ci=s.ci,&co=s.co;cz.clear();ci.clear();co.clear();
    for(int v:zs)if(is_core(v))cz.push_back(core_id[v]);
    for(int v:cut.cut_in)if(is_core(v))ci.push_back(core_id[v]);
    for(int v:cut.cut_out)if(is_core(v))co.push_back(core_id[v]);
    auto &c=s.cs;
    if(++c.z.epoch==0){std::fill(c.z.mark.begin(),c.z.mark.end(),0);c.z.epoch=1;}
    for(int v:cz)c.z.mark[v]=c.z.epoch;
    if(++c.cut_epoch==0){std::fill(c.cut_in.begin(),c.cut_in.end(),0);std::fill(c.cut_out.begin(),c.cut_out.end(),0);c.cut_epoch=1;}
    for(int v:ci)c.cut_in[v]=c.cut_epoch;
    for(int v:co)c.cut_out[v]=c.cut_epoch;
    c.cut_in_list=ci;c.cut_out_list=co;
    if(++c.forbid_epoch==0){std::fill(c.forbid.begin(),c.forbid.end(),0);c.forbid_epoch=1;}
    c.reset_blocked();
    for(int v:cz){inner.add_blocked(c,state(v,kUp));inner.add_blocked(c,state(v,kEnter));}
    for(int v:ci){inner.add_blocked(c,state(v,kLeave));inner.add_blocked(c,state(v,kEnter));}
    for(int v:co){inner.add_blocked(c,state(v,kUp));inner.add_blocked(c,state(v,kDown));}
    // Blocked shortcut arcs: witnessed by a peeled conditioned or cut vertex.
    if(++s.arc_epoch==0){std::fill(s.arc_block.begin(),s.arc_block.end(),0);s.arc_epoch=1;}
    c.arc_epoch=s.arc_epoch;
    auto block_from=[&](int v){if(is_core(v))return;const int q=pidx(v);
      for(uint32_t i=cw_off[q];i<cw_off[q+1];++i){uint32_t id=cw_arc[i];if(s.arc_block[id]==s.arc_epoch)continue;
        if(arc_blocked(s,arcs[id].wit)){s.arc_block[id]=s.arc_epoch;const Arc &a=arcs[id];
          inner.add_mask_only(c,cstate(a.from));inner.add_mask_only(c,cstate(a.to));}}};
    for(int v:zs)block_from(v);
    for(int v:cut.cut_in)block_from(v);
    for(int v:cut.cut_out)block_from(v);
    c.nodes.clear();c.masks_ready=false;
  }

  // Every parent of v is conditioned or has its out-edges cut (v has at least one parent).
  bool parents_blocked(const Scratch &s,int v)const {
    if(g->pa[v].empty())return false;
    for(int p:g->pa[v])if(!s.in_z(p)&&!s.out_cut(p))return false;
    return true;
  }
  uint32_t live_start(const Scratch &s,int i)const {
    int v=s.nodes[i];
    if(s.conditioned[i])return s.in_cut(v)?UINT32_MAX:state(v,kLeave);
    bool in=s.in_cut(v)||s.upblocked[i],out=s.out_cut(v);
    if(in&&out)return UINT32_MAX;
    return out?state(v,kLeave):in?state(v,kDown):state(v,kUp);
  }
  uint32_t live_end(const Scratch &s,int j)const {
    int v=s.nodes[j];
    if(s.conditioned[j])return s.in_cut(v)?UINT32_MAX:state(v,kEnter);
    bool in=s.in_cut(v)||s.upblocked[j],out=s.out_cut(v);
    if(in&&out)return UINT32_MAX;
    return out?state(v,kEnter):in?state(v,kUp):state(v,kDown);
  }

  // ---- projections ----------------------------------------------------------
  // Forward projection of state a: all states reachable from a by T-paths whose
  // interior avoids forbidden states and whose vertices are peeled (contracted
  // local arcs), stopping at core states. The DFS visits the states of the
  // portal chain only, so it is O(depth · width). Core states are recorded and
  // not expanded; `out` lists every reached state (a included).
  void project(Scratch &s,uint32_t a,bool forward,std::vector<uint32_t> &out)const {
    out.clear();out.push_back(a);
    const int va=vertex_of(a);
    if(is_core(va))return;   // a core start/end stands for itself; the labels cover its internal moves
    const int qa=pidx(va);
    if(use_anc_masks&&!chain_wit[qa]&&forbidden_mask(s,va)==0) {
      // Fast path: nothing forbidden on the chain and no witnessed arcs, so the
      // F-free masks are the exact projection.
      bool clean=true;for(uint32_t i=ch_off[qa];i<ch_off[qa+1]&&clean;++i)clean=forbidden_mask(s,ch_v[i])==0;
      if(clean) {
        static const uint8_t self_f[4]={1u<<kUp|1u<<kLeave|1u<<kDown,1u<<kLeave,1u<<kDown,1u<<kEnter|1u<<kDown};
        static const uint8_t self_b[4]={1u<<kUp,1u<<kLeave|1u<<kUp,1u<<kDown|1u<<kUp|1u<<kEnter,1u<<kEnter};
        const uint32_t kd=a&3;uint8_t sm=(forward?self_f:self_b)[kd];
        for(uint32_t k=0;k<4;++k)if((sm>>k&1)&&k!=kd)out.push_back(state(va,k));
        for(uint32_t i=ch_off[qa];i<ch_off[qa+1];++i){uint8_t m=uint8_t(((forward?ch_fwd:ch_bwd)[i]>>(4*kd))&15);
          for(uint32_t k=0;k<4;++k)if(m>>k&1)out.push_back(state(ch_v[i],k));}
        if(st_)++st_->proj_fast;
        return;
      }
    }
    if(st_)++st_->proj_dfs;
    auto seen=[&](uint32_t x){return std::find(out.begin(),out.end(),x)!=out.end();};
    auto &st=s.dfs;st.assign(1,a);
    while(!st.empty()) {
      uint32_t w=st.back();st.pop_back();
      if(w!=a&&forbidden(s,w))continue;             // may be an end, never an interior
      int v=vertex_of(w);
      auto push=[&](uint32_t x){if(!seen(x)){out.push_back(x);st.push_back(x);}};
      // internal moves of v (also at core portals: the junction test needs the closure there)
      uint32_t k=w&3;
      if(forward){if(k==kUp){push(state(v,kLeave));push(state(v,kDown));}else if(k==kEnter)push(state(v,kDown));}
      else {if(k==kDown){push(state(v,kUp));push(state(v,kEnter));}else if(k==kLeave)push(state(v,kUp));}
      if(is_core(v))continue;                       // core states are not expanded beyond their own vertex
      const int q=pidx(v);
      // direct arcs to v's portals (derived from the edge direction)
      for(uint32_t i=ch_off[q];i<ch_off[q+1];++i){if(!(ch_flag[i]&1))continue;const int p=ch_v[i];
        if(ch_flag[i]&2){   // p is a parent of v: P_v -> U_p, D_p -> E_v
          if(forward){if(k==kLeave)push(state(p,kUp));}else{if(k==kEnter)push(state(p,kDown));}
        } else {            // p is a child of v: D_v -> E_p, P_p -> U_v
          if(forward){if(k==kDown)push(state(p,kEnter));}else{if(k==kUp)push(state(p,kLeave));}
        }}
      // local shortcut arcs
      for(uint32_t i=lsc_off[q];i<lsc_off[q+1];++i){const Arc &ar=arcs[lsc_arc[i]];
        if(forward){if(ar.from==w&&!arc_blocked(s,ar.wit))push(ar.to);}
        else {if(ar.to==w&&!arc_blocked(s,ar.wit))push(ar.from);}}
    }
  }
  // a reaches t in T' avoiding forbidden interior states, given projections.
  bool connect_proj(Scratch &s,uint32_t a,uint32_t t,const std::vector<uint32_t> &fa,const std::vector<uint32_t> &bt,ZSkeletonStats *st,bool core_route=true,int i=-1,int j=-1)const {
    // Junction: a state in both projections that is allowed as interior (or is a or t).
    for(uint32_t x:bt)if(std::find(fa.begin(),fa.end(),x)!=fa.end()&&(x==a||x==t||!forbidden(s,x))){if(st)++st->fringe_resolved;return true;}
    if(!core_route)return false;   // no projected end state of t is label-reachable from a (probe miss)
    // Core route between projected core states. Only *dominant* states are
    // tried: in the forward projection a state reached from another usable
    // state of the same vertex by internal moves (P, D after U; D after E) adds
    // no path, since that state can continue through it as an interior state;
    // symmetrically backwards (U, E before D; U before P). A node whose
    // projection keeps exactly one state reuses that state's row/column bits.
    auto usable=[&](uint32_t x,uint32_t end){return is_core(vertex_of(x))&&(x==end||!forbidden(s,x));};
    auto has=[&](const std::vector<uint32_t> &L,uint32_t x,uint32_t end){return std::find(L.begin(),L.end(),x)!=L.end()&&usable(x,end);};
    auto dominated_fwd=[&](uint32_t x){const int v=vertex_of(x);const uint32_t k=x&3;
      if(!prune_dominated)return false;
      if((k==kLeave||k==kDown)&&x!=a&&has(fa,state(v,kUp),a))return true;
      if(k==kDown&&x!=a&&has(fa,state(v,kEnter),a))return true;
      return false;};
    auto dominated_bwd=[&](uint32_t y){const int v=vertex_of(y);const uint32_t k=y&3;
      if(!prune_dominated)return false;
      if((k==kUp||k==kEnter)&&y!=t&&has(bt,state(v,kDown),t))return true;
      if(k==kUp&&y!=t&&has(bt,state(v,kLeave),t))return true;
      return false;};
    size_t nx=0,ny=0;for(uint32_t x:fa)nx+=usable(x,a)&&!dominated_fwd(x);for(uint32_t y:bt)ny+=usable(y,t)&&!dominated_bwd(y);
    bool any=false;
    for(uint32_t x:fa){if(!usable(x,a)||dominated_fwd(x))continue;
      for(uint32_t y:bt){if(!usable(y,t)||dominated_bwd(y)||vertex_of(y)==vertex_of(x))continue;
        if(!any){any=true;if(st)++st->core_routes;}
        const uint32_t cx=cstate(x),cy=cstate(y);
        if(!inner.reach(s.cs,cx,cy,st)){if(st)++st->neg_reach;continue;}
        inner.ensure_masks(s.cs);const size_t W=s.cs.words;
        const uint64_t *r=nullptr,*c=nullptr;
        if(i>=0&&nx==1){if(s.rowc.size()<size_t(s.nodes.size())*W)s.rowc.resize(size_t(s.nodes.size())*W);r=&s.rowc[size_t(i)*W];
          if(!s.rowc_ready[i]){s.rowc_ready[i]=1;inner.row_bits(s.cs,cx,&s.rowc[size_t(i)*W],st);}}
        if(j>=0&&ny==1){if(s.colc.size()<size_t(s.nodes.size())*W)s.colc.resize(size_t(s.nodes.size())*W);c=&s.colc[size_t(j)*W];
          if(!s.colc_ready[j]){s.colc_ready[j]=1;inner.col_bits(s.cs,cy,&s.colc[size_t(j)*W],st);}}
        if(inner.connect(s.cs,cx,cy,r,c,st,true))return true;}}
    return false;
  }
  uint32_t cstate(uint32_t full)const{return state(core_id[vertex_of(full)],full&3);}
  bool connect_states(Scratch &s,uint32_t a,uint32_t t,ZSkeletonStats *st)const {
    project(s,a,true,s.tmp_fwd);project(s,t,false,s.tmp_bwd);
    return connect_proj(s,a,t,s.tmp_fwd,s.tmp_bwd,st);
  }

  // ---- edge test ---------------------------------------------------------------
  // Edge test, classified for the statistics by its costliest step: a search,
  // a certificate, an up-closure lookup, label refutations only, or the
  // projections alone (a junction, no live state, or no usable core state).
  bool edge(Scratch &s,int i,int j,ZSkeletonStats *st)const {
    if(!st)return edge_test(s,i,j,nullptr);
    const uint64_t n0=st->neg_reach,c0=st->cert_free,q0=st->searches;
    s.clo_hit=s.no_live=false;
    const bool r=edge_test(s,i,j,st);
    if(st->searches>q0)++st->et_searched;
    else if(st->cert_free>c0)++st->et_certified;
    else if(s.clo_hit)++st->et_closure;
    else if(st->neg_reach>n0&&!s.no_live)++st->et_refuted;
    else ++st->et_proj;
    return r;
  }
  bool edge_test(Scratch &s,int i,int j,ZSkeletonStats *st)const {
    if(st)++st->edge_checks;
    const uint32_t a=live_start(s,i),t=live_end(s,j);
    if(a==UINT32_MAX||t==UINT32_MAX){if(st)++st->neg_reach;s.no_live=true;return false;}
    if(!s.fwd_ready[i]){s.fwd_ready[i]=1;project(s,a,true,s.fwd[i]);}
    if(!s.bwd_ready[j]){s.bwd_ready[j]=1;project(s,t,false,s.bwd[j]);}
    // Probe: is any projected end state of j label-reachable from a projected
    // start state of i? If not, only a junction inside the fringe can connect them.
    bool hit=true;
    if(s.probing){hit=false;for(uint32_t b:s.probe_idx[j])if(s.rowacc[b/64]>>(b%64)&1){hit=true;break;}if(hit&&st)++st->candidates;}
    if(closure_budget&&s.nodes.size()>=closure_min_h&&is_core(s.nodes[i])&&up_closure(s,i,a,st))return edge_via_closure(s,i,j,a,t,st,hit);
    return connect_proj(s,a,t,s.fwd[i],s.bwd[j],st,hit,i,j);
  }
  // Row bits of node i: OR over its projected start states of the mask-set
  // states reachable from them (one label scan each).
  void probe_rows(Scratch &s,int i,uint32_t a,ZSkeletonStats *st)const {
    const size_t W=s.cs.words;s.rowacc.assign(W,0);s.rowtmp.resize(W);
    for(uint32_t x:s.fwd[i]){if(!is_core(vertex_of(x))||(x!=a&&forbidden(s,x)))continue;
      inner.row_bits(s.cs,cstate(x),s.rowtmp.data(),st);for(size_t q=0;q<W;++q)s.rowacc[q]|=s.rowtmp[q];
      if(st)++st->probe_scans;}
  }
  // Bounded DFS over the U/P core states reachable from a: the up-phase of
  // every trek from a. Turns (U_r -> D_r) and shortcut arcs into the D/E layer
  // give the down-phase start states. A complete closure turns all edge tests
  // from node i into down-phase tests: one closure instead of one search per
  // target, which is what large conditioning sets around a need.
  // Attempt the closure only when the start looks walled in: at least half of
  // (the first 16 of) its parents are forbidden. Random conditioning sets
  // rarely pass this test, so they never pay for a failed attempt.
  bool closure_worthwhile(const Scratch &s,int v)const {
    size_t n=0,f=0;for(int p:g->pa[v]){if(++n>16)break;f+=forbidden(s,state(p,kUp));}
    return n==0||2*f>=n;
  }
  bool up_closure(Scratch &s,int i,uint32_t a,ZSkeletonStats *st)const {
    if(s.clo_state[i])return s.clo_state[i]==1;
    if(!closure_worthwhile(s,s.nodes[i])){s.clo_state[i]=2;return false;}
    if(++s.up_epoch==0){std::fill(s.upvis.begin(),s.upvis.end(),0);s.up_epoch=1;}
    auto &up=s.clo_up[i];auto &down=s.clo_down[i];up.clear();down.clear();
    const uint32_t ep=s.up_epoch;const size_t budget=closure_budget+closure_per_node*s.nodes.size();
    if((a&3)==kDown||(a&3)==kEnter){down.push_back(a);s.clo_state[i]=1;if(st)++st->closures;return true;}   // no up-phase
    auto visit=[&](uint32_t u)->bool{   // u: full state id in the U/P layer; false when over budget
      if(forbidden(s,u))return true;      // conditioned / out-cut parents end the up-phase; never end states
      uint32_t c=cstate(u);if(s.upvis[c]==ep)return true;
      s.upvis[c]=ep;up.push_back(u);s.upstack.push_back(u);return up.size()<=budget;};
    s.upstack.assign(1,a);s.upvis[cstate(a)]=ep;up.push_back(a);
    while(!s.upstack.empty()) {
      uint32_t w=s.upstack.back();s.upstack.pop_back();int v=vertex_of(w);
      if(w!=a&&forbidden(s,w))continue;
      if((w&3)==kUp) {
        down.push_back(state(v,kDown));
        if(!visit(state(v,kLeave))){s.clo_state[i]=2;return false;}
      } else {   // P_v -> U_p for core parents p (pendant parents are detours)
        if(!s.in_cut(v))for(int p:g->pa[v]){if(!is_core(p)||s.out_cut(p))continue;if(!visit(state(p,kUp))){s.clo_state[i]=2;return false;}}
      }
      // open shortcut arcs out of w: into the U/P layer they continue the up-phase, into D/E they are turns
      for(const auto &e:sc.out[cstate(w)]){if(!arc_open(s,e.arc))continue;uint32_t to=state(core_vertex[vertex_of(e.to)],e.to&3);
        if((to&3)==kUp||(to&3)==kLeave){if(!visit(to)){s.clo_state[i]=2;return false;}}
        else down.push_back(to);}   // may be forbidden: still a valid end state (filtered at use)
    }
    if(down.size()>closure_max_down){s.clo_state[i]=2;return false;}   // complete but too wide to pay off
    s.clo_state[i]=1;if(st)++st->closures;return true;
  }
  bool edge_via_closure(Scratch &s,int i,int j,uint32_t a,uint32_t t,ZSkeletonStats *st,bool hit)const {
    const auto &bt=s.bwd[j];const int vi=s.nodes[i];
    for(uint32_t x:bt)if(vertex_of(x)==vi)return connect_proj(s,a,t,s.fwd[i],bt,st,hit,i,j);   // junction inside v_i's own states
    if(!hit)return false;
    const auto &up=s.clo_up[i];
    for(uint32_t tau:bt) {
      if(!is_core(vertex_of(tau))||(tau!=t&&forbidden(s,tau)))continue;
      if((tau&3)==kUp||(tau&3)==kLeave){if(std::find(up.begin(),up.end(),tau)!=up.end()){s.clo_hit=true;return true;}continue;}
      for(uint32_t d:s.clo_down[i]){if(d==tau){s.clo_hit=true;return true;}if(forbidden(s,d))continue;if(inner.connect(s.cs,cstate(d),cstate(tau),nullptr,nullptr,st))return true;}
    }
    return false;
  }

  bool run_h(Scratch &s,ZSkeletonStats *st)const {
    st_=st;s.cs.stats=st;
    const int h=int(s.nodes.size());
    if(st)st->h_nodes+=h;
    if(s.fwd.size()<size_t(h)){s.fwd.resize(h);s.bwd.resize(h);s.clo_up.resize(h);s.clo_down.resize(h);}
    s.fwd_ready.assign(h,0);s.bwd_ready.assign(h,0);s.reached.assign(h,0);s.clo_state.assign(h,0);
    s.rowc_ready.assign(h,0);s.colc_ready.assign(h,0);
    s.upblocked.assign(h,0);
    if(parents_blocked_live)for(int i=0;i<h;++i)if(!s.conditioned[i])s.upblocked[i]=parents_blocked(s,s.nodes[i]);
    s.cs.nodes.clear();s.cs.masks_ready=false;
    // Output-sensitive mode: every node's projected end states become probes
    // in the mask set, so one scan per discovered node finds all candidates.
    s.probing=probe_threshold&&size_t(h)>=probe_threshold;
    if(s.probing){
      if(s.probe_idx.size()<size_t(h))s.probe_idx.resize(h);
      for(int j=0;j<h;++j){s.probe_idx[j].clear();const uint32_t t=live_end(s,j);if(t==UINT32_MAX)continue;
        s.bwd_ready[j]=1;project(s,t,false,s.bwd[j]);
        for(uint32_t y:s.bwd[j])if(is_core(vertex_of(y))&&(y==t||!forbidden(s,y)))s.probe_idx[j].push_back(inner.add_probe(s.cs,cstate(y)));}
      inner.ensure_masks(s.cs);
    }
    uint8_t src=1;
    if(surgery_first) {
      auto cut=[&](int side){for(int i=0;i<h;++i)if(s.side[i]==side&&(s.in_cut(s.nodes[i])||s.out_cut(s.nodes[i])||s.upblocked[i]))return true;return false;};
      if(cut(2)&&!cut(1))src=2;
    }
    const uint8_t dst=3-src;
    s.stack_nodes.clear();
    for(int i=0;i<h;++i)if(s.side[i]==src){s.reached[i]=1;s.stack_nodes.push_back(i);}
    bool separated=true;
    while(!s.stack_nodes.empty()&&separated) {
      int i=s.stack_nodes.back();s.stack_nodes.pop_back();
      if(s.probing){const uint32_t a=live_start(s,i);if(a==UINT32_MAX)continue;
        if(!s.fwd_ready[i]){s.fwd_ready[i]=1;project(s,a,true,s.fwd[i]);}probe_rows(s,i,a,st);}
      for(int pass=0;pass<2&&separated;++pass)for(int j=0;j<h&&separated;++j) {
        if(s.reached[j]||(pass==0)!=(s.side[j]==dst))continue;
        if(edge(s,i,j,st)){
#ifdef CT_DEBUG
          std::cerr<<"edge "<<s.nodes[i]<<" -> "<<s.nodes[j]<<" (closure state "<<int(s.clo_state[i])<<")\n";
#endif
          if(s.side[j]==dst)separated=false;s.reached[j]=1;s.stack_nodes.push_back(j);}
      }
    }
    return separated;
  }
};
} // namespace rs
