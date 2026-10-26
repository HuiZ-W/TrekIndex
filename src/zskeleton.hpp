#pragma once
#include "polytree.hpp"
#include <iostream>

namespace rs {
// Graph-only arbitrary-Z index by the trek-DAG reduction (paper, Section 4); the query
// graph H_Z is called the Z-skeleton in the code.
//
// Trek DAG T, four states per vertex v, never materialized:
//   U=4v   arrived from a child (or start of an unconditioned endpoint)
//   P=4v+1 leave v upward; start of a conditioned endpoint (no switch at v)
//   D=4v+2 move down from v; accepting state of an unconditioned endpoint
//   E=4v+3 entered v from a parent; accepting state of a conditioned endpoint
// Arcs: U_v->P_v, U_v->D_v, P_v->U_p (p parent), E_v->D_v, D_v->E_c (c child).
// T is acyclic. A T-path is a collider-free trek; Z vertices may be only its
// start/end. x and y are d-connected given Z iff they are connected in the
// graph H_Z over {x,y} ∪ Z whose edges are such Z-avoiding treks.
// Labels index T once, independent of Z; Z acts as blocked states at query time.
enum : uint32_t { kUp=0, kLeave=1, kDown=2, kEnter=3 };
inline uint32_t state(int v,uint32_t kind){return 4*uint32_t(v)+kind;}
inline int vertex_of(uint32_t s){return int(s/4);}

template<class F> void trek_out(const Graph &g,uint32_t s,F f) {
  int v=vertex_of(s);
  switch(s&3) {
    case kUp: f(state(v,kLeave)); f(state(v,kDown)); break;
    case kLeave: for(int p:g.pa[v]) f(state(p,kUp)); break;
    case kDown: for(int c:g.ch[v]) f(state(c,kEnter)); break;
    default: f(state(v,kDown));
  }
}
template<class F> void trek_in(const Graph &g,uint32_t s,F f) {
  int v=vertex_of(s);
  switch(s&3) {
    case kUp: for(int c:g.ch[v]) f(state(c,kLeave)); break;
    case kLeave: f(state(v,kUp)); break;
    case kDown: f(state(v,kUp)); f(state(v,kEnter)); break;
    default: for(int p:g.pa[v]) f(state(p,kDown));
  }
}
inline size_t trek_out_degree(const Graph &g,uint32_t s) {
  int v=vertex_of(s);
  switch(s&3){case kUp:return 2;case kLeave:return g.pa[v].size();case kDown:return g.ch[v].size();default:return 1;}
}
inline size_t trek_in_degree(const Graph &g,uint32_t s) {
  int v=vertex_of(s);
  switch(s&3){case kUp:return g.ch[v].size();case kLeave:return 1;case kDown:return 2;default:return g.pa[v].size();}
}

// Extra arcs on T between core states, produced by contracting peeled vertices
// (witness contraction, paper Section 5). Each arc has a global id; a query may block
// arcs whose witness vertices are conditioned or cut. Indexed by state.
struct Shortcuts {
  struct E { uint32_t to,arc; };
  std::vector<std::vector<E>> out,in;
  size_t count=0;
  bool empty()const{return count==0;}
};

// Canonical hierarchical 2-hop reachability labels on T (pruned landmark
// labeling). Hub k is in Lout(s) iff k is the minimum-rank state among all
// states lying on s-to-hub paths, and symmetrically for Lin. Every reachable
// pair shares its minimum-rank between-state, so reach(s,t) is one merge.
struct TrekLabels {
  uint32_t states=0;
  std::vector<uint32_t> hub_state;              // rank -> state
  std::vector<uint64_t> out_off,in_off;
  std::vector<uint32_t> out_hub,in_hub;         // sorted hub ranks
  size_t bfs_visits=0;
  double build_ms=0;
  // A topological order of T: all U/P states children-before-parents (reverse
  // topological order of G), then all E/D states in topological order of G.
  static std::vector<uint32_t> topological(const Graph &g) {
    std::vector<uint32_t> t;t.reserve(4*size_t(g.n));
    for(auto it=g.order.rbegin();it!=g.order.rend();++it){t.push_back(state(*it,kUp));t.push_back(state(*it,kLeave));}
    for(int v:g.order){t.push_back(state(v,kEnter));t.push_back(state(v,kDown));}
    return t;
  }
  // Hub order. "degree": state in-degree x out-degree (distribution labeling,
  // Jin & Wang 2013). "butterfly": estimated |ancestors| x |descendants| in T,
  // the pair-coverage score of TOL's Butterfly order and O'Reach supports;
  // set sizes come from k-min sketches (one sweep per direction).
  // Hub order keys on T, or on T_w = T[core] + shortcut arcs when `sc` is set
  // (the contracted arcs count as ordinary arcs: a state with many shortcut
  // arcs is a hub of T_w even if its G-degree is small).
  static std::vector<uint64_t> order_keys(const Graph &g,const std::string &order,uint32_t states,const Shortcuts *sc=nullptr) {
    std::vector<uint64_t> key(states);
    auto outd=[&](uint32_t s){return trek_out_degree(g,s)+(sc?sc->out[s].size():0);};
    auto ind=[&](uint32_t s){return trek_in_degree(g,s)+(sc?sc->in[s].size():0);};
    auto each_out=[&](uint32_t s,auto f){trek_out(g,s,f);if(sc)for(const auto &e:sc->out[s])f(e.to);};
    auto each_in=[&](uint32_t s,auto f){trek_in(g,s,f);if(sc)for(const auto &e:sc->in[s])f(e.to);};
    if(order=="degree") {        // G-degree of the trek state (shortcut arcs ignored)
      for(uint32_t s=0;s<states;++s)key[s]=uint64_t(trek_in_degree(g,s)+1)*(trek_out_degree(g,s)+1);
      return key;
    }
    if(order=="degree-sc") {     // degree in T_w, shortcut arcs included
      for(uint32_t s=0;s<states;++s)key[s]=uint64_t(ind(s)+1)*(outd(s)+1);
      return key;
    }
    if(order!="butterfly")throw std::invalid_argument("hub order must be degree|degree-sc|butterfly");
    constexpr int K=8;
    // Topological order of the state graph (Kahn), so shortcut arcs are respected.
    std::vector<uint32_t> topo;topo.reserve(states);
    {std::vector<uint32_t> indeg(states);for(uint32_t s=0;s<states;++s)indeg[s]=uint32_t(ind(s));
     for(uint32_t s=0;s<states;++s)if(!indeg[s])topo.push_back(s);
     for(size_t i=0;i<topo.size();++i)each_out(topo[i],[&](uint32_t w){if(--indeg[w]==0)topo.push_back(w);});
     if(topo.size()!=states)throw std::logic_error("state graph is not acyclic");}
    std::vector<float> rnd(size_t(states)*K),mins(size_t(states)*K);
    std::mt19937_64 rng(0x9e3779b97f4a7c15ULL);std::uniform_real_distribution<float> u(0.f,1.f);
    for(auto &r:rnd)r=u(rng);
    auto estimate=[&](uint32_t s){double sum=0;for(int k=0;k<K;++k)sum+=mins[size_t(s)*K+k];return std::max(1.0,K/std::max(sum,1e-9)-1);};
    std::vector<double> desc(states);
    // Descendants: reverse topological order, min over successors.
    for(auto it=topo.rbegin();it!=topo.rend();++it) {
      uint32_t s=*it;float *m=&mins[size_t(s)*K];std::copy(&rnd[size_t(s)*K],&rnd[size_t(s)*K]+K,m);
      each_out(s,[&](uint32_t w){const float *o=&mins[size_t(w)*K];for(int k=0;k<K;++k)m[k]=std::min(m[k],o[k]);});
    }
    for(uint32_t s=0;s<states;++s)desc[s]=estimate(s);
    for(uint32_t s:topo) {
      float *m=&mins[size_t(s)*K];std::copy(&rnd[size_t(s)*K],&rnd[size_t(s)*K]+K,m);
      each_in(s,[&](uint32_t w){const float *o=&mins[size_t(w)*K];for(int k=0;k<K;++k)m[k]=std::min(m[k],o[k]);});
    }
    for(uint32_t s=0;s<states;++s){double score=desc[s]*estimate(s);key[s]=uint64_t(std::min(score,1.8e19));}
    return key;
  }
  std::string order_name="degree";
  double order_ms=0;
  bool compressed=false;
  size_t num_entries=0;
  std::vector<uint8_t> out_bytes,in_bytes;      // delta-LEB128 rows when compressed
  // Bit-parallel hubs (PLL, Akiba et al. SIGMOD 2013, adapted to reachability):
  // the bp_hubs top-ranked states are not stored in the label rows; instead
  // out_bp[s] has bit k iff s reaches hub k and in_bp[s] has bit k iff hub k
  // reaches s (full reachability, not pruned). reach(s, t) also holds when
  // out_bp[s] & in_bp[t] != 0. Exact for every query layer.
  int bp_hubs=0;
  std::vector<uint64_t> out_bp,in_bp;
  // FELINE-style dominance filter: two topological orders of the state graph
  // (Kahn with FIFO and with LIFO tie-breaking). s reaches t only if
  // ord1[s] < ord1[t] and ord2[s] < ord2[t]; a failed test is an O(1) negative.
  bool use_topo=false;
  std::vector<uint32_t> ord1,ord2;
  mutable uint64_t topo_pruned=0,topo_tests=0;
  // Single-threaded build, matching the single-threaded baselines.
  static inline bool topo_filter=false;   // set before construction (experiment switch)
  explicit TrekLabels(const Graph &g,bool progress=false,const std::string &order="degree",bool compress=false,const Shortcuts *sc=nullptr,int bp=0):order_name(order) {
    auto start=Clock::now();
    if(uint64_t(g.n)*4>UINT32_MAX-1)throw std::length_error("trek state id capacity");
    states=4*uint32_t(g.n);
    auto key=order_keys(g,order,states,sc);
    hub_state.resize(states);std::iota(hub_state.begin(),hub_state.end(),0u);
    std::stable_sort(hub_state.begin(),hub_state.end(),[&](uint32_t a,uint32_t b){return key[a]>key[b];});
    order_ms=elapsed(start);
    if(progress)std::cerr<<"{\"event\":\"hub_order\",\"order\":\""<<order<<"\",\"ms\":"<<order_ms<<"}\n";
    std::vector<uint32_t> rank(states);
    for(uint32_t k=0;k<states;++k)rank[hub_state[k]]=k;
    key.clear();key.shrink_to_fit();
    std::vector<std::vector<uint32_t>> lout(states),lin(states);
    std::vector<uint32_t> hub_mark(states,UINT32_MAX),seen(states,UINT32_MAX),queue;
    // Bit-parallel hubs: full forward and backward reachability of the top hubs.
    bp_hubs=int(std::min<uint32_t>(uint32_t(std::max(bp,0)),std::min<uint32_t>(64,states)));
    if(bp_hubs) {
      out_bp.assign(states,0);in_bp.assign(states,0);
      for(int k=0;k<bp_hubs;++k) {
        const uint32_t h=hub_state[k];const uint64_t bit=uint64_t(1)<<k;
        for(int direction=0;direction<2;++direction) {
          const uint32_t tag=2*uint32_t(k)+uint32_t(direction);
          auto &mark=direction==0?in_bp:out_bp;   // forward: descendants of h get in_bp; backward: ancestors get out_bp
          queue.assign(1,h);seen[h]=tag;
          for(size_t i=0;i<queue.size();++i) {
            auto v=queue[i];++bfs_visits;mark[v]|=bit;
            auto visit=[&](uint32_t w){if(seen[w]!=tag){seen[w]=tag;queue.push_back(w);}};
            if(direction==0)trek_out(g,v,visit);else trek_in(g,v,visit);
            if(sc){for(const auto &e:(direction==0?sc->out:sc->in)[v])visit(e.to);}
          }
        }
      }
      std::fill(seen.begin(),seen.end(),UINT32_MAX);
    }
    auto last=Clock::now();
    for(uint32_t k=uint32_t(bp_hubs);k<states;++k) {
      if(progress && elapsed(last)>10000) {
        std::cerr<<"{\"event\":\"pll_progress\",\"hub_rank\":"<<k<<",\"states\":"<<states
          <<",\"entries\":"<<num_entries<<",\"bfs_visits\":"<<bfs_visits<<",\"elapsed_ms\":"<<elapsed(start)<<"}\n";
        last=Clock::now();
      }
      const uint32_t h=hub_state[k];
      for(int direction=0;direction<2;++direction) {
        // forward fills Lin of descendants; backward fills Lout of ancestors.
        auto &mine=direction==0?lout[h]:lin[h];
        auto &target=direction==0?lin:lout;
        const uint32_t tag=2*k+direction;
        for(auto r:mine)hub_mark[r]=tag;
        queue.assign(1,h);seen[h]=tag;
        for(size_t i=0;i<queue.size();++i) {
          auto v=queue[i];++bfs_visits;
          if(v!=h) {
            bool covered=false;
            if(bp_hubs)covered=direction==0?(out_bp[h]&in_bp[v])!=0:(out_bp[v]&in_bp[h])!=0;
            for(auto r:target[v])if(covered||hub_mark[r]==tag){covered=true;break;}
            if(covered)continue;
          }
          target[v].push_back(k);++num_entries;
          auto visit=[&](uint32_t w){if(seen[w]!=tag && rank[w]>k){seen[w]=tag;queue.push_back(w);}};
          if(direction==0)trek_out(g,v,visit);else trek_in(g,v,visit);
          if(sc){for(const auto &e:(direction==0?sc->out:sc->in)[v])visit(e.to);}
        }
      }
    }
    auto flatten=[&](std::vector<std::vector<uint32_t>> &rows,std::vector<uint64_t> &off,std::vector<uint32_t> &flat,std::vector<uint8_t> &bytes) {
      off.resize(size_t(states)+1);
      if(!compress) {
        size_t total=0;
        for(uint32_t s=0;s<states;++s){off[s]=total;total+=rows[s].size();}
        off[states]=total;flat.resize(total);
        for(uint32_t s=0;s<states;++s){std::copy(rows[s].begin(),rows[s].end(),flat.begin()+off[s]);std::vector<uint32_t>().swap(rows[s]);}
      } else {
        for(uint32_t s=0;s<states;++s) {
          off[s]=bytes.size();uint32_t prev=0;
          for(auto r:rows[s]){uint32_t d=r-prev;prev=r;do{uint8_t b=d&127;d>>=7;bytes.push_back(uint8_t(b|(d?128:0)));}while(d);}
          std::vector<uint32_t>().swap(rows[s]);
        }
        off[states]=bytes.size();bytes.shrink_to_fit();
      }
    };
    compressed=compress;
    flatten(lout,out_off,out_hub,out_bytes);flatten(lin,in_off,in_hub,in_bytes);
    if(topo_filter) {
      use_topo=true;ord1.assign(states,0);ord2.assign(states,0);
      std::vector<uint32_t> indeg(states),base(states);
      for(uint32_t v=0;v<states;++v){uint32_t d=0;trek_in(g,v,[&](uint32_t){++d;});if(sc)d+=uint32_t(sc->in[v].size());base[v]=d;}
      for(int pass=0;pass<2;++pass) {
        indeg=base;std::vector<uint32_t> q;q.reserve(states);
        for(uint32_t v=0;v<states;++v)if(!indeg[v])q.push_back(v);
        uint32_t next=0;auto &ord=pass==0?ord1:ord2;
        if(pass==0){for(size_t i=0;i<q.size();++i){uint32_t v=q[i];ord[v]=next++;
            auto f=[&](uint32_t w){if(--indeg[w]==0)q.push_back(w);};trek_out(g,v,f);if(sc)for(const auto &e:sc->out[v])f(e.to);}}
        else {std::vector<uint32_t> st(q.rbegin(),q.rend());while(!st.empty()){uint32_t v=st.back();st.pop_back();ord[v]=next++;
            auto f=[&](uint32_t w){if(--indeg[w]==0)st.push_back(w);};trek_out(g,v,f);if(sc)for(const auto &e:sc->out[v])f(e.to);}}
      }
    }
    build_ms=elapsed(start);
  }
  // Sequential reader of one label row (sorted hub ranks), raw or compressed.
  struct Cursor {
    const uint32_t *p=nullptr,*e=nullptr;const uint8_t *b=nullptr,*be=nullptr;uint32_t prev=0;bool raw=true;
    bool next(uint32_t &r) {
      if(raw){if(p==e)return false;r=*p++;return true;}
      if(b==be)return false;
      uint32_t v=0;int sh=0;uint8_t x;do{x=*b++;v|=uint32_t(x&127)<<sh;sh+=7;}while(x&128);
      prev+=v;r=prev;return true;
    }
  };
  Cursor out(uint32_t s)const{return cursor(s,out_off,out_hub,out_bytes);}
  Cursor in(uint32_t s)const{return cursor(s,in_off,in_hub,in_bytes);}
  // Row length in storage units (entries, or bytes when compressed): the
  // adaptive engine's work unit for scanning this row.
  uint64_t out_len(uint32_t s)const{return out_off[s+1]-out_off[s];}
  uint64_t in_len(uint32_t s)const{return in_off[s+1]-in_off[s];}
  bool reach(uint32_t s,uint32_t t)const {
    if(s==t)return true;
    if(use_topo){++topo_tests;if(ord1[s]>=ord1[t]||ord2[s]>=ord2[t]){++topo_pruned;return false;}}
    if(bp_hubs&&(out_bp[s]&in_bp[t]))return true;
    auto a=out(s),b=in(t);uint32_t x,y;
    if(!a.next(x)||!b.next(y))return false;
    for(;;) {
      if(x==y)return true;
      if(x<y){if(!a.next(x))return false;}
      else if(!b.next(y))return false;
    }
  }
  size_t entries()const{return num_entries;}
  size_t bytes()const {
    return sizeof(*this)+4*(hub_state.capacity()+out_hub.capacity()+in_hub.capacity())+
           out_bytes.capacity()+in_bytes.capacity()+8*(out_off.capacity()+in_off.capacity()+out_bp.capacity()+in_bp.capacity());
  }
 private:
  Cursor cursor(uint32_t s,const std::vector<uint64_t> &off,const std::vector<uint32_t> &raw,const std::vector<uint8_t> &bytes)const {
    Cursor c;c.raw=!compressed;
    if(c.raw){c.p=raw.data()+off[s];c.e=raw.data()+off[s+1];}
    else {c.b=bytes.data()+off[s];c.be=bytes.data()+off[s+1];}
    return c;
  }
};

struct BudgetExceeded {};
struct ZSkeletonStats {
  uint64_t queries=0,h_nodes=0,edge_checks=0,neg_reach=0,cert_free=0,
           searches=0,search_found=0,search_states=0,label_merges=0,
           adaptive_index=0,adaptive_traversal=0,adaptive_rounds=0,
           fringe_resolved=0,core_routes=0,closures=0,candidates=0,probe_scans=0,proj_fast=0,proj_dfs=0,mask_builds=0,mask_entries=0;   // core–tree engine; mask cost
  // core–tree engine: H_Z edge tests by their costliest step (counters only; answers unchanged)
  uint64_t et_proj=0,et_refuted=0,et_certified=0,et_searched=0,et_closure=0;
};

// Graph surgery applied at query time: G with every edge INTO cut_in and every
// edge OUT OF cut_out removed (Pearl's G_{\overline{A}\underline{B}}). The
// index is never rebuilt. In T, cutting the in-edges of v disables the
// out-arcs of P_v and the in-arcs of E_v; cutting the out-edges of v disables
// the in-arcs of U_v and the out-arcs of D_v. A path uses no removed arc iff
// those states are not interior, P_v (cut_in) is not its start, and E_v
// (cut_in) is not its end. The surgered trek DAG T' is a subgraph of T, so
// reachability in T still gives exact negatives.
struct Surgery { std::vector<int> cut_in,cut_out; };

// Exact queries. Each H edge is a reachability test in T' where no interior
// state may belong to a Z vertex or be a surgery-blocked state. Interior-
// blocked states B = {U_f,E_f : f in Z} ∪ surgery states; P_f/D_f of Z
// vertices are reachable only through U_f/E_f.
// 1. not reachable in T at all            -> no edge (monotone, exact)
// 2. no blocked state between start, end  -> edge (exact)
// 3. otherwise exact bidirectional search in T minus forbidden states.
// (A hub certificate was evaluated and removed: it resolved <=1.2% of checks
// and never paid for itself.)
// Per query, every label hub rank r gets two bit masks: which blocked states
// have r in Lin (mask_in) or Lout (mask_out). One scan of a state's label then
// yields reachability to/from ALL blocked states at once.
struct ZSkeleton {
  friend struct CoreTreeZSkeleton;   // core–tree driver reuses the edge machinery on the core
  const Graph *g;
  const Shortcuts *sc;               // contracted arcs between core states, or nullptr
  TrekLabels labels;
  // Ablation switches (all on by default; answers stay exact either way).
  bool bidirectional_search=true,live_states=true,balance_by_arcs=true;
  explicit ZSkeleton(const Graph &graph,bool progress=false,const std::string &order="degree",bool compress=false,const Shortcuts *shortcuts=nullptr,int bp=0)
    :g(&graph),sc(shortcuts&&!shortcuts->empty()?shortcuts:nullptr),labels(graph,progress,order,compress,sc,bp) {}

  struct Scratch {
    QueryContext z;                        // Z membership (or avoided vertices)
    std::vector<uint32_t> cut_in,cut_out;  // surgery membership, epoch cut_epoch
    std::vector<int> cut_in_list,cut_out_list;
    std::vector<uint32_t> forbid;          // per state: in B, epoch forbid_epoch
    uint32_t cut_epoch=0,forbid_epoch=0;
    std::vector<uint32_t> visit,visit_back;// search epochs per state
    std::vector<uint32_t> mark_src,mark_dst;// epochs per hub rank: Lout(a), Lin(t)
    uint32_t visit_epoch=0,mark_epoch=0;
    // Sparse masks: slot_in[r] / slot_out[r] point into a pool of W-word
    // bitsets over the blocked states, valid while slot_epoch matches. Only
    // ranks in some blocked state's label get a slot, so memory is
    // O(states) + O(sum of blocked label lengths * W), never states * W.
    std::vector<uint32_t> slot_in,slot_out,slot_in_ep,slot_out_ep;
    std::vector<uint64_t> pool_in,pool_out;
    uint32_t slot_epoch=0;
    std::vector<int> nodes;                // H vertex ids
    std::vector<uint8_t> conditioned,side; // per H node: in Z?; 1=X side,2=Y side
    std::vector<uint32_t> blocked;         // mask-set states: interior-blocked, blocked-arc ends, and probes
    std::vector<uint8_t> blocked_cert;     // per mask-set entry: 1 = counts for the certificate, 0 = probe only
    std::vector<uint64_t> cert_mask;       // W words: certificate-relevant bits
    std::vector<uint32_t> bidx,bidx_ep;uint32_t bidx_epoch=0;   // state -> index in `blocked` (epoch-stamped)
    std::vector<uint8_t> row_ready,col_ready,reached;
    std::vector<uint64_t> row,col,rel,nb,tmp_row,tmp_col;
    std::vector<int> stack_nodes;
    std::vector<uint32_t> front[2],next;
    std::vector<uint64_t> front_bits[2],next_bits;
    size_t words=1;
    const uint32_t *arc_block=nullptr;uint32_t arc_epoch=0;   // blocked shortcut arcs (epoch-stamped), set by the core–tree driver
    uint64_t src_bp=0,dst_bp=0;                 // bit-parallel labels of the marked endpoints
    ZSkeletonStats *stats=nullptr;              // optional counters for mask builds (set per query)
    std::vector<uint64_t> blk_in_bp,blk_out_bp; // bit-parallel labels of the mask-set states
    bool masks_ready=false;              // masks are built at the first edge that passes the T test
    uint64_t work=0,budget=UINT64_MAX;   // label entries + arcs scanned
    explicit Scratch(int n):z(n),cut_in(n),cut_out(n),forbid(4*size_t(n)),visit(4*size_t(n)),visit_back(4*size_t(n)),
      mark_src(4*size_t(n)),mark_dst(4*size_t(n)),slot_in(4*size_t(n)),slot_out(4*size_t(n)),slot_in_ep(4*size_t(n)),slot_out_ep(4*size_t(n)),
      bidx(4*size_t(n)),bidx_ep(4*size_t(n)) {}
    // Fresh mask set for a query.
    void reset_blocked(){blocked.clear();blocked_cert.clear();if(++bidx_epoch==0){std::fill(bidx_ep.begin(),bidx_ep.end(),0);bidx_epoch=1;}}
    // Index of a state in the mask set, or UINT32_MAX.
    uint32_t mask_index(uint32_t sigma)const{return bidx_ep[sigma]==bidx_epoch?bidx[sigma]:UINT32_MAX;}
    bool in_cut(int v)const{return cut_in[v]==cut_epoch;}
    bool out_cut(int v)const{return cut_out[v]==cut_epoch;}
    size_t bytes()const {
      return sizeof(*this)+4*(z.mark.capacity()+cut_in.capacity()+cut_out.capacity()+forbid.capacity()+visit.capacity()+
             visit_back.capacity()+mark_src.capacity()+mark_dst.capacity()+blocked.capacity()+
             slot_in.capacity()+slot_out.capacity()+slot_in_ep.capacity()+slot_out_ep.capacity())+
             8*(pool_in.capacity()+pool_out.capacity()+row.capacity()+col.capacity());
    }
  };

  bool separated(int x,int y,const std::vector<int> &zs,Scratch &s,ZSkeletonStats *st=nullptr)const {
    return separated_general({x},{y},zs,Surgery{},s,st);
  }
  // Set endpoints: X and Y nonempty and disjoint, both disjoint from Z.
  bool separated_sets(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,Scratch &s,ZSkeletonStats *st=nullptr)const {
    return separated_general(xs,ys,zs,Surgery{},s,st);
  }
  // X ⊥ Y | Z in the surgered graph.
  bool separated_general(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut,Scratch &s,ZSkeletonStats *st=nullptr)const {
    mark_query(xs,ys,zs,cut,s);
    return answer_marked(xs,ys,zs,s,st);
  }
  // Validate and mark Z and the surgery. Both adaptive arms read these marks.
  void mark_query(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut,Scratch &s)const {
    mark_marks(*g,xs,ys,zs,cut,s);
  }
  // Index-free marking, shared with traversal-only baselines.
  static void mark_marks(const Graph &graph,const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut,Scratch &s) {
    const Graph *g=&graph;
    if(xs.empty()||ys.empty())throw std::invalid_argument("empty endpoint set");
    if(int(s.z.mark.size())!=g->n)throw std::invalid_argument("scratch size");
    if(++s.z.epoch==0){std::fill(s.z.mark.begin(),s.z.mark.end(),0);s.z.epoch=1;}
    for(int v:zs){if(v<0||v>=g->n)throw std::invalid_argument("conditioning vertex");s.z.mark[v]=s.z.epoch;}
    if(++s.forbid_epoch==0){std::fill(s.forbid.begin(),s.forbid.end(),0);s.forbid_epoch=1;}
    for(const auto *set:{&xs,&ys})for(int v:*set) {
      if(v<0||v>=g->n||s.z.contains(v)||s.forbid[4*size_t(v)]==s.forbid_epoch)
        throw std::invalid_argument("endpoint sets must be disjoint and outside Z");
      s.forbid[4*size_t(v)]=s.forbid_epoch;  // duplicate detection only
    }
    if(++s.cut_epoch==0){std::fill(s.cut_in.begin(),s.cut_in.end(),0);std::fill(s.cut_out.begin(),s.cut_out.end(),0);s.cut_epoch=1;}
    for(int v:cut.cut_in){if(v<0||v>=g->n)throw std::invalid_argument("surgery vertex");s.cut_in[v]=s.cut_epoch;}
    for(int v:cut.cut_out){if(v<0||v>=g->n)throw std::invalid_argument("surgery vertex");s.cut_out[v]=s.cut_epoch;}
    s.cut_in_list=cut.cut_in;s.cut_out_list=cut.cut_out;
  }
  // v is a proper descendant of u in G (one label merge).
  bool is_descendant(int u,int v)const {
    if(u<0||v<0||u>=g->n||v>=g->n)throw std::invalid_argument("vertex");
    return u!=v && labels.reach(state(u,kDown),state(v,kDown));
  }
  bool is_descendant(int u,int v,Scratch &)const{return is_descendant(u,v);}
  // Some w in targets has a directed path u -> ... -> w whose interior avoids
  // `avoid` (equivalently: u is an ancestor of w in G with in-edges of `avoid`
  // removed). u, targets and avoid must be pairwise disjoint.
  bool reaches_avoiding(int u,const std::vector<int> &targets,const std::vector<int> &avoid,Scratch &s,ZSkeletonStats *st=nullptr)const {
    if(targets.empty())return false;
    mark_query({u},targets,avoid,Surgery{},s);
    if(++s.forbid_epoch==0){std::fill(s.forbid.begin(),s.forbid.end(),0);s.forbid_epoch=1;}
    s.reset_blocked();
    for(int v:avoid)add_blocked(s,state(v,kUp)),add_blocked(s,state(v,kEnter));
    s.masks_ready=false;
    try {
      const uint32_t a=state(u,kDown);bool found=false;
      for(int w:targets){if(connect(s,a,state(w,kEnter),nullptr,nullptr,st)){found=true;break;}}
      clear_masks(s);return found;
    } catch(const BudgetExceeded &){clear_masks(s);throw;}
  }
  // Batched form for rule 3: out[i] = reaches_avoiding(us[i], targets, avoid),
  // with one marking pass and one mask build shared by all sources.
  void reaches_avoiding_many(const std::vector<int> &us,const std::vector<int> &targets,const std::vector<int> &avoid,
                             std::vector<uint8_t> &out,Scratch &s,ZSkeletonStats *st=nullptr)const {
    out.assign(us.size(),0);
    if(targets.empty()||us.empty())return;
    // mark_query needs disjoint endpoint sets; mark with the first source and
    // reuse the marks (they cover Z = avoid and no surgery).
    mark_query({us[0]},targets,avoid,Surgery{},s);
    if(++s.forbid_epoch==0){std::fill(s.forbid.begin(),s.forbid.end(),0);s.forbid_epoch=1;}
    s.reset_blocked();
    for(int v:avoid)add_blocked(s,state(v,kUp)),add_blocked(s,state(v,kEnter));
    s.masks_ready=false;
    try {
      for(size_t i=0;i<us.size();++i){const uint32_t a=state(us[i],kDown);
        for(int w:targets)if(connect(s,a,state(w,kEnter),nullptr,nullptr,st)){out[i]=1;break;}}
      clear_masks(s);
    } catch(const BudgetExceeded &){clear_masks(s);throw;}
  }
  // Answer with marks already set by mark_query (used by the adaptive engine).
  bool answer_marked(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,Scratch &s,ZSkeletonStats *st)const {
    if(st)++st->queries;
    // Fresh blocked-set epoch: an aborted adaptive round may have left marks.
    if(++s.forbid_epoch==0){std::fill(s.forbid.begin(),s.forbid.end(),0);s.forbid_epoch=1;}
    s.nodes.clear();s.conditioned.clear();s.side.clear();s.reset_blocked();
    for(int v:xs){s.nodes.push_back(v);s.conditioned.push_back(0);s.side.push_back(1);}
    for(int v:ys){s.nodes.push_back(v);s.conditioned.push_back(0);s.side.push_back(2);}
    for(int v:zs) {   // duplicate Z entries are legal input; keep each vertex once
      if(s.forbid[state(v,kUp)]==s.forbid_epoch)continue;
      s.nodes.push_back(v);s.conditioned.push_back(1);s.side.push_back(0);
      add_blocked(s,state(v,kUp));add_blocked(s,state(v,kEnter));
    }
    for(int v:s.cut_in_list){add_blocked(s,state(v,kLeave));add_blocked(s,state(v,kEnter));}
    for(int v:s.cut_out_list){add_blocked(s,state(v,kUp));add_blocked(s,state(v,kDown));}
    return run_h(s,st);
  }

 private:
  bool forbidden(const Scratch &s,uint32_t w)const{return s.z.contains(vertex_of(w))||s.forbid[w]==s.forbid_epoch;}
  // Mask-set entries. `cert` entries count for the certificate (forbidden
  // interior states, and the ends of blocked shortcut arcs); probe entries
  // only answer "is this state label-reachable from a" in one scan.
  uint32_t add_mask(Scratch &s,uint32_t sigma,bool cert)const {
    uint32_t i=s.mask_index(sigma);
    if(i==UINT32_MAX){i=uint32_t(s.blocked.size());s.blocked.push_back(sigma);s.blocked_cert.push_back(cert);s.bidx[sigma]=i;s.bidx_ep[sigma]=s.bidx_epoch;}
    else if(cert)s.blocked_cert[i]=1;
    return i;
  }
  void add_blocked(Scratch &s,uint32_t sigma)const{if(s.forbid[sigma]!=s.forbid_epoch){s.forbid[sigma]=s.forbid_epoch;add_mask(s,sigma,true);}}
  // The state joins the certificate's set (so that a blocked shortcut arc at it
  // is never silently assumed usable) but stays usable as an interior state.
  void add_mask_only(Scratch &s,uint32_t sigma)const{add_mask(s,sigma,true);}
  uint32_t add_probe(Scratch &s,uint32_t sigma)const{return add_mask(s,sigma,false);}
  bool arc_open(const Scratch &s,uint32_t arc)const{return !s.arc_block||s.arc_block[arc]!=s.arc_epoch;}
  size_t out_degree(uint32_t w)const{return trek_out_degree(*g,w)+(sc?sc->out[w].size():0);}
  size_t in_degree(uint32_t w)const{return trek_in_degree(*g,w)+(sc?sc->in[w].size():0);}
  bool reach(Scratch &s,uint32_t a,uint32_t b,ZSkeletonStats *st)const {
    if(st)++st->label_merges;
    charge(s,labels.out_len(a)+labels.in_len(b));
    return labels.reach(a,b);
  }
  static void charge(Scratch &s,uint64_t units){s.work+=units;if(s.work>s.budget)throw BudgetExceeded{};}

  // OR of per-rank masks over a state's outgoing (to blocked) or incoming
  // (from blocked) label; optionally also tests one marked endpoint label.
  bool scan_out(Scratch &s,uint32_t w,uint64_t *bits,bool test_dst,ZSkeletonStats *st)const {
    if(st)++st->label_merges;
    const size_t W=s.words;std::fill(bits,bits+W,0);bool hit=false;
    charge(s,labels.out_len(w));
    auto c=labels.out(w);
    for(uint32_t r;c.next(r);) {
      if(test_dst && s.mark_dst[r]==s.mark_epoch)hit=true;
      if(s.slot_in_ep[r]!=s.slot_epoch)continue;
      const uint64_t *m=&s.pool_in[size_t(s.slot_in[r])*W];for(size_t q=0;q<W;++q)bits[q]|=m[q];
    }
    if(labels.bp_hubs) {   // via a bit-parallel hub: w reaches hub k, hub k reaches the target
      const uint64_t ob=labels.out_bp[w];
      if(ob) {
        if(test_dst&&(ob&s.dst_bp))hit=true;
        for(size_t b=0;b<s.blk_in_bp.size();++b)if(ob&s.blk_in_bp[b])bits[b/64]|=uint64_t(1)<<(b%64);
      }
    }
    return hit;
  }
  bool scan_in(Scratch &s,uint32_t w,uint64_t *bits,bool test_src,ZSkeletonStats *st)const {
    if(st)++st->label_merges;
    const size_t W=s.words;std::fill(bits,bits+W,0);bool hit=false;
    charge(s,labels.in_len(w));
    auto c=labels.in(w);
    for(uint32_t r;c.next(r);) {
      if(test_src && s.mark_src[r]==s.mark_epoch)hit=true;
      if(s.slot_out_ep[r]!=s.slot_epoch)continue;
      const uint64_t *m=&s.pool_out[size_t(s.slot_out[r])*W];for(size_t q=0;q<W;++q)bits[q]|=m[q];
    }
    if(labels.bp_hubs) {
      const uint64_t ib=labels.in_bp[w];
      if(ib) {
        if(test_src&&(ib&s.src_bp))hit=true;
        for(size_t b=0;b<s.blk_out_bp.size();++b)if(ib&s.blk_out_bp[b])bits[b/64]|=uint64_t(1)<<(b%64);
      }
    }
    return hit;
  }
  // Charge the mask cost before building them, so an over-budget round
  // aborts without paying it.
  void begin_masks(Scratch &s)const {
    s.words=std::max<size_t>(1,(s.blocked.size()+63)/64);
    uint64_t cost=0;
    for(auto sigma:s.blocked)cost+=labels.out_len(sigma)+labels.in_len(sigma);
    charge(s,cost);
    if(s.stats){++s.stats->mask_builds;s.stats->mask_entries+=cost;}
    const size_t W=s.words;
    if(++s.slot_epoch==0){std::fill(s.slot_in_ep.begin(),s.slot_in_ep.end(),0);std::fill(s.slot_out_ep.begin(),s.slot_out_ep.end(),0);s.slot_epoch=1;}
    s.pool_in.clear();s.pool_out.clear();
    auto slot=[&](std::vector<uint32_t> &sl,std::vector<uint32_t> &ep,std::vector<uint64_t> &pool,uint32_t r)->uint64_t* {
      if(ep[r]!=s.slot_epoch){ep[r]=s.slot_epoch;sl[r]=uint32_t(pool.size()/W);pool.resize(pool.size()+W,0);}
      return &pool[size_t(sl[r])*W];
    };
    s.blk_in_bp.clear();s.blk_out_bp.clear();
    if(labels.bp_hubs)for(auto sigma:s.blocked){s.blk_in_bp.push_back(labels.in_bp[sigma]);s.blk_out_bp.push_back(labels.out_bp[sigma]);}
    for(size_t b=0;b<s.blocked.size();++b) {
      auto sigma=s.blocked[b];const uint64_t bit=uint64_t(1)<<(b%64);const size_t q=b/64;
      {auto c=labels.in(sigma);for(uint32_t r;c.next(r);)slot(s.slot_in,s.slot_in_ep,s.pool_in,r)[q]|=bit;}
      {auto c=labels.out(sigma);for(uint32_t r;c.next(r);)slot(s.slot_out,s.slot_out_ep,s.pool_out,r)[q]|=bit;}
    }
  }
  // First edge that needs the blocked states: build masks now.
  void ensure_masks(Scratch &s)const {
    if(s.masks_ready)return;
    begin_masks(s);
    s.cert_mask.assign(s.words,0);
    for(size_t b=0;b<s.blocked.size();++b)if(s.blocked_cert[b])s.cert_mask[b/64]|=uint64_t(1)<<(b%64);
    const size_t h=s.nodes.size();
    s.row.assign(h*s.words,0);s.col.assign(h*s.words,0);s.masks_ready=true;
  }
  // Masks are invalidated by the next begin_masks (epoch); nothing to clear.
  void clear_masks(Scratch &)const {}
  // Blocked states reachable from a / reaching t; a path's own start and end
  // are never "between" themselves.
  void row_bits(Scratch &s,uint32_t a,uint64_t *r,ZSkeletonStats *st)const {
    scan_out(s,a,r,false,st);
    for(size_t b=0;b<s.blocked.size();++b)if(s.blocked[b]==a)r[b/64]&=~(uint64_t(1)<<(b%64));
  }
  void col_bits(Scratch &s,uint32_t t,uint64_t *c,ZSkeletonStats *st)const {
    scan_in(s,t,c,false,st);
    for(size_t b=0;b<s.blocked.size();++b)if(s.blocked[b]==t)c[b/64]&=~(uint64_t(1)<<(b%64));
  }
  void ensure_row(Scratch &s,int i,ZSkeletonStats *st)const {
    if(s.row_ready[i])return;
    s.row_ready[i]=1;row_bits(s,live_start(s,i),&s.row[size_t(i)*s.words],st);
  }
  void ensure_col(Scratch &s,int j,ZSkeletonStats *st)const {
    if(s.col_ready[j])return;
    s.col_ready[j]=1;col_bits(s,live_end(s,j),&s.col[size_t(j)*s.words],st);
  }
  void mark_endpoints(Scratch &s,uint32_t a,uint32_t t)const {
    charge(s,labels.out_len(a)+labels.in_len(t));
    if(++s.mark_epoch==0){std::fill(s.mark_src.begin(),s.mark_src.end(),0);std::fill(s.mark_dst.begin(),s.mark_dst.end(),0);s.mark_epoch=1;}
    {auto c=labels.out(a);for(uint32_t r;c.next(r);)s.mark_src[r]=s.mark_epoch;}
    {auto c=labels.in(t);for(uint32_t r;c.next(r);)s.mark_dst[r]=s.mark_epoch;}
    if(labels.bp_hubs){s.src_bp=labels.out_bp[a];s.dst_bp=labels.in_bp[t];}
  }

  // Under surgery an endpoint's start/end may keep a single live arc. Every
  // T' path then passes that neighbour, so test from it: the T-reachability
  // filter becomes much sharper. U_v -> {P_v, D_v}: D_v is dead if v's out-
  // edges are cut, P_v if its in-edges are cut. D_v <- {U_v, E_v}: U_v is
  // unenterable if out-edges are cut, E_v if in-edges are cut.
  // Returns UINT32_MAX when nothing survives.
  uint32_t live_start(const Scratch &s,int i)const {
    int v=s.nodes[i];
    if(s.conditioned[i])return s.in_cut(v)?UINT32_MAX:state(v,kLeave);
    bool in=s.in_cut(v),out=s.out_cut(v);
    if(!live_states)return state(v,kUp);
    if(in&&out)return UINT32_MAX;
    return out?state(v,kLeave):in?state(v,kDown):state(v,kUp);
  }
  uint32_t live_end(const Scratch &s,int j)const {
    int v=s.nodes[j];
    if(s.conditioned[j])return s.in_cut(v)?UINT32_MAX:state(v,kEnter);
    bool in=s.in_cut(v),out=s.out_cut(v);
    if(!live_states)return state(v,kDown);
    if(in&&out)return UINT32_MAX;
    return out?state(v,kEnter):in?state(v,kUp):state(v,kDown);
  }
  bool edge(Scratch &s,int i,int j,ZSkeletonStats *st)const {
    if(st)++st->edge_checks;
    const uint32_t a=live_start(s,i),t=live_end(s,j);
    if(a==UINT32_MAX||t==UINT32_MAX){if(st)++st->neg_reach;return false;}
    if(!reach(s,a,t,st)){if(st)++st->neg_reach;return false;}
    ensure_masks(s);
    ensure_row(s,i,st);ensure_col(s,j,st);
    return connect(s,a,t,&s.row[size_t(i)*s.words],&s.col[size_t(j)*s.words],st,true);
  }
  // Reachability a -> t in T' given the configured blocked states. Row/column
  // bits are computed here when not supplied. `checked` skips the T test.
  bool connect(Scratch &s,uint32_t a,uint32_t t,const uint64_t *r,const uint64_t *c,ZSkeletonStats *st,bool checked=false)const {
    if(!checked) {
      if(st)++st->edge_checks;
      if(!reach(s,a,t,st)){if(st)++st->neg_reach;return false;}
    }
    if(!r||!c)ensure_masks(s);
    if(!r){s.tmp_row.resize(s.words);row_bits(s,a,s.tmp_row.data(),st);r=s.tmp_row.data();}
    if(!c){s.tmp_col.resize(s.words);col_bits(s,t,s.tmp_col.data(),st);c=s.tmp_col.data();}
    const size_t W=s.words;
    s.rel.assign(W,0);bool any=false;
    for(size_t q=0;q<W;++q){s.rel[q]=r[q]&c[q]&s.cert_mask[q];any|=s.rel[q]!=0;}
    // Blocked shortcut arcs have both ends in the mask set, so any one on an
    // a–t path puts a state strictly between a and t — except an arc a -> t
    // itself, which the certificate must not overlook.
    if(!any&&sc&&s.arc_block)for(const auto &e:sc->out[a])if(e.to==t&&!arc_open(s,e.arc)){any=true;break;}
    if(!any){if(st)++st->cert_free;return true;}
    if(st)++st->searches;
    bool found=search(s,a,t,st);
    if(st&&found)++st->search_found;
    return found;
  }

  // Bidirectional exact search in T minus every forbidden interior state.
  // Always expand the smaller frontier. Forward states must still reach t and
  // backward states must still be reachable from a (label pruning). Each side
  // carries the subset of relevant blocked states still between it and the
  // other end; an empty subset certifies a blocked-free completion.
  bool search(Scratch &s,uint32_t a,uint32_t t,ZSkeletonStats *st)const {
    if(++s.visit_epoch==0){std::fill(s.visit.begin(),s.visit.end(),0);std::fill(s.visit_back.begin(),s.visit_back.end(),0);s.visit_epoch=1;}
    mark_endpoints(s,a,t);
    const uint32_t ep=s.visit_epoch;const size_t W=s.words;
    auto &fs=s.front[0],&bs=s.front[1];auto &fb=s.front_bits[0],&bb=s.front_bits[1];
    fs.assign(1,a);bs.assign(1,t);fb=s.rel;bb=s.rel;
    s.visit[a]=ep;s.visit_back[t]=ep;s.nb.resize(W);
    // Expand the side whose next level scans fewer arcs (total out-degree of
    // the forward frontier vs total in-degree of the backward frontier).
    uint64_t fcost=out_degree(a),bcost=in_degree(t);
    while(!fs.empty()&&!bs.empty()) {
      const bool fwd=!bidirectional_search||(balance_by_arcs?fcost<=bcost:fs.size()<=bs.size());
      auto &cur=fwd?fs:bs;auto &cur_bits=fwd?fb:bb;
      auto &mine=fwd?s.visit:s.visit_back;auto &other=fwd?s.visit_back:s.visit;
      s.next.clear();s.next_bits.clear();bool found=false;
      for(size_t i=0;i<cur.size()&&!found;++i) {
        if(st)++st->search_states;
        const size_t base=i*W;
        auto expand=[&](uint32_t w) {
          if(found)return;
          charge(s,1);
          if(w!=a&&w!=t&&forbidden(s,w))return;
          if(other[w]==ep){found=true;return;}
          if(mine[w]==ep)return;
          mine[w]=ep;
          bool ok=fwd?scan_out(s,w,s.nb.data(),true,st):scan_in(s,w,s.nb.data(),true,st);
          if(!ok)return;
          bool any=false;
          for(size_t q=0;q<W;++q){s.nb[q]&=cur_bits[base+q];any|=s.nb[q]!=0;}
          if(!any){found=true;return;}
          s.next.push_back(w);s.next_bits.insert(s.next_bits.end(),s.nb.begin(),s.nb.end());
        };
        if(fwd)trek_out(*g,cur[i],expand);else trek_in(*g,cur[i],expand);
        if(sc){for(const auto &e:(fwd?sc->out:sc->in)[cur[i]]){if(found)break;if(arc_open(s,e.arc))expand(e.to);}}
      }
      if(found)return true;
      cur.swap(s.next);cur_bits.swap(s.next_bits);
      uint64_t c=0;for(auto w:cur)c+=fwd?out_degree(w):in_degree(w);
      (fwd?fcost:bcost)=c;
    }
    return false;
  }

  bool run_h(Scratch &s,ZSkeletonStats *st)const {
    const int h=int(s.nodes.size());
    if(st)st->h_nodes+=h;
    // Masks are deferred to the first edge that passes the T-reachability test:
    // queries decided by negatives alone never pay for them.
    s.masks_ready=false;
    try {
    s.row_ready.assign(h,0);s.col_ready.assign(h,0);s.reached.assign(h,0);
    s.stack_nodes.clear();
    for(int i=0;i<h;++i)if(s.side[i]==1){s.reached[i]=1;s.stack_nodes.push_back(i);}
    bool separated=true;
    while(!s.stack_nodes.empty() && separated) {
      int i=s.stack_nodes.back();s.stack_nodes.pop_back();
      // Try Y-side nodes first so connected answers stop early.
      for(int pass=0;pass<2 && separated;++pass)for(int j=0;j<h && separated;++j) {
        if(s.reached[j]||(pass==0)!=(s.side[j]==2))continue;
        if(edge(s,i,j,st)) {
          if(s.side[j]==2)separated=false;
          s.reached[j]=1;s.stack_nodes.push_back(j);
        }
      }
    }
    clear_masks(s);
    return separated;
    } catch(const BudgetExceeded &) {clear_masks(s);throw;}
  }
};

// Ablation of the index: the trek-DAG reduction with plain searches in T, no
// labels, no certificate, no decomposition. H_Z is explored as in
// ZSkeleton::run_h (Y-side nodes first), and each edge test is a bidirectional
// search in T whose interior avoids every state of a Z vertex ("trek-edge").
// In mode Sweep ("trek-sweep") the searches of a query form one forward
// traversal of T instead: the search from each newly reached H node collects
// every end state it meets and skips states that an earlier search of the
// same query visited, whose end states are already known. In mode Bidir
// ("trek-bidir") one bidirectional search runs from U_x to D_y in T without
// the states of Z vertices, plus an arc E_z -> P_z for each z in Z (a path
// that enters z from a parent leaves it upward: a collider at z); this is the
// path concatenation of H_Z, i.e. Bayes-Ball on the four states of T.
struct TrekSearch {
  enum Mode { Edge, Sweep, Bidir };
  const Graph *g;
  Mode mode=Edge;
  QueryContext z;
  std::vector<uint32_t> visit,visit_back,owner,owner_ep;
  uint32_t epoch=0,query_epoch=0;
  std::vector<int> nodes,stack_nodes;
  std::vector<uint8_t> conditioned,side,reached;
  std::vector<uint32_t> fs,bs,next;
  uint64_t edge_tests=0,states=0;
  TrekSearch(const Graph &graph,Mode m):g(&graph),mode(m),z(graph.n),visit(4*size_t(graph.n)),
    visit_back(4*size_t(graph.n)),owner(4*size_t(graph.n)),owner_ep(4*size_t(graph.n)) {}
  static uint32_t start_state(int v,bool cond){return state(v,cond?kLeave:kUp);}
  static uint32_t end_state(int v,bool cond){return state(v,cond?kEnter:kDown);}
  bool separated(int x,int y,const std::vector<int> &zs) {
    if(++z.epoch==0){std::fill(z.mark.begin(),z.mark.end(),0);z.epoch=1;}
    for(int v:zs){if(v<0||v>=g->n)throw std::invalid_argument("conditioning vertex");z.mark[v]=z.epoch;}
    if(x<0||y<0||x>=g->n||y>=g->n||x==y||z.contains(x)||z.contains(y))throw std::invalid_argument("endpoints");
    if(mode==Bidir)return !bidir(x,y);
    if(++query_epoch==0){std::fill(owner_ep.begin(),owner_ep.end(),0);query_epoch=1;}
    nodes.clear();conditioned.clear();side.clear();
    auto add=[&](int v,uint8_t cond,uint8_t sd) {   // each vertex once; owner maps its end state to its H node
      uint32_t t=end_state(v,cond);
      if(owner_ep[t]==query_epoch)return;
      owner_ep[t]=query_epoch;owner[t]=uint32_t(nodes.size());
      nodes.push_back(v);conditioned.push_back(cond);side.push_back(sd);
    };
    add(x,0,1);add(y,0,2);for(int v:zs)add(v,1,0);
    const int h=int(nodes.size());
    reached.assign(h,0);stack_nodes.clear();
    for(int i=0;i<h;++i)if(side[i]==1){reached[i]=1;stack_nodes.push_back(i);}
    if(mode==Sweep)return sweep();
    while(!stack_nodes.empty()) {
      int i=stack_nodes.back();stack_nodes.pop_back();
      for(int pass=0;pass<2;++pass)for(int j=0;j<h;++j) {
        if(reached[j]||(pass==0)!=(side[j]==2))continue;
        if(search(start_state(nodes[i],conditioned[i]),end_state(nodes[j],conditioned[j]))) {
          if(side[j]==2)return false;
          reached[j]=1;stack_nodes.push_back(j);
        }
      }
    }
    return true;
  }
  size_t bytes()const{return sizeof(*this)+4*(z.mark.capacity()+visit.capacity()+visit_back.capacity()+owner.capacity()+owner_ep.capacity());}
 private:
  // a -> t in T with no interior state of a Z vertex; the smaller frontier (in arcs) is expanded.
  bool search(uint32_t a,uint32_t t) {
    ++edge_tests;
    if(++epoch==0){std::fill(visit.begin(),visit.end(),0);std::fill(visit_back.begin(),visit_back.end(),0);epoch=1;}
    fs.assign(1,a);bs.assign(1,t);visit[a]=epoch;visit_back[t]=epoch;
    uint64_t fcost=trek_out_degree(*g,a),bcost=trek_in_degree(*g,t);
    while(!fs.empty()&&!bs.empty()) {
      const bool fwd=fcost<=bcost;
      auto &cur=fwd?fs:bs;auto &mine=fwd?visit:visit_back;auto &other=fwd?visit_back:visit;
      next.clear();bool found=false;uint64_t c=0;
      for(size_t i=0;i<cur.size()&&!found;++i) {
        ++states;
        auto expand=[&](uint32_t w) {
          if(found)return;
          if(w!=a&&w!=t&&z.contains(vertex_of(w)))return;
          if(other[w]==epoch){found=true;return;}
          if(mine[w]==epoch)return;
          mine[w]=epoch;next.push_back(w);c+=fwd?trek_out_degree(*g,w):trek_in_degree(*g,w);
        };
        if(fwd)trek_out(*g,cur[i],expand);else trek_in(*g,cur[i],expand);
      }
      if(found)return true;
      cur.swap(next);(fwd?fcost:bcost)=c;
    }
    return false;
  }
  // U_x reaches D_y in T_Z; the smaller frontier (in arcs) is expanded.
  bool bidir(int x,int y) {
    ++edge_tests;
    if(++epoch==0){std::fill(visit.begin(),visit.end(),0);std::fill(visit_back.begin(),visit_back.end(),0);epoch=1;}
    const uint32_t a=state(x,kUp),t=state(y,kDown);
    fs.assign(1,a);bs.assign(1,t);visit[a]=epoch;visit_back[t]=epoch;
    uint64_t fcost=trek_out_degree(*g,a),bcost=trek_in_degree(*g,t);
    auto in_z=[&](uint32_t w){return z.contains(vertex_of(w));};
    while(!fs.empty()&&!bs.empty()) {
      const bool fwd=fcost<=bcost;
      auto &cur=fwd?fs:bs;auto &mine=fwd?visit:visit_back;auto &other=fwd?visit_back:visit;
      next.clear();bool found=false;uint64_t c=0;
      auto add=[&](uint32_t w) {
        if(found)return;
        if(other[w]==epoch){found=true;return;}
        if(mine[w]==epoch)return;
        mine[w]=epoch;next.push_back(w);c+=fwd?trek_out_degree(*g,w):trek_in_degree(*g,w);
      };
      for(size_t i=0;i<cur.size()&&!found;++i) {
        ++states;
        const uint32_t s=cur[i];
        // Of a Z vertex, only E_z (entered from a parent) and P_z (left upward) are used, joined by E_z -> P_z.
        if(fwd) {
          if(in_z(s)&&(s&3)==kEnter)add(state(vertex_of(s),kLeave));
          else trek_out(*g,s,[&](uint32_t w){if(!in_z(w)||(w&3)==kEnter)add(w);});
        } else {
          if(in_z(s)&&(s&3)==kLeave)add(state(vertex_of(s),kEnter));
          else trek_in(*g,s,[&](uint32_t w){if(!in_z(w)||(w&3)==kLeave)add(w);});
        }
      }
      if(found)return true;
      cur.swap(next);(fwd?fcost:bcost)=c;
    }
    return false;
  }
  // One traversal of T per query; the end state of an unreached H node makes it reached.
  bool sweep() {
    if(++epoch==0){std::fill(visit.begin(),visit.end(),0);epoch=1;}
    while(!stack_nodes.empty()) {
      int i=stack_nodes.back();stack_nodes.pop_back();
      const uint32_t a=start_state(nodes[i],conditioned[i]);
      if(visit[a]==epoch)continue;
      ++edge_tests;
      visit[a]=epoch;fs.assign(1,a);
      while(!fs.empty()) {
        uint32_t u=fs.back();fs.pop_back();++states;
        bool connected=false;
        trek_out(*g,u,[&](uint32_t w) {
          if(connected)return;
          if(owner_ep[w]==query_epoch) {
            uint32_t j=owner[w];
            if(!reached[j]){if(side[j]==2){connected=true;return;}reached[j]=1;stack_nodes.push_back(int(j));}
          }
          if(z.contains(vertex_of(w))||visit[w]==epoch)return;
          visit[w]=epoch;fs.push_back(w);
        });
        if(connected)return false;
      }
    }
    return true;
  }
};

// Native bidirectional local-rule traversal (same rules as BidirectionalCIfly),
// for endpoint sets and surgery, with a work budget in arcs scanned. Reads Z
// and surgery marks from a ZSkeleton::Scratch. Exact when it finishes.
struct BudgetedTraversal {
  const Graph *g;
  std::vector<uint32_t> fseen,bseen;
  std::vector<int> forward,backward,next;
  uint32_t epoch=0;
  explicit BudgetedTraversal(const Graph &gg):g(&gg),fseen(2*size_t(gg.n)),bseen(2*size_t(gg.n)) {}
  // Returns 0 separated, 1 connected, 2 budget exhausted.
  // The budget is checked on every arc, so a hub expansion cannot overshoot.
  int run(const std::vector<int> &xs,const std::vector<int> &ys,const ZSkeleton::Scratch &s,uint64_t budget) {
    if(++epoch==0){std::fill(fseen.begin(),fseen.end(),0);std::fill(bseen.begin(),bseen.end(),0);epoch=1;}
    const QueryContext &zc=s.z;
    auto keep=[&](int p,int c){return !s.in_cut(c)&&!s.out_cut(p);};   // edge p -> c survives
    uint64_t work=0;
    forward.clear();backward.clear();
    for(int x:xs){forward.push_back(2*x);fseen[2*x]=epoch;}
    for(int y:ys){backward.push_back(2*y);backward.push_back(2*y+1);bseen[2*y]=bseen[2*y+1]=epoch;}
    while(!forward.empty()&&!backward.empty()) {
      const bool f=forward.size()<=backward.size();
      auto &front=f?forward:backward;auto &seen=f?fseen:bseen;auto &other=f?bseen:fseen;
      next.clear();
      int status=-1;   // -1 continue, 1 met, 2 over budget
      auto push=[&](int st){
        if(++work>budget){status=2;return;}
        if(other[st]==epoch){status=1;return;}
        if(seen[st]!=epoch){seen[st]=epoch;next.push_back(st);}
      };
      for(size_t i=0;i<front.size()&&status<0;++i) {
        int v=front[i]/2,d=front[i]%2;
        const auto &pa=g->pa[v],&ch=g->ch[v];
        if(f) {
          if(!zc.contains(v)){
            if(!d)for(size_t k=0;k<pa.size()&&status<0;++k){if(keep(pa[k],v))push(2*pa[k]);}
            for(size_t k=0;k<ch.size()&&status<0;++k){if(keep(v,ch[k]))push(2*ch[k]+1);}
          } else if(d)for(size_t k=0;k<pa.size()&&status<0;++k){if(keep(pa[k],v))push(2*pa[k]);}
        } else if(!d){for(size_t k=0;k<ch.size()&&status<0;++k){if(keep(v,ch[k]))push(2*ch[k]+int(zc.contains(ch[k])));}}
        else for(size_t k=0;k<pa.size()&&status<0;++k)if(!zc.contains(pa[k])&&keep(pa[k],v)){push(2*pa[k]);if(status<0)push(2*pa[k]+1);}
      }
      if(status>0)return status;
      front.swap(next);
    }
    return 0;
  }
};

// Exact adaptive engine: a tiny traversal probe, then a doubling race. Round r
// gives the index and then traversal B0*2^r work units each, restarting both.
// Both arms are exact. A query costs O(min(index work, traversal work)) plus
// the probe, up to a constant factor.
struct AdaptiveZSkeleton {
  const ZSkeleton *index;
  BudgetedTraversal traversal;
  uint64_t probe_budget=256,initial_budget=65536;
  bool probe_per_z=false;
  uint64_t large_h=64;    // above this |H_Z|, the probe grows to h(h-1)/2
  AdaptiveZSkeleton(const ZSkeleton &idx,const Graph &g):index(&idx),traversal(g) {}
  bool separated(int x,int y,const std::vector<int> &zs,ZSkeleton::Scratch &s,ZSkeletonStats *st=nullptr) {
    return separated_general({x},{y},zs,Surgery{},s,st);
  }
  bool separated_general(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut,
                         ZSkeleton::Scratch &s,ZSkeletonStats *st=nullptr) {
    // Validate and mark Z and surgery once; both arms read the same marks.
    index->mark_query(xs,ys,zs,cut,s);
    // When Z hems x in, traversal ends within a few arcs. The index must
    // consider up to h(h-1)/2 H_Z pairs (h = |X u Y u Z|), a lower bound on
    // its work, so the probe is at least that: large-Z queries go straight to
    // traversal, small-Z queries keep the small probe.
    const uint64_t h=xs.size()+ys.size()+zs.size();
    uint64_t probe=h>large_h?std::max<uint64_t>(probe_budget,h*(h-1)/2):probe_budget;
    if(probe_per_z)probe=std::max<uint64_t>(probe,probe_budget*zs.size());
    int r=traversal.run(xs,ys,s,probe);
    if(r!=2){if(st)++st->adaptive_traversal;return r==0;}
    for(uint64_t budget=initial_budget;;budget*=2) {
      if(st)++st->adaptive_rounds;
      s.work=0;s.budget=budget;
      try {
        bool answer=index->answer_marked(xs,ys,zs,s,st);
        s.budget=UINT64_MAX;if(st)++st->adaptive_index;return answer;
      } catch(const BudgetExceeded &) {}
      s.budget=UINT64_MAX;
      r=traversal.run(xs,ys,s,budget);
      if(r!=2){if(st)++st->adaptive_traversal;return r==0;}
    }
  }
};
} // namespace rs
