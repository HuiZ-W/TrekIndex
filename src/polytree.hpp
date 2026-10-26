#pragma once
#include "universal.hpp"
namespace rs {
struct PolyScratch {
  std::vector<uint32_t> activated;
  uint32_t epoch = 0;
  explicit PolyScratch(int n) : activated(n) {}
  void next() {
    if (++epoch == 0) {
      std::fill(activated.begin(), activated.end(), 0);
      epoch = 1;
    }
  }
};
// Exact graph-only oracle for DAGs whose undirected skeleton is a forest.
// The constructor rejects all other DAGs, rather than approximating them.
struct PolytreeIndex {
  int n, levels = 1;
  bool collider_free = true;
  std::vector<int> parent, depth, component, tin, tout, down_prefix, collider_prefix;
  std::vector<uint8_t> down, weight;
  std::vector<std::vector<int>> up;
  explicit PolytreeIndex(const Graph &g)
      : n(g.n), parent(n,-1), depth(n), component(n,-1), tin(n), tout(n),
        down_prefix(n), collider_prefix(n), down(n), weight(n) {
    for (int v=0; v<n; ++v) if (g.pa[v].size()>1) collider_free=false;
    int clock=0, max_depth=0;
    for (int root : g.order) if (component[root]<0) {
      component[root]=root; parent[root]=root; tin[root]=clock++;
      std::vector<std::pair<int,size_t>> stack{{root,0}};
      while (!stack.empty()) {
        int v=stack.back().first;
        size_t &i=stack.back().second;
        if (i==g.pa[v].size()+g.ch[v].size()) {
          tout[v]=clock; stack.pop_back(); continue;
        }
        bool forward=i>=g.pa[v].size();
        int w=forward ? g.ch[v][i-g.pa[v].size()] : g.pa[v][i];
        ++i;
        if (w==parent[v]) continue;
        if (component[w]>=0) throw std::invalid_argument("polytree requires an acyclic undirected skeleton");
        parent[w]=v; component[w]=root; depth[w]=depth[v]+1;
        max_depth=std::max(max_depth,depth[w]);
        down[w]=forward;
        down_prefix[w]=down_prefix[v]+down[w];
        // Weight on tree edge (v,w) counts a collider at v using its own parent
        // and w. Corrections at the query LCA replace that parent by the other arm.
        weight[w]=down[v] && !down[w];
        collider_prefix[w]=collider_prefix[v]+weight[w];
        tin[w]=clock++; stack.emplace_back(w,0);
      }
    }
    while ((uint64_t(1)<<levels)<=uint64_t(max_depth)) ++levels;
    up.assign(levels,std::vector<int>(n));
    if (n) up[0]=parent;
    for (int j=1;j<levels;++j)
      for (int v=0;v<n;++v) up[j][v]=up[j-1][up[j-1][v]];
  }
  bool ancestor(int u,int v) const { return tin[u]<=tin[v] && tout[v]<=tout[u]; }
  int lift(int v,int distance) const {
    for (int j=0;distance;distance>>=1,++j) if (distance&1) v=up[j][v];
    return v;
  }
  int lca(int x,int y) const {
    if (component[x]!=component[y]) return -1;
    if (ancestor(x,y)) return x;
    if (ancestor(y,x)) return y;
    for (int j=levels-1;j>=0;--j) if (!ancestor(up[j][x],y)) x=up[j][x];
    return parent[x];
  }
  int arm(int l,int x) const { return l==x ? -1 : lift(x,depth[x]-depth[l]-1); }
  bool directed_reach(int x,int y) const {
    int l=lca(x,y);
    return l>=0 && down_prefix[x]==down_prefix[l] &&
           down_prefix[y]-down_prefix[l]==depth[y]-depth[l];
  }
  bool on_path(int z,int x,int y,int l) const {
    return depth[z]>=depth[l] && (ancestor(z,x)||ancestor(z,y));
  }
  int count_colliders(int x,int y,int l,int cx,int cy) const {
    int c=collider_prefix[x]+collider_prefix[y]-2*collider_prefix[l];
    if (cx>=0) c-=weight[cx];
    if (cy>=0) c-=weight[cy];
    if (cx>=0 && cy>=0 && !down[cx] && !down[cy]) ++c;
    return c;
  }
  bool collider_on_path(int p,int x,int y,int l,int cx,int cy) const {
    if (p==x || p==y) return false;
    if (p==l) return cx>=0 && cy>=0 && !down[cx] && !down[cy];
    if (!down[p]) return false;
    int endpoint=ancestor(p,x) ? x : y;
    return !down[arm(p,endpoint)];
  }
  bool separated(int x,int y,const std::vector<int> &z,PolyScratch &scratch) const {
    if (x<0 || y<0 || x>=n || y>=n || x==y || int(scratch.activated.size())!=n)
      throw std::invalid_argument("endpoints or scratch size");
    // Read/validate every supplied Z entry, even when a structural shortcut applies.
    for (int v:z) if (v<0 || v>=n || v==x || v==y)
      throw std::invalid_argument("conditioning vertex");
    if (component[x]!=component[y]) return true;
    int l=lca(x,y),cx=-1,cy=-1,colliders=0;
    if (!collider_free) {
      cx=arm(l,x);cy=arm(l,y);
      colliders=count_colliders(x,y,l,cx,cy);
    }
    if (colliders==0) {
      for (int v:z) if (on_path(v,x,y,l)) return true;
      return false;
    }
    // In a polytree, a single Z vertex can activate at most one path collider.
    if (size_t(colliders)>z.size()) return true;
    scratch.next(); int active=0;
    for (int v:z) {
      if (component[v]!=component[x]) continue;
      int p=l,a=lca(x,v),b=lca(y,v);
      if (depth[a]>depth[p]) p=a;
      if (depth[b]>depth[p]) p=b;
      bool collider=collider_on_path(p,x,y,l,cx,cy);
      if (p==v && !collider) return true;
      if (collider && scratch.activated[p]!=scratch.epoch && directed_reach(p,v)) {
        scratch.activated[p]=scratch.epoch; ++active;
      }
    }
    return active!=colliders;
  }
  size_t bytes() const {
    size_t b=sizeof(*this)+down.capacity()+weight.capacity();
    for (auto *v:{&parent,&depth,&component,&tin,&tout,&down_prefix,&collider_prefix}) b+=4*v->capacity();
    b+=up.capacity()*sizeof(std::vector<int>);
    for (auto &v:up) b+=4*v.capacity();
    return b;
  }
};

// Exact bidirectional BFS on the implicit CIfly-like guarded state graph.
// Backward search starts at both accepting states. Each round expands a whole
// smaller frontier, and applies guards at the original transition's source.
struct BidirectionalCIfly {
  const Graph *g;
  QueryContext context;
  std::vector<uint32_t> forward_seen,backward_seen;
  std::vector<int> forward,backward,next;
  explicit BidirectionalCIfly(const Graph &gg)
      : g(&gg),context(gg.n),forward_seen(2*gg.n),backward_seen(2*gg.n) {}
  bool separated(int x,int y,const std::vector<int> &z){return run<false>(x,y,z);}
  // Separation in G with the out-edges of x removed (back-door graph G_{\underline x}).
  bool separated_cut_out(int x,int y,const std::vector<int> &z){return run<true>(x,y,z);}
  template<bool Cut> bool run(int x,int y,const std::vector<int> &z) {
    if (context.epoch==UINT32_MAX) {
      std::fill(forward_seen.begin(),forward_seen.end(),0);
      std::fill(backward_seen.begin(),backward_seen.end(),0);
    }
    context.set(z,x,y);
    forward.assign(1,2*x); backward={2*y,2*y+1};
    forward_seen[2*x]=context.epoch;
    backward_seen[2*y]=backward_seen[2*y+1]=context.epoch;
    while (!forward.empty() && !backward.empty()) {
      bool f=forward.size()<=backward.size(), found=false;
      auto &front=f?forward:backward;
      auto &seen=f?forward_seen:backward_seen;
      auto &other=f?backward_seen:forward_seen;
      next.clear();
      auto push=[&](int s) {
        if (other[s]==context.epoch) { found=true; return; }
        if (seen[s]!=context.epoch) {seen[s]=context.epoch;next.push_back(s);}
      };
      for (int s:front) {
        int v=s/2,d=s%2;
        // Cut: the edges x -> c do not exist (x = the source).
        const bool vcut=Cut&&v==x;
        if (f) {
          if (!context.contains(v)) {
            if (!d) for (int p:g->pa[v]) {if(Cut&&p==x)continue;push(2*p);}
            if(!vcut)for (int c:g->ch[v]) push(2*c+1);
          } else if (d) for (int p:g->pa[v]) {if(Cut&&p==x)continue;push(2*p);}
        } else if (!d) {
          if(!vcut)for (int c:g->ch[v]) push(2*c+int(context.contains(c)));
        } else {
          for (int p:g->pa[v]) if (!context.contains(p)) {if(Cut&&p==x)continue;push(2*p);push(2*p+1);}
        }
        if (found) return false;
      }
      front.swap(next);
    }
    return true;
  }
  size_t bytes() const {
    return sizeof(*this)+4*(context.mark.capacity()+forward_seen.capacity()+backward_seen.capacity()+
                            forward.capacity()+backward.capacity()+next.capacity());
  }
};
// Bidirectional traversal restricted to the ancestral set An(X ∪ Y ∪ Z)
// (Lauritzen et al. 1990: X ⊥ Y | Z in G iff in G[An(X ∪ Y ∪ Z)]). The
// ancestral set is computed per query by a reverse DFS; the traversal uses the
// same local rules as BidirectionalCIfly and never leaves the ancestral set.
// This is the textbook pruning used by networkx's d_separated.
struct AncestralBidirectional {
  const Graph *g;
  QueryContext context;
  std::vector<uint32_t> anc,forward_seen,backward_seen;uint32_t anc_epoch=0;
  std::vector<int> stack,forward,backward,next;
  size_t last_ancestral=0;   // size of the last ancestral set (diagnostics)
  explicit AncestralBidirectional(const Graph &gg)
      : g(&gg),context(gg.n),anc(gg.n),forward_seen(2*gg.n),backward_seen(2*gg.n) {}
  bool separated(int x,int y,const std::vector<int> &z){return run<false>(x,y,z);}
  // Separation in G with the out-edges of x removed; the ancestral set is taken in that graph.
  bool separated_cut_out(int x,int y,const std::vector<int> &z){return run<true>(x,y,z);}
  template<bool Cut> bool run(int x,int y,const std::vector<int> &z) {
    if (context.epoch==UINT32_MAX) {
      std::fill(forward_seen.begin(),forward_seen.end(),0);
      std::fill(backward_seen.begin(),backward_seen.end(),0);
    }
    context.set(z,x,y);
    if (++anc_epoch==0) {std::fill(anc.begin(),anc.end(),0);anc_epoch=1;}
    {PhaseTimer pt(phase_times().anc_ns);
    stack.clear();
    auto seed=[&](int v){if(anc[v]!=anc_epoch){anc[v]=anc_epoch;stack.push_back(v);}};
    seed(x);seed(y);for(int v:z)seed(v);
    size_t count=stack.size();
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      for(int p:g->pa[v])if(!(Cut&&p==x)&&anc[p]!=anc_epoch){anc[p]=anc_epoch;stack.push_back(p);++count;}}
    last_ancestral=count;}
    auto in=[&](int v){return anc[v]==anc_epoch;};
    forward.assign(1,2*x); backward={2*y,2*y+1};
    forward_seen[2*x]=context.epoch;
    backward_seen[2*y]=backward_seen[2*y+1]=context.epoch;
    while (!forward.empty() && !backward.empty()) {
      bool f=forward.size()<=backward.size(), found=false;
      auto &front=f?forward:backward;
      auto &seen=f?forward_seen:backward_seen;
      auto &other=f?backward_seen:forward_seen;
      next.clear();
      auto push=[&](int s) {
        if (!in(s/2)) return;
        if (other[s]==context.epoch) { found=true; return; }
        if (seen[s]!=context.epoch) {seen[s]=context.epoch;next.push_back(s);}
      };
      for (int s:front) {
        int v=s/2,d=s%2;
        // Cut: the edges x -> c do not exist (x = the source).
        const bool vcut=Cut&&v==x;
        if (f) {
          if (!context.contains(v)) {
            if (!d) for (int p:g->pa[v]) {if(Cut&&p==x)continue;push(2*p);}
            if(!vcut)for (int c:g->ch[v]) push(2*c+1);
          } else if (d) for (int p:g->pa[v]) {if(Cut&&p==x)continue;push(2*p);}
        } else if (!d) {
          if(!vcut)for (int c:g->ch[v]) push(2*c+int(context.contains(c)));
        } else {
          for (int p:g->pa[v]) if (!context.contains(p)) {if(Cut&&p==x)continue;push(2*p);push(2*p+1);}
        }
        if (found) return false;
      }
      front.swap(next);
    }
    return true;
  }
};
// Darwiche's pruning test (Modeling and Reasoning with Bayesian Networks, 2009), the method of
// networkx's d_separated: delete the leaves outside X ∪ Y ∪ Z until none is left, which keeps exactly
// An(X ∪ Y ∪ Z), delete the edges out of Z, and test whether X and Y are disconnected in the
// skeleton of what remains. One reverse DFS for the ancestral set, then one undirected DFS from x.
struct DarwichePruning {
  const Graph *g;
  std::vector<uint32_t> anc,zmark,seen;uint32_t epoch=0;
  std::vector<int> stack;
  size_t last_ancestral=0;   // size of the last ancestral set (diagnostics)
  explicit DarwichePruning(const Graph &gg):g(&gg),anc(gg.n),zmark(gg.n),seen(gg.n) {}
  bool separated(int x,int y,const std::vector<int> &z) {
    if (++epoch==0) {
      std::fill(anc.begin(),anc.end(),0);std::fill(zmark.begin(),zmark.end(),0);
      std::fill(seen.begin(),seen.end(),0);epoch=1;
    }
    for (int v:z) zmark[v]=epoch;
    {PhaseTimer pt(phase_times().anc_ns);
    stack.clear();
    auto seed=[&](int v){if(anc[v]!=epoch){anc[v]=epoch;stack.push_back(v);}};
    seed(x);seed(y);for(int v:z)seed(v);
    size_t count=stack.size();
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      for(int p:g->pa[v])if(anc[p]!=epoch){anc[p]=epoch;stack.push_back(p);++count;}}
    last_ancestral=count;}
    stack.assign(1,x);seen[x]=epoch;
    while (!stack.empty()) {
      int v=stack.back();stack.pop_back();
      // edge p -> v survives unless p is in Z; every parent of an ancestral vertex is ancestral
      for (int p:g->pa[v]) if (zmark[p]!=epoch&&seen[p]!=epoch) {
        if (p==y) return false;
        seen[p]=epoch;stack.push_back(p);
      }
      // edges v -> c survive unless v is in Z; c must lie in the ancestral set
      if (zmark[v]!=epoch) for (int c:g->ch[v]) if (anc[c]==epoch&&seen[c]!=epoch) {
        if (c==y) return false;
        seen[c]=epoch;stack.push_back(c);
      }
    }
    return true;
  }
};

// Geiger, Verma and Pearl (Networks 1990), the classical d-separation procedure: mark An(Z), then search
// from x over pairs (vertex, whether the last edge points into it). Two consecutive edges are legal at v
// if v is a collider in An(Z) (Z included) or a non-collider outside Z. One-sided, as published; stops
// when y is reached. separated_cut_out removes the edges out of x first (back-door mode).
struct GeigerReachable {
  const Graph *g;
  std::vector<uint32_t> anc,zmark,seen;uint32_t epoch=0;   // seen: state 2v (last edge leaves v), 2v+1 (enters v)
  std::vector<int> stack;
  explicit GeigerReachable(const Graph &gg):g(&gg),anc(gg.n),zmark(gg.n),seen(2*size_t(gg.n)) {}
  bool separated(int x,int y,const std::vector<int> &z){return run<false>(x,y,z);}
  bool separated_cut_out(int x,int y,const std::vector<int> &z){return run<true>(x,y,z);}
  template<bool CUT> bool run(int x,int y,const std::vector<int> &z) {
    if (++epoch==0) {
      std::fill(anc.begin(),anc.end(),0);std::fill(zmark.begin(),zmark.end(),0);
      std::fill(seen.begin(),seen.end(),0);epoch=1;
    }
    for (int v:z) zmark[v]=epoch;
    {PhaseTimer pt(phase_times().anc_ns);
    stack.clear();
    for (int v:z) if (anc[v]!=epoch) {anc[v]=epoch;stack.push_back(v);}
    while (!stack.empty()) {int v=stack.back();stack.pop_back();
      for (int p:g->pa[v]) if ((!CUT||p!=x)&&anc[p]!=epoch) {anc[p]=epoch;stack.push_back(p);}}}
    stack.clear();
    bool met=false;
    auto push=[&](int s){if(s/2==y)met=true;else if(seen[s]!=epoch){seen[s]=epoch;stack.push_back(s);}};
    for (int p:g->pa[x]) {push(2*p);if(met)return false;}       // every first edge is legal
    if (!CUT) for (int c:g->ch[x]) {push(2*c+1);if(met)return false;}
    while (!stack.empty()) {
      int s=stack.back();stack.pop_back();int v=s/2;bool into=s&1;
      // to a parent: v is a collider if the last edge enters v, else a non-collider
      if (into?anc[v]==epoch:zmark[v]!=epoch) for (int p:g->pa[v]) if (!CUT||p!=x) {push(2*p);if(met)return false;}
      // to a child: v is a non-collider either way
      if (zmark[v]!=epoch&&(!CUT||v!=x)) for (int c:g->ch[v]) {push(2*c+1);if(met)return false;}
    }
    return true;
  }
};

// The d-collision graph method for d-separation on causal property graphs (Neo4j UDF variant with BFS),
// ported to C++. Per query, An(Z) by a BFS from Z; the d-collision graph is computed on the fly: the
// edges out of Z are dropped, the edges into an unconditioned collider (a vertex outside Z ∪ An(Z) with
// at least two parents outside Z) are directed, and all other edges are undirected. ANTERIORS from x
// follows undirected edges and directed edges against their direction; POSTERIORS from x and its
// anteriors follows undirected edges and directed edges along their direction. y is d-connected to x iff
// one of the two searches reaches it. With `full`, both searches run to completion, as in the original,
// which returns every vertex d-separated from x; otherwise they stop when they reach y.
// separated_cut_out removes the edges out of x first (back-door mode).
struct DCollision {
  const Graph *g;
  bool full=false;
  std::vector<uint32_t> cond,zmark,colep,ant,post;std::vector<uint8_t> col;uint32_t epoch=0;
  std::vector<int> queue;
  int cut_x=-1;   // the vertex whose out-edges are removed, or -1
  DCollision(const Graph &gg,bool complete):g(&gg),full(complete),cond(gg.n),zmark(gg.n),colep(gg.n),
    ant(gg.n),post(gg.n),col(gg.n) {}
  bool kept(int p)const{return zmark[p]!=epoch&&p!=cut_x;}   // an edge out of p remains
  bool collider(int v) {   // unconditioned collider of G with the edges out of Z removed
    if (cond[v]==epoch) return false;
    if (colep[v]!=epoch) {colep[v]=epoch;int k=0;for(int p:g->pa[v])if(kept(p)&&++k>=2)break;col[v]=k>=2;}
    return col[v];
  }
  bool separated(int x,int y,const std::vector<int> &z){cut_x=-1;return run(x,y,z);}
  bool separated_cut_out(int x,int y,const std::vector<int> &z){cut_x=x;return run(x,y,z);}
  bool run(int x,int y,const std::vector<int> &z) {
    if (++epoch==0) {
      for (auto *a:{&cond,&zmark,&colep,&ant,&post}) std::fill(a->begin(),a->end(),0);
      epoch=1;
    }
    for (int v:z) zmark[v]=epoch;
    {PhaseTimer pt(phase_times().anc_ns);
    queue.clear();
    for (int v:z) if (cond[v]!=epoch) {cond[v]=epoch;queue.push_back(v);}
    for (size_t h=0;h<queue.size();++h) for (int p:g->pa[queue[h]]) if (p!=cut_x&&cond[p]!=epoch) {cond[p]=epoch;queue.push_back(p);}}
    bool reached=false;
    // ANTERIORS: to a parent p by p -> v (undirected, or directed towards v); to a child c by v -> c only
    // if the edge is undirected.
    queue.assign(1,x);ant[x]=epoch;
    for (size_t h=0;h<queue.size()&&(full||!reached);++h) {
      int v=queue[h];
      auto visit=[&](int w){if(ant[w]!=epoch){ant[w]=epoch;queue.push_back(w);if(w==y)reached=true;}};
      for (int p:g->pa[v]) if (kept(p)) visit(p);
      if (kept(v)) for (int c:g->ch[v]) if (!collider(c)) visit(c);
    }
    if (reached&&!full) return false;
    // POSTERIORS from x and the anteriors outside Z: to a child c by v -> c (undirected, or directed away
    // from v); to a parent p by p -> v only if the edge is undirected.
    std::vector<int> &q=queue;size_t w_out=0;
    for (size_t i=0;i<q.size();++i) if (zmark[q[i]]!=epoch) {post[q[i]]=epoch;q[w_out++]=q[i];}
    q.resize(w_out);
    bool reached2=false;
    for (size_t h=0;h<q.size()&&(full||!reached2);++h) {
      int v=q[h];
      auto visit=[&](int w){if(post[w]!=epoch){post[w]=epoch;q.push_back(w);if(w==y)reached2=true;}};
      if (kept(v)) for (int c:g->ch[v]) visit(c);
      if (!collider(v)) for (int p:g->pa[v]) if (kept(p)) visit(p);
    }
    return !(reached||reached2);
  }
};
} // namespace rs
