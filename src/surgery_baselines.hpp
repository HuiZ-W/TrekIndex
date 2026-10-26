#pragma once
#include "polytree.hpp"
#include "universal.hpp"
#include "reach_labels.hpp"
#include <unordered_set>
namespace rs {
// The four d-separation baselines in their general form, for the causal-criteria
// experiments: endpoint sets X and Y, conditioning set Z, and graph surgery (the
// in-edges of cut_in and the out-edges of cut_out are removed; the sets may
// overlap). Each struct follows its single-pair original line by line
// (BidirectionalCIfly, AncestralBidirectional, StreamingCIfly, RPReachable); the
// only changes are the multi-vertex start/target sets and the edge filter
// keep(p, c) = c's in-edges and p's out-edges are not cut. The single-pair
// originals stay untouched for the d-separation experiments.
struct SurgeryMarks {
  std::vector<uint32_t> z,in,out,x,y;uint32_t epoch=0;
  explicit SurgeryMarks(int n):z(n),in(n),out(n),x(n),y(n) {}
  void set(int n,const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    if(xs.empty()||ys.empty())throw std::invalid_argument("empty endpoint set");
    if(++epoch==0){for(auto *v:{&z,&in,&out,&x,&y})std::fill(v->begin(),v->end(),0);epoch=1;}
    auto chk=[&](int v){if(v<0||v>=n)throw std::invalid_argument("vertex out of range");};
    for(int v:zs){chk(v);z[v]=epoch;}
    for(int v:ci){chk(v);in[v]=epoch;}
    for(int v:co){chk(v);out[v]=epoch;}
    for(int v:xs){chk(v);if(z[v]==epoch||x[v]==epoch)throw std::invalid_argument("endpoint sets must be disjoint and outside Z");x[v]=epoch;}
    for(int v:ys){chk(v);if(z[v]==epoch||x[v]==epoch||y[v]==epoch)throw std::invalid_argument("endpoint sets must be disjoint and outside Z");y[v]=epoch;}
  }
  bool inz(int v)const{return z[v]==epoch;}
  bool iny(int v)const{return y[v]==epoch;}
  bool keep(int p,int c)const{return in[c]!=epoch&&out[p]!=epoch;}   // edge p -> c survives
};

// Bidirectional traversal with the CIfly/Bayes-ball local rules (BidirectionalCIfly).
struct BidirectionalGeneral {
  const Graph *g;SurgeryMarks k;
  std::vector<uint32_t> fseen,bseen;uint32_t ep=0;std::vector<int> forward,backward,next;
  explicit BidirectionalGeneral(const Graph &gg):g(&gg),k(gg.n),fseen(2*size_t(gg.n)),bseen(2*size_t(gg.n)) {}
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){std::fill(fseen.begin(),fseen.end(),0);std::fill(bseen.begin(),bseen.end(),0);ep=1;}
    forward.clear();backward.clear();
    for(int x:xs){forward.push_back(2*x);fseen[2*x]=ep;}
    for(int y:ys){backward.push_back(2*y);backward.push_back(2*y+1);bseen[2*y]=bseen[2*y+1]=ep;}
    while(!forward.empty()&&!backward.empty()) {
      bool f=forward.size()<=backward.size(),found=false;
      auto &front=f?forward:backward;auto &seen=f?fseen:bseen;auto &other=f?bseen:fseen;
      next.clear();
      auto push=[&](int s){if(other[s]==ep){found=true;return;}if(seen[s]!=ep){seen[s]=ep;next.push_back(s);}};
      for(int s:front) {
        int v=s/2,d=s%2;
        if(f) {
          if(!k.inz(v)) {
            if(!d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
            for(int c:g->ch[v])if(k.keep(v,c))push(2*c+1);
          } else if(d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
        } else if(!d) {
          for(int c:g->ch[v])if(k.keep(v,c))push(2*c+int(k.inz(c)));
        } else {
          for(int p:g->pa[v])if(!k.inz(p)&&k.keep(p,v)){push(2*p);push(2*p+1);}
        }
        if(found)return false;
      }
      front.swap(next);
    }
    return true;
  }
};

// The same traversal restricted to the ancestral set An(X ∪ Y ∪ Z) of the surgered
// graph (AncestralBidirectional; Lauritzen et al. 1990).
struct AncestralGeneral {
  const Graph *g;SurgeryMarks k;
  std::vector<uint32_t> anc,fseen,bseen;uint32_t ep=0,aep=0;std::vector<int> stack,forward,backward,next;
  explicit AncestralGeneral(const Graph &gg):g(&gg),k(gg.n),anc(gg.n),fseen(2*size_t(gg.n)),bseen(2*size_t(gg.n)) {}
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){std::fill(fseen.begin(),fseen.end(),0);std::fill(bseen.begin(),bseen.end(),0);ep=1;}
    if(++aep==0){std::fill(anc.begin(),anc.end(),0);aep=1;}
    {PhaseTimer pt(phase_times().anc_ns);
    stack.clear();
    auto seed=[&](int v){if(anc[v]!=aep){anc[v]=aep;stack.push_back(v);}};
    for(int v:xs)seed(v);for(int v:ys)seed(v);for(int v:zs)seed(v);
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      for(int p:g->pa[v])if(k.keep(p,v)&&anc[p]!=aep){anc[p]=aep;stack.push_back(p);}}}
    auto in=[&](int v){return anc[v]==aep;};
    forward.clear();backward.clear();
    for(int x:xs){forward.push_back(2*x);fseen[2*x]=ep;}
    for(int y:ys){backward.push_back(2*y);backward.push_back(2*y+1);bseen[2*y]=bseen[2*y+1]=ep;}
    while(!forward.empty()&&!backward.empty()) {
      bool f=forward.size()<=backward.size(),found=false;
      auto &front=f?forward:backward;auto &seen=f?fseen:bseen;auto &other=f?bseen:fseen;
      next.clear();
      auto push=[&](int s){if(!in(s/2))return;if(other[s]==ep){found=true;return;}if(seen[s]!=ep){seen[s]=ep;next.push_back(s);}};
      for(int s:front) {
        int v=s/2,d=s%2;
        if(f) {
          if(!k.inz(v)) {
            if(!d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
            for(int c:g->ch[v])if(k.keep(v,c))push(2*c+1);
          } else if(d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
        } else if(!d) {
          for(int c:g->ch[v])if(k.keep(v,c))push(2*c+int(k.inz(c)));
        } else {
          for(int p:g->pa[v])if(!k.inz(p)&&k.keep(p,v)){push(2*p);push(2*p+1);}
        }
        if(found)return false;
      }
      front.swap(next);
    }
    return true;
  }
};

// CIfly's d-separation rule table, depth-first from X with workspace reuse (StreamingCIfly).
struct CIflyGeneral {
  const Graph *g;SurgeryMarks k;std::vector<uint32_t> seen;uint32_t ep=0;std::vector<int> todo;
  explicit CIflyGeneral(const Graph &gg):g(&gg),k(gg.n),seen(2*size_t(gg.n)) {}
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){std::fill(seen.begin(),seen.end(),0);ep=1;}
    todo.clear();
    auto push=[&](int s){if(seen[s]!=ep){seen[s]=ep;todo.push_back(s);}};
    for(int x:xs)push(2*x);
    while(!todo.empty()) {
      int s=todo.back();todo.pop_back();int v=s/2,d=s%2;
      if(k.iny(v))return false;
      if(!k.inz(v)) {
        if(!d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
        for(int c:g->ch[v])if(k.keep(v,c))push(2*c+1);
      } else if(d)for(int p:g->pa[v])if(k.keep(p,v))push(2*p);
    }
    return true;
  }
};

// Relevant-path separation (Butz et al., PGM 2016, Algorithm 1; RPReachable, Boolean early exit)
// on the surgered graph.
struct RPGeneral {
  const Graph *g;SurgeryMarks k;
  std::vector<uint32_t> ancestor_z,relevant,visited;uint32_t ep=0;std::vector<int> todo,stack;
  explicit RPGeneral(const Graph &gg):g(&gg),k(gg.n),ancestor_z(gg.n),relevant(gg.n),visited(2*size_t(gg.n)) {}
  void ancestors(std::vector<uint32_t> &mark,const std::vector<int> &seeds) {
    stack.clear();
    for(int v:seeds)if(mark[v]!=ep){mark[v]=ep;stack.push_back(v);}
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      for(int p:g->pa[v])if(k.keep(p,v)&&mark[p]!=ep){mark[p]=ep;stack.push_back(p);}}
  }
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){for(auto *p:{&ancestor_z,&relevant,&visited})std::fill(p->begin(),p->end(),0);ep=1;}
    {PhaseTimer pt(phase_times().anc_ns);
    ancestors(ancestor_z,zs);                                     // paper line 3
    ancestors(relevant,zs);ancestors(relevant,xs);ancestors(relevant,ys);}   // line 4
    todo.clear();for(int x:xs){visited[2*x]=ep;todo.push_back(2*x);}
    while(!todo.empty()) {
      int s=todo.back();todo.pop_back();int v=s/2,d=s%2;
      if(!k.inz(v)&&k.iny(v))return false;
      auto push=[&](int w,int dir){int t=2*w+dir;if(relevant[w]==ep&&visited[t]!=ep){visited[t]=ep;todo.push_back(t);}};
      if(!d&&!k.inz(v)) {
        for(int p:g->pa[v])if(k.keep(p,v))push(p,0);
        for(int c:g->ch[v])if(k.keep(v,c))push(c,1);
      } else if(d) {
        if(!k.inz(v))for(int c:g->ch[v])if(k.keep(v,c))push(c,1);
        if(ancestor_z[v]==ep)for(int p:g->pa[v])if(k.keep(p,v))push(p,0);
      }
    }
    return true;
  }
};

// Geiger et al.'s procedure (GeigerReachable) for endpoint sets and surgery: An(Z) in the surgered graph,
// then a one-sided search over (vertex, whether the last edge enters it) from X.
struct GeigerGeneral {
  const Graph *g;SurgeryMarks k;
  std::vector<uint32_t> anc,seen;uint32_t ep=0;std::vector<int> stack;
  explicit GeigerGeneral(const Graph &gg):g(&gg),k(gg.n),anc(gg.n),seen(2*size_t(gg.n)) {}
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){std::fill(anc.begin(),anc.end(),0);std::fill(seen.begin(),seen.end(),0);ep=1;}
    {PhaseTimer pt(phase_times().anc_ns);
    stack.clear();for(int z:zs)if(anc[z]!=ep){anc[z]=ep;stack.push_back(z);}
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      for(int p:g->pa[v])if(k.keep(p,v)&&anc[p]!=ep){anc[p]=ep;stack.push_back(p);}}}
    stack.clear();bool met=false;
    auto push=[&](int s){if(k.iny(s/2))met=true;else if(seen[s]!=ep){seen[s]=ep;stack.push_back(s);}};
    for(int x:xs) {
      for(int p:g->pa[x])if(k.keep(p,x)){push(2*p);if(met)return false;}
      for(int c:g->ch[x])if(k.keep(x,c)){push(2*c+1);if(met)return false;}
    }
    while(!stack.empty()) {
      int s=stack.back();stack.pop_back();int v=s/2;bool into=s&1;
      if(into?anc[v]==ep:!k.inz(v))for(int p:g->pa[v])if(k.keep(p,v)){push(2*p);if(met)return false;}
      if(!k.inz(v))for(int c:g->ch[v])if(k.keep(v,c)){push(2*c+1);if(met)return false;}
    }
    return true;
  }
};

// The d-collision graph method (DCollision, UDF variant with BFS, stopping at Y) for endpoint sets and
// surgery: the d-collision graph of the surgered graph, computed on the fly.
struct DCollisionGeneral {
  const Graph *g;SurgeryMarks k;
  std::vector<uint32_t> cond,colep,ant,post;std::vector<uint8_t> col;uint32_t ep=0;std::vector<int> q;
  explicit DCollisionGeneral(const Graph &gg):g(&gg),k(gg.n),cond(gg.n),colep(gg.n),ant(gg.n),post(gg.n),col(gg.n) {}
  bool kept(int p,int c)const{return !k.inz(p)&&k.keep(p,c);}   // edge p -> c of the surgered graph, not out of Z
  bool collider(int v) {
    if(cond[v]==ep)return false;
    if(colep[v]!=ep){colep[v]=ep;int n=0;for(int p:g->pa[v])if(kept(p,v)&&++n>=2)break;col[v]=n>=2;}
    return col[v];
  }
  bool separated(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co) {
    k.set(g->n,xs,ys,zs,ci,co);
    if(++ep==0){for(auto *a:{&cond,&colep,&ant,&post})std::fill(a->begin(),a->end(),0);ep=1;}
    {PhaseTimer pt(phase_times().anc_ns);
    q.clear();for(int z:zs)if(cond[z]!=ep){cond[z]=ep;q.push_back(z);}
    for(size_t h=0;h<q.size();++h)for(int p:g->pa[q[h]])if(k.keep(p,q[h])&&cond[p]!=ep){cond[p]=ep;q.push_back(p);}}
    bool met=false;
    q.clear();for(int x:xs)if(ant[x]!=ep){ant[x]=ep;q.push_back(x);}
    for(size_t h=0;h<q.size()&&!met;++h) {                 // ANTERIORS
      int v=q[h];
      auto visit=[&](int w){if(ant[w]!=ep){ant[w]=ep;q.push_back(w);if(k.iny(w))met=true;}};
      for(int p:g->pa[v])if(kept(p,v))visit(p);
      for(int c:g->ch[v])if(kept(v,c)&&!collider(c))visit(c);
    }
    if(met)return false;
    size_t w_out=0;for(size_t i=0;i<q.size();++i)if(!k.inz(q[i])){post[q[i]]=ep;q[w_out++]=q[i];}
    q.resize(w_out);
    for(size_t h=0;h<q.size()&&!met;++h) {                 // POSTERIORS
      int v=q[h];
      auto visit=[&](int w){if(post[w]!=ep){post[w]=ep;q.push_back(w);if(k.iny(w))met=true;}};
      for(int c:g->ch[v])if(kept(v,c))visit(c);
      if(!collider(v))for(int p:g->pa[v])if(kept(p,v))visit(p);
    }
    return !met;
  }
};

// The nine causal queries answered with a baseline as the separation engine; ancestry
// and interception tests by graph search (identical for every baseline). Mirrors
// CausalQueries (src/causal.hpp) and TraversalCausal (src/zskeleton.cpp).
template<class Sep> struct BaselineCausal {
  const Graph *g;Sep sep_;std::vector<uint32_t> seen;uint32_t ep=0;std::vector<int> stack;
  explicit BaselineCausal(const Graph &gg):g(&gg),sep_(gg),seen(gg.n) {}
  // +R: the directed-path (Reach) tests of a criterion from the reachability labels of G where plain
  // reachability decides them (descendants, no forbidden set), and the labels as a refutation before
  // the search otherwise (paths that avoid a set, Z' of rule 3). The d-separation searches, including
  // their ancestor marking, are those of the baseline, unchanged.
  const ReachLabels *R=nullptr;
  void use_labels(const ReachLabels &r){R=&r;}
  uint32_t fresh(){if(++ep==0){std::fill(seen.begin(),seen.end(),0);ep=1;}return ep;}
  bool sep(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &ci,const std::vector<int> &co){PhaseTimer pt(phase_times().sep_ns);return sep_.separated(xs,ys,zs,ci,co);}
  bool descendant(int u,int v) {
    if(R)return u!=v&&R->reach(u,v);
    auto e=fresh();stack.assign(1,u);
    while(!stack.empty()){int a=stack.back();stack.pop_back();for(int c:g->ch[a]){if(c==v)return true;if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}}
    return false;
  }
  bool reach_avoiding(int u,const std::vector<int> &targets,const std::vector<int> &avoid) {
    if(R){bool any=false;for(int t:targets)if(t!=u&&R->reach(u,t)){any=true;break;}
      if(!any)return false;if(avoid.empty())return true;}
    auto e=fresh();for(int a:avoid)seen[a]=e;std::unordered_set<int> t(targets.begin(),targets.end());stack.assign(1,u);
    while(!stack.empty()){int a=stack.back();stack.pop_back();for(int c:g->ch[a]){if(t.count(c))return true;if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}}
    return false;
  }
  bool backdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs) {
    if(R){for(int x:xs)for(int z:zs)if(z!=x&&R->reach(x,z))return false;return sep(xs,ys,zs,{},xs);}
    auto e=fresh();std::unordered_set<int> zset(zs.begin(),zs.end());
    stack.clear();for(int x:xs)for(int c:g->ch[x])if(seen[c]!=e){seen[c]=e;stack.push_back(c);}
    while(!stack.empty()){int v=stack.back();stack.pop_back();if(zset.count(v))return false;
      for(int c:g->ch[v])if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}
    return sep(xs,ys,zs,{},xs);
  }
  bool rule(int which,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &xs,const std::vector<int> &ws) {
    auto cond=xs;cond.insert(cond.end(),ws.begin(),ws.end());
    if(which==1)return sep(ys,zs,cond,xs,{});
    if(which==2)return sep(ys,zs,cond,xs,zs);
    auto e=fresh();std::unordered_set<int> xset(xs.begin(),xs.end());
    // Z' = the z not in An(W) of G without the in-edges of X. +R: a z that reaches no w in G is in Z'
    // without a search; the backward search runs only if some z is left undecided.
    std::vector<int> cut=xs,open;
    if(R){for(int z:zs){if(xset.count(z)){cut.push_back(z);continue;}
        bool r=false;for(int w:ws)if(z==w||R->reach(z,w)){r=true;break;}
        if(r)open.push_back(z);else cut.push_back(z);}
      if(open.empty())return sep(ys,zs,cond,cut,{});}
    stack.clear();for(int w:ws){seen[w]=e;stack.push_back(w);}
    while(!stack.empty()){int v=stack.back();stack.pop_back();if(xset.count(v))continue;
      for(int p:g->pa[v])if(seen[p]!=e){seen[p]=e;stack.push_back(p);}}
    if(R){for(int z:open)if(seen[z]!=e)cut.push_back(z);}
    else for(int z:zs)if(seen[z]!=e||xset.count(z))cut.push_back(z);
    return sep(ys,zs,cond,cut,{});
  }
  bool frontdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ms) {
    for(int x:xs)if(reach_avoiding(x,ys,ms))return false;
    return ms.empty()||(sep(xs,ms,{},{},xs)&&sep(ms,ys,xs,{},ms));
  }
  bool sequential_backdoor(const std::vector<int> &xs,const std::vector<std::vector<int>> &ws,const std::vector<int> &ys) {
    std::vector<int> cond;
    for(size_t i=0;i<xs.size();++i) {
      for(size_t j=i;j<xs.size();++j)for(int w:ws[i])if(descendant(xs[j],w))return false;
      cond.insert(cond.end(),ws[i].begin(),ws[i].end());
      if(!sep(ys,{xs[i]},cond,std::vector<int>(xs.begin()+long(i)+1,xs.end()),{xs[i]}))return false;
      cond.push_back(xs[i]);
    }
    return true;
  }
  bool mediation(const std::vector<int> &xs,const std::vector<int> &ms,const std::vector<int> &ys,const std::vector<int> &ws) {
    for(int w:ws){for(int x:xs)if(descendant(x,w))return false;for(int m:ms)if(descendant(m,w))return false;}
    auto xw=xs;xw.insert(xw.end(),ws.begin(),ws.end());
    return sep(ms,ys,xw,{},ms)&&sep(xs,ms,ws,{},xs)&&sep(xs,ys,ws,ms,xs);
  }
  bool selection_recoverable(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ps) {
    std::unordered_set<int> in_x(xs.begin(),xs.end()),in_y(ys.begin(),ys.end()),done;std::vector<int> ends;
    for(int p:ps){if(in_y.count(p))return false;if(!in_x.count(p)&&done.insert(p).second)ends.push_back(p);}
    return ends.empty()||sep(ys,ends,xs,{},{});
  }
};
} // namespace rs
