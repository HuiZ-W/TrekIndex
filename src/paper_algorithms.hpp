#pragma once
#include "universal.hpp"
namespace rs {
// Butz et al., PGM 2016, Algorithm 1. Their I(X,Y,Z) uses Y as
// conditioning set; our API is separated(x,y,z).
struct RPReachable {
  const Graph &g;
  QueryContext ctx;
  std::vector<uint32_t> ancestor_z, relevant, visited, reached;
  std::vector<int> todo, stack;
  bool early;
  RPReachable(const Graph &gg,bool stop) : g(gg),ctx(g.n),ancestor_z(g.n),
    relevant(g.n),visited(2*g.n),reached(g.n),early(stop) {}
  // cut: a vertex whose out-edges are removed (-1: none).
  void ancestors(std::vector<uint32_t> &mark,const std::vector<int> &seeds,int cut=-1) {
    stack.clear();
    for(int v:seeds) if(mark[v]!=ctx.epoch){mark[v]=ctx.epoch;stack.push_back(v);}
    while(!stack.empty()) {
      int v=stack.back();stack.pop_back();
      for(int p:g.pa[v]) if(p!=cut&&mark[p]!=ctx.epoch){mark[p]=ctx.epoch;stack.push_back(p);}
    }
  }
  bool separated(int x,int y,const std::vector<int> &z){return run<false>(x,y,z);}
  // RP on G with the out-edges of x removed (back-door graph G_{\underline x}).
  bool separated_cut_out(int x,int y,const std::vector<int> &z){return run<true>(x,y,z);}
  template<bool Cut> bool run(int x,int y,const std::vector<int> &z) {
    const int cut=Cut?x:-1;
    if(ctx.epoch==UINT32_MAX) for(auto *p:{&ancestor_z,&relevant,&visited,&reached})
      std::fill(p->begin(),p->end(),0);
    ctx.set(z,x,y);
    ancestors(ancestor_z,z,cut); // paper line 3
    ancestors(relevant,z,cut);ancestors(relevant,{x,y},cut); // line 4
    todo.assign(1,2*x);visited[2*x]=ctx.epoch;
    while(!todo.empty()) {
      int s=todo.back();todo.pop_back();int v=s/2,d=s%2;
      if(!ctx.contains(v)) {reached[v]=ctx.epoch;if(early && v==y)return false;}
      auto push=[&](int w,int dir){int t=2*w+dir;
        if(relevant[w]==ctx.epoch && visited[t]!=ctx.epoch){visited[t]=ctx.epoch;todo.push_back(t);}};
      if(!d && !ctx.contains(v)) {
        for(int p:g.pa[v])if(!(Cut&&p==x))push(p,0);
        if(!(Cut&&v==x))for(int c:g.ch[v])push(c,1);
      } else if(d) {
        if(!ctx.contains(v)&&!(Cut&&v==x))for(int c:g.ch[v])push(c,1);
        if(ancestor_z[v]==ctx.epoch)for(int p:g.pa[v])if(!(Cut&&p==x))push(p,0);
      }
    }
    return reached[y]!=ctx.epoch;
  }
  size_t bytes() const {return sizeof(*this)+4*(ctx.mark.capacity()+ancestor_z.capacity()+
    relevant.capacity()+visited.capacity()+reached.capacity()+todo.capacity()+stack.capacity());}
};

// C++ port of the official CIfly 0.1.3 reach.rs BFS, specialized ONLY
// to its published DAG d-separation rule table (one color, two arrival types).
// Per-query visited/output arrays, FIFO and full output materialization retained.
// First matching rule wins: down -> up uses in-Z; other pairs use not-in-Z.
struct CIflyPaperCPP {
  const Graph &g;
  std::vector<size_t> offset;
  std::vector<std::pair<int,int>> neighbors;
  bool early;
  CIflyPaperCPP(const Graph &gg,bool stop=false):g(gg),offset(g.n+1),early(stop) {
    for(int v=0;v<g.n;++v)offset[v+1]=offset[v]+g.pa[v].size()+g.ch[v].size();
    neighbors.resize(offset.back());auto cursor=offset;
    // Same CSR scatter as instance.rs Graph::new, using sorted edge list.
    for(int u=0;u<g.n;++u)for(int v:g.ch[u]){
      neighbors[cursor[u]++]={v,1};neighbors[cursor[v]++]={u,0};
    }
  }
  bool separated(int x,int y,const std::vector<int>&zs)const{
    auto z=g.context(zs);g.validate(x,y,z);
    std::vector<uint8_t> seen(2*g.n),added(g.n);
    std::deque<std::pair<int,int>> queue;std::vector<int> result;
    queue.emplace_back(x,0);seen[2*x]=1;
    while(!queue.empty()){
      auto [v,d]=queue.front();queue.pop_front();
      if(!added[v]){result.push_back(v);added[v]=1;if(early&&v==y)return false;}
      for(size_t i=offset[v];i<offset[v+1];++i){auto [w,e]=neighbors[i];
        if(!seen[2*w+e] && ((d==1&&e==0)?bool(z[v]):!z[v])){
          seen[2*w+e]=1;queue.emplace_back(w,e);
        }
      }
    }
    return std::find(result.begin(),result.end(),y)==result.end();
  }
  size_t bytes()const{return sizeof(*this)+sizeof(size_t)*offset.capacity()+sizeof(neighbors[0])*neighbors.capacity();}
};

// Classical control: unique tree path + ancestor intervals/binary-lifting LCA.
// Independently implemented; no projection/collider-count trick from PolytreeIndex.
// This composition is a control baseline, NOT attributed to a published d-sep index.
struct TreePathBaseline {
  const Graph &g;
  QueryContext ctx;
  std::vector<int> parent,depth,component,tin,tout,prefix,path,right;
  std::vector<uint8_t> down;
  std::vector<std::vector<int>> jump;
  bool collider_free=true, classic=false;
  TreePathBaseline(const Graph &gg,bool optimized=false):g(gg),ctx(g.n),parent(g.n,-1),depth(g.n),component(g.n,-1),
    tin(g.n),tout(g.n),prefix(g.n),down(g.n) {
    classic=optimized;
    for(auto &p:g.pa)if(p.size()>1)collider_free=false;
    int timer=0,maxdepth=0;
    for(int ri=0;ri<g.n;++ri){int root=classic?g.order[ri]:ri;if(parent[root]<0) {
      parent[root]=root;component[root]=root;tin[root]=timer++;
      std::vector<std::pair<int,size_t>> st{{root,0}};
      while(!st.empty()) {
        int v=st.back().first;size_t &i=st.back().second;
        if(i==g.pa[v].size()+g.ch[v].size()){tout[v]=timer;st.pop_back();continue;}
        bool f=i>=g.pa[v].size();int w=f?g.ch[v][i-g.pa[v].size()]:g.pa[v][i];++i;
        if(w==parent[v])continue;
        if(parent[w]>=0)throw std::invalid_argument("tree baseline requires forest skeleton");
        parent[w]=v;component[w]=root;depth[w]=depth[v]+1;down[w]=f;
        prefix[w]=prefix[v]+f;tin[w]=timer++;maxdepth=std::max(maxdepth,depth[w]);st.emplace_back(w,0);
      }
    }
    }
    int levels=1;while((uint64_t(1)<<levels)<=uint64_t(maxdepth))++levels;
    jump.assign(levels,parent);
    for(int k=1;k<levels;++k)for(int v=0;v<g.n;++v)jump[k][v]=jump[k-1][jump[k-1][v]];
  }
  bool anc(int a,int b)const{return tin[a]<=tin[b] && tout[b]<=tout[a];}
  int lca(int a,int b)const {
    if(component[a]!=component[b])return -1;
    if(anc(a,b))return a;if(anc(b,a))return b;
    for(int k=int(jump.size())-1;k>=0;--k)if(!anc(jump[k][a],b))a=jump[k][a];
    return parent[a];
  }
  bool reaches(int a,int b)const {
    int l=lca(a,b);return l>=0 && prefix[a]==prefix[l] && prefix[b]-prefix[l]==depth[b]-depth[l];
  }
  bool arrow(int a,int b)const{return parent[b]==a?bool(down[b]):!down[a];}
  bool separated(int x,int y,const std::vector<int>&z) {
    if(classic && collider_free){
      if(x<0||y<0||x>=g.n||y>=g.n||x==y)throw std::invalid_argument("endpoints");
      for(int v:z)if(v<0||v>=g.n||v==x||v==y)throw std::invalid_argument("Z");
    }else ctx.set(z,x,y);
    int l=lca(x,y);if(l<0)return true;
    if(collider_free) {
      for(int v:z)if(component[v]==component[x] && depth[v]>=depth[l] && (anc(v,x)||anc(v,y)))return true;
      return false;
    }
    path.clear();right.clear();
    for(int v=x;v!=l;v=parent[v])path.push_back(v);
    path.push_back(l);
    for(int v=y;v!=l;v=parent[v])right.push_back(v);
    path.insert(path.end(),right.rbegin(),right.rend());
    for(size_t i=1;i+1<path.size();++i) {
      int v=path[i];bool collider=arrow(path[i-1],v)&&arrow(path[i+1],v);
      if(!collider){if(ctx.contains(v))return true;}
      else {bool active=false;for(int w:z)if(reaches(v,w)){active=true;break;}if(!active)return true;}
    }
    return false;
  }
  size_t bytes()const{
    size_t b=sizeof(*this)+down.capacity()+4*ctx.mark.capacity()+jump.capacity()*sizeof(std::vector<int>);
    for(auto *v:{&parent,&depth,&component,&tin,&tout,&prefix,&path,&right})b+=4*v->capacity();
    for(auto &v:jump)b+=4*v.capacity();return b;
  }
};
} // namespace rs
