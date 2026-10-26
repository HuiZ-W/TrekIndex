#include "causal.hpp"
#include "coretree.hpp"
#include "backdoor_baselines.hpp"
#include "surgery_baselines.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
using namespace rs;
// "+R" baselines (method name ending in "-r"): the reachability labels of G, built once per process
// on first use; their build time is reported as the build time of every "-r" method.
static bool plus_r(const std::string &m){return m.size()>2&&m.compare(m.size()-2,2,"-r")==0;}
static std::string base_of(const std::string &m){return plus_r(m)?m.substr(0,m.size()-2):m;}
static const ReachLabels &reach_labels(const Graph &g){
  static std::unique_ptr<ReachLabels> rl;
  if(!rl){rl=std::make_unique<ReachLabels>(g);
    std::cerr<<"{\"event\":\"reach_labels_built\",\"ms\":"<<rl->build_ms<<",\"entries\":"<<rl->entries()<<",\"bytes\":"<<rl->bytes()<<"}\n";}
  return *rl;}
// One process: load G, build the graph-only index once, generate (or read) a
// seeded arbitrary-Z workload, then time every method on identical requests.
// Every answer is cross-checked against the first traversal baseline.
struct Request { std::vector<int> z; int x, y; };
struct Options {
  std::string graph, methods="bidirectional,zskel", zmode="random", zsizes="4", pairmode="random", workload_in, workload_out;
  size_t queries=1000; int zsize=4; uint64_t seed=1; bool progress=true,emit_answers=false,bidi=true,live=true,arcs=true; uint64_t probe=256; bool probe_per_z=false,compress=false; int warmup=0; uint64_t large_h=64,initial_budget=65536; std::string hub_order="degree"; size_t closure_budget=8; int width=0; size_t probe_threshold=0; bool ancmask=true,surgery_first=true,dominated=true,parent_live=false;   // width 0 = automatic
  int bp=0;   // bit-parallel top hubs in the labels (0..64)
  bool warmup_disjoint=false;   // warm up on a different request set (seed + 1000) instead of the timed one
};
Graph read_binary(const std::string &path) {
  std::ifstream in(path,std::ios::binary);char magic[8];uint64_t n=0,m=0;
  if(!in.read(magic,8)||std::string(magic,8)!="FDGBIN01")throw std::invalid_argument("binary magic");
  if(!in.read(reinterpret_cast<char*>(&n),8)||!in.read(reinterpret_cast<char*>(&m),8))throw std::invalid_argument("header");
  // Graph ids are int; trek-state ids (4n) are checked where an index is built.
  if(n>uint64_t(std::numeric_limits<int>::max()/2))throw std::length_error("vertex count");
  std::vector<uint64_t> buf(2*m);
  if(!in.read(reinterpret_cast<char*>(buf.data()),std::streamsize(16*m)))throw std::invalid_argument("truncated edges");
  std::vector<std::pair<int,int>> es(m);
  for(uint64_t i=0;i<m;++i){if(buf[2*i]>=n||buf[2*i+1]>=n)throw std::invalid_argument("endpoint");es[i]={int(buf[2*i]),int(buf[2*i+1])};}
  buf.clear();buf.shrink_to_fit();
  return Graph(int(n),es);
}
uint64_t peak_rss() {
  rusage r{};getrusage(RUSAGE_SELF,&r);
#ifdef __APPLE__
  return r.ru_maxrss;
#else
  std::ifstream in("/proc/self/status");std::string line;
  while(std::getline(in,line))if(line.rfind("VmHWM:",0)==0){std::istringstream v(line.substr(6));uint64_t kib;if(v>>kib)return kib*1024;}
  return uint64_t(r.ru_maxrss)*1024;
#endif
}
std::vector<Request> generate(const Graph &g,const Options &o) {
  std::mt19937_64 rng(o.seed);
  std::vector<int> live;for(int v=0;v<g.n;++v)if(!g.pa[v].empty()||!g.ch[v].empty())live.push_back(v);
  if(live.size()<2)throw std::invalid_argument("graph has fewer than two non-isolated vertices");
  if(o.zsize>g.n-2)throw std::invalid_argument("zsize "+std::to_string(o.zsize)+" exceeds n-2 for this graph");
  auto pick=[&]{return live[rng()%live.size()];};
  auto step=[&](int v,bool up_only){
    size_t a=g.pa[v].size(),b=up_only?0:g.ch[v].size();
    if(!a&&!b)return v;
    auto r=rng()%(a+b);return r<a?g.pa[v][r]:g.ch[v][r-a];
  };
  auto walk=[&](int v,int len,bool up_only){for(int i=0;i<len;++i)v=step(v,up_only);return v;};
  std::vector<Request> qs(o.queries);
  std::unordered_set<int> used;
  for(auto &q:qs) {
    q.x=pick();q.y=q.x;
    if(o.pairmode=="walk")for(int t=0;t<16&&q.y==q.x;++t)q.y=walk(q.x,1+int(rng()%6),false);
    while(q.y==q.x)q.y=pick();
    used.clear();used.insert(q.x);used.insert(q.y);q.z.clear();
    int budget=64*o.zsize+64;
    while(int(q.z.size())<o.zsize && budget-->0) {
      int v;
      // Draw the walk length before the start vertex, in two statements: the
      // order of evaluation of function arguments is unspecified in C++, and
      // this is the order g++ used for every saved workload.
      if(o.zmode=="near"){int len=1+int(rng()%4);int from=rng()%2?q.x:q.y;v=walk(from,len,false);}
      else if(o.zmode=="ancestor"){int len=1+int(rng()%8);int from=rng()%2?q.x:q.y;v=walk(from,len,true);}
      else v=int(rng()%uint64_t(g.n));
      if(used.insert(v).second)q.z.push_back(v);
    }
    while(int(q.z.size())<o.zsize){int v=int(rng()%uint64_t(g.n));if(used.insert(v).second)q.z.push_back(v);}
  }
  return qs;
}
std::vector<Request> read_workload(const std::string &path) {
  std::ifstream in(path);size_t count;if(!(in>>count))throw std::invalid_argument("workload count");
  std::vector<Request> qs(count);
  for(auto &q:qs){int k;if(!(in>>k)||k<0)throw std::invalid_argument("context size");q.z.resize(k);for(auto &v:q.z)in>>v;if(!(in>>q.x>>q.y))throw std::invalid_argument("endpoints");}
  return qs;
}
void write_workload(const std::string &path,const std::vector<Request> &qs) {
  std::ofstream out(path);out<<qs.size()<<'\n';
  for(auto &q:qs){out<<q.z.size();for(int v:q.z)out<<' '<<v;out<<' '<<q.x<<' '<<q.y<<'\n';}
}
double pct(std::vector<double> v,double p){if(v.empty())return 0;std::sort(v.begin(),v.end());return v[std::min(v.size()-1,size_t(p*double(v.size())))];}

// Ablation switches from the command line; defaults are the full method.
void configure(ZSkeleton &z,const Options &o){z.balance_by_arcs=o.arcs;z.bidirectional_search=o.bidi;z.live_states=o.live;}

// Exact traversal-only implementation of the same causal criteria: the
// baseline a system without an index would run. Descendant / ancestor side
// conditions by BFS on G; separation by budget-free bidirectional traversal
// with surgery.
struct TraversalCausal {
  const Graph *g;
  ZSkeleton::Scratch s;
  BudgetedTraversal t;
  std::vector<uint32_t> seen;uint32_t ep=0;std::vector<int> stack;
  explicit TraversalCausal(const Graph &gg):g(&gg),s(gg.n),t(gg),seen(gg.n) {}
  uint32_t fresh(){if(++ep==0){std::fill(seen.begin(),seen.end(),0);ep=1;}return ep;}
  bool sep(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs,const Surgery &cut) {
    ZSkeleton::mark_marks(*g,xs,ys,zs,cut,s);return t.run(xs,ys,s,UINT64_MAX)==0;
  }
  BackdoorResult backdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs) {
    BackdoorResult r;auto e=fresh();std::unordered_set<int> zset(zs.begin(),zs.end());
    stack.clear();for(int x:xs){for(int c:g->ch[x])if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      if(zset.count(v)){r.descendant_in_z=true;r.z=v;return r;}
      for(int c:g->ch[v])if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}
    r.open_backdoor_path=!sep(xs,ys,zs,Surgery{{},xs});r.valid=!r.open_backdoor_path;return r;
  }
  bool rule(int which,const std::vector<int> &ys,const std::vector<int> &zs,const std::vector<int> &xs,const std::vector<int> &ws) {
    auto cond=xs;cond.insert(cond.end(),ws.begin(),ws.end());
    if(which==1)return sep(ys,zs,cond,Surgery{xs,{}});
    if(which==2)return sep(ys,zs,cond,Surgery{xs,zs});
    // Ancestors of W in G with in-edges of X removed: reverse BFS avoiding X.
    auto e=fresh();std::unordered_set<int> xset(xs.begin(),xs.end());
    stack.clear();for(int w:ws){seen[w]=e;stack.push_back(w);}
    while(!stack.empty()){int v=stack.back();stack.pop_back();
      if(xset.count(v))continue;
      for(int p:g->pa[v])if(seen[p]!=e){seen[p]=e;stack.push_back(p);}}
    auto cut=xs;for(int z:zs)if(seen[z]!=e||xset.count(z))cut.push_back(z);
    return sep(ys,zs,cond,Surgery{cut,{}});
  }
  // The four further criteria of CausalQueries, by graph search and bidirectional traversal.
  bool descendant(int u,int v) {   // v is a proper descendant of u
    auto e=fresh();stack.assign(1,u);
    while(!stack.empty()){int a=stack.back();stack.pop_back();for(int c:g->ch[a]){if(c==v)return true;if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}}
    return false;
  }
  bool reach_avoiding(int u,const std::vector<int> &targets,const std::vector<int> &avoid) {
    auto e=fresh();for(int a:avoid)seen[a]=e;std::unordered_set<int> t(targets.begin(),targets.end());stack.assign(1,u);
    while(!stack.empty()){int a=stack.back();stack.pop_back();for(int c:g->ch[a]){if(t.count(c))return true;if(seen[c]!=e){seen[c]=e;stack.push_back(c);}}}
    return false;
  }
  bool frontdoor(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ms) {
    for(int x:xs)if(reach_avoiding(x,ys,ms))return false;
    return ms.empty()||(sep(xs,ms,{},Surgery{{},xs})&&sep(ms,ys,xs,Surgery{{},ms}));
  }
  bool sequential_backdoor(const std::vector<int> &xs,const std::vector<std::vector<int>> &ws,const std::vector<int> &ys) {
    std::vector<int> cond;
    for(size_t k=0;k<xs.size();++k) {
      for(size_t j=k;j<xs.size();++j)for(int w:ws[k])if(descendant(xs[j],w))return false;
      cond.insert(cond.end(),ws[k].begin(),ws[k].end());
      if(!sep(ys,{xs[k]},cond,Surgery{std::vector<int>(xs.begin()+long(k)+1,xs.end()),{xs[k]}}))return false;
      cond.push_back(xs[k]);
    }
    return true;
  }
  bool mediation(const std::vector<int> &xs,const std::vector<int> &ms,const std::vector<int> &ys,const std::vector<int> &ws) {
    for(int w:ws){for(int x:xs)if(descendant(x,w))return false;for(int m:ms)if(descendant(m,w))return false;}
    auto xw=xs;xw.insert(xw.end(),ws.begin(),ws.end());
    return sep(ms,ys,xw,Surgery{{},ms})&&sep(xs,ms,ws,Surgery{{},xs})&&sep(xs,ys,ws,Surgery{ms,xs});
  }
  bool selection_recoverable(const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &ps) {
    std::unordered_set<int> in_x(xs.begin(),xs.end()),in_y(ys.begin(),ys.end()),done;std::vector<int> ends;
    for(int p:ps){if(in_y.count(p))return false;if(!in_x.count(p)&&done.insert(p).second)ends.push_back(p);}
    return ends.empty()||sep(ys,ends,xs,Surgery{});
  }
};

int causal_main(const Options &o,const std::vector<std::string> &kinds) {
  auto list=[](const std::string &v){std::vector<std::string> r;std::stringstream ss(v);std::string x;while(std::getline(ss,x,','))r.push_back(x);return r;};
  auto methods=list(o.methods),zmodes=list(o.zmode);std::vector<int> zsizes;for(auto &v:list(o.zsizes))zsizes.push_back(std::stoi(v));
  auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
  std::cerr<<"{\"event\":\"graph_loaded\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"ms\":"<<load_ms<<"}\n";
  std::unique_ptr<ZSkeleton> idx;double build_ms=0;
  for(auto &m:methods)if((m=="zskel"||m=="zskel-adaptive")&&!idx){auto b=Clock::now();idx=std::make_unique<ZSkeleton>(g,o.progress,o.hub_order,o.compress,nullptr,o.bp);configure(*idx,o);build_ms=elapsed(b);
    std::cerr<<"{\"event\":\"index_built\",\"ms\":"<<build_ms<<",\"entries\":"<<idx->labels.entries()<<",\"peak_rss\":"<<peak_rss()<<"}\n";}
  std::unique_ptr<ZSkeleton::Scratch> scratch;std::unique_ptr<AdaptiveZSkeleton> adaptive;
  if(idx){scratch=std::make_unique<ZSkeleton::Scratch>(g.n);adaptive=std::make_unique<AdaptiveZSkeleton>(*idx,g);}
  std::unique_ptr<CausalQueries<const ZSkeleton>> pure;std::unique_ptr<CausalQueries<AdaptiveZSkeleton>> ad;
  if(idx){pure=std::make_unique<CausalQueries<const ZSkeleton>>(*idx,*idx);ad=std::make_unique<CausalQueries<AdaptiveZSkeleton>>(*idx,*adaptive);}
  // Core–tree index (labels on the pendant-peeled core).
  std::unique_ptr<CoreTreeZSkeleton> ct;std::unique_ptr<CoreTreeZSkeleton::Scratch> cts;std::unique_ptr<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>> qc;
  for(auto &m:methods)if((m=="ctskel"||m=="ctskel-sf")&&!ct){auto b=Clock::now();ct=std::make_unique<CoreTreeZSkeleton>(g,o.progress,o.hub_order,o.compress,o.width,o.bp);configure(ct->inner,o);ct->closure_budget=o.closure_budget;ct->use_anc_masks=o.ancmask;ct->surgery_first=o.surgery_first;ct->prune_dominated=o.dominated;ct->parents_blocked_live=o.parent_live;
    cts=std::make_unique<CoreTreeZSkeleton::Scratch>(ct->make_scratch());qc=std::make_unique<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>>(*ct,*ct);
    if(!idx)build_ms=elapsed(b);
    std::cerr<<"{\"event\":\"index_built\",\"method\":\"ctskel\",\"ms\":"<<elapsed(b)<<",\"entries\":"<<ct->inner.labels.entries()<<",\"core_n\":"<<ct->core.n<<",\"peak_rss\":"<<peak_rss()<<"}\n";}
  TraversalCausal trav(g);
  // The four d-separation baselines with graph surgery and endpoint sets (src/surgery_baselines.hpp).
  std::unique_ptr<BaselineCausal<BidirectionalGeneral>> cb_bi;std::unique_ptr<BaselineCausal<AncestralGeneral>> cb_an;
  std::unique_ptr<BaselineCausal<CIflyGeneral>> cb_cf;std::unique_ptr<BaselineCausal<RPGeneral>> cb_rp;
  std::unique_ptr<BaselineCausal<BidirectionalGeneral>> cr_bi;std::unique_ptr<BaselineCausal<AncestralGeneral>> cr_an;
  std::unique_ptr<BaselineCausal<CIflyGeneral>> cr_cf;std::unique_ptr<BaselineCausal<RPGeneral>> cr_rp;
  std::unique_ptr<BaselineCausal<GeigerGeneral>> cb_ge;std::unique_ptr<BaselineCausal<DCollisionGeneral>> cb_dc;   // classical baselines
  for(auto &m:methods) {
    if(m=="geiger"&&!cb_ge)cb_ge=std::make_unique<BaselineCausal<GeigerGeneral>>(g);
    if(m=="dcollision"&&!cb_dc)cb_dc=std::make_unique<BaselineCausal<DCollisionGeneral>>(g);
    if(m=="bidirectional"&&!cb_bi)cb_bi=std::make_unique<BaselineCausal<BidirectionalGeneral>>(g);
    if(m=="bidirectional-anc"&&!cb_an)cb_an=std::make_unique<BaselineCausal<AncestralGeneral>>(g);
    if(m=="cifly"&&!cb_cf)cb_cf=std::make_unique<BaselineCausal<CIflyGeneral>>(g);
    if(m=="rp"&&!cb_rp)cb_rp=std::make_unique<BaselineCausal<RPGeneral>>(g);
    if(m=="bidirectional-r"&&!cr_bi){cr_bi=std::make_unique<BaselineCausal<BidirectionalGeneral>>(g);cr_bi->use_labels(reach_labels(g));}
    if(m=="bidirectional-anc-r"&&!cr_an){cr_an=std::make_unique<BaselineCausal<AncestralGeneral>>(g);cr_an->use_labels(reach_labels(g));}
    if(m=="cifly-r"&&!cr_cf){cr_cf=std::make_unique<BaselineCausal<CIflyGeneral>>(g);cr_cf->use_labels(reach_labels(g));}
    if(m=="rp-r"&&!cr_rp){cr_rp=std::make_unique<BaselineCausal<RPGeneral>>(g);cr_rp->use_labels(reach_labels(g));}
  }
  for(auto &kind:kinds)for(auto &zmode:zmodes)for(int k:zsizes) {
    if(zmode=="parents"&&kind!="backdoor")continue;
    if(zmode=="children"&&kind!="frontdoor")continue;
    // Requests: x, y from the pair mode; extra = k vertices, plus one for the rules
    // (z), mediation (m) and the sequential back-door (the second treatment).
    // Front-door: M = k vertices, or with Z mode "children" up to k children of x.
    // Sequential back-door: treatments (x, z0), W_1 and W_2 = halves of the rest.
    // Selection: the k vertices are the parents of the selection node.
    const bool plus1=kind!="backdoor"&&kind!="frontdoor"&&kind!="selection";
    Options w=o;w.zmode=zmode=="parents"||zmode=="children"?"random":zmode;w.zsize=k+(plus1?1:0);
    struct Q{std::vector<int> X,Y,Z,W;};
    auto build_queries=[&](const Options &wo){std::vector<Q> out;
      for(auto &r:generate(g,wo)) {
        Q q;q.X={r.x};q.Y={r.y};
        if(kind=="backdoor") {
          if(zmode=="parents"){for(int p:g.pa[r.x])if(p!=r.y&&int(q.Z.size())<k)q.Z.push_back(p);}
          else q.Z=r.z;
        } else if(kind=="frontdoor") {
          if(zmode=="children"){for(int c:g.ch[r.x])if(c!=r.y&&int(q.Z.size())<k)q.Z.push_back(c);}
          else q.Z=r.z;
        } else if(kind=="selection")q.Z=r.z;
        else if(kind=="seqbackdoor") {
          q.X={r.x,r.z[0]};const size_t half=1+(r.z.size()-1)/2;
          q.Z.assign(r.z.begin()+1,r.z.begin()+long(half));q.W.assign(r.z.begin()+long(half),r.z.end());
        } else {q.Z={r.z[0]};q.W.assign(r.z.begin()+1,r.z.end());}
        out.push_back(q);
      }
      return out;};
    std::vector<Q> qs=build_queries(w),wqs;
    if(o.warmup_disjoint){Options w2=w;w2.seed=o.seed+1000;wqs=build_queries(w2);}
    const std::vector<Q> &warm=o.warmup_disjoint?wqs:qs;
    std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"kind\":\""<<kind<<"\",\"zmode\":\""<<zmode
      <<"\",\"k\":"<<k<<",\"search\":\""<<(o.bidi?"bi":"uni")<<"\",\"live\":"<<o.live<<",\"balance\":\""<<(o.arcs?"arcs":"states")<<"\""<<",\"hub_order\":\""<<o.hub_order<<"\",\"surgery_first\":"<<o.surgery_first<<",\"pairmode\":\""<<o.pairmode<<"\",\"queries\":"<<qs.size()<<",\"warmup\":"<<o.warmup<<",\"build_ms\":"<<build_ms<<",\"methods\":[";
    std::vector<uint8_t> ref;bool first=true;
    for(auto &m:methods) {
      std::vector<uint8_t> ans(qs.size());std::vector<double> lat(qs.size());ZSkeletonStats st;
      int which=kind=="rule1"?1:kind=="rule2"?2:3;
      // "ctskel-sf": the same core-tree index with the surgery-first exploration start.
      if(ct)ct->surgery_first=m=="ctskel-sf"||(m=="ctskel"&&o.surgery_first);
      auto further=[&](auto &C,auto &S,const Q &q,ZSkeletonStats &st)->bool{
        if(kind=="frontdoor")return C.frontdoor(q.X,q.Y,q.Z,S,&st).valid;
        if(kind=="seqbackdoor")return C.sequential_backdoor(q.X,{q.Z,q.W},q.Y,S,&st);
        if(kind=="mediation")return C.mediation(q.X,q.Z,q.Y,q.W,S,&st);
        return C.selection_recoverable(q.X,q.Y,q.Z,S,&st);
      };
      const bool is_further=kind=="frontdoor"||kind=="seqbackdoor"||kind=="mediation"||kind=="selection";
      auto via_baseline=[&](auto &B,const Q &q)->bool{
        if(kind=="backdoor")return B.backdoor(q.X,q.Y,q.Z);
        if(kind=="frontdoor")return B.frontdoor(q.X,q.Y,q.Z);
        if(kind=="seqbackdoor")return B.sequential_backdoor(q.X,{q.Z,q.W},q.Y);
        if(kind=="mediation")return B.mediation(q.X,q.Z,q.Y,q.W);
        if(kind=="selection")return B.selection_recoverable(q.X,q.Y,q.Z);
        return B.rule(which,q.Y,q.Z,q.X,q.W);
      };
      auto answer=[&](const Q &q,ZSkeletonStats &st)->bool{
        if(m=="bidirectional")return via_baseline(*cb_bi,q);
        if(m=="bidirectional-anc")return via_baseline(*cb_an,q);
        if(m=="cifly")return via_baseline(*cb_cf,q);
        if(m=="rp")return via_baseline(*cb_rp,q);
        if(m=="geiger")return via_baseline(*cb_ge,q);
        if(m=="dcollision")return via_baseline(*cb_dc,q);
        if(m=="bidirectional-r")return via_baseline(*cr_bi,q);
        if(m=="bidirectional-anc-r")return via_baseline(*cr_an,q);
        if(m=="cifly-r")return via_baseline(*cr_cf,q);
        if(m=="rp-r")return via_baseline(*cr_rp,q);
        if(is_further) {
          if(m=="traversal") {
            if(kind=="frontdoor")return trav.frontdoor(q.X,q.Y,q.Z);
            if(kind=="seqbackdoor")return trav.sequential_backdoor(q.X,{q.Z,q.W},q.Y);
            if(kind=="mediation")return trav.mediation(q.X,q.Z,q.Y,q.W);
            return trav.selection_recoverable(q.X,q.Y,q.Z);
          }
          if(m=="ctskel"||m=="ctskel-sf")return further(*qc,*cts,q,st);
          return m=="zskel"?further(*pure,*scratch,q,st):further(*ad,*scratch,q,st);
        }
        if(m=="ctskel"||m=="ctskel-sf") {
          if(kind=="backdoor")return qc->backdoor(q.X,q.Y,q.Z,*cts,&st).valid;
          return which==1?qc->rule1(q.Y,q.Z,q.X,q.W,*cts,&st):which==2?qc->rule2(q.Y,q.Z,q.X,q.W,*cts,&st):qc->rule3(q.Y,q.Z,q.X,q.W,*cts,&st);
        }
        if(kind=="backdoor") return m=="traversal"?trav.backdoor(q.X,q.Y,q.Z).valid:m=="zskel"?pure->backdoor(q.X,q.Y,q.Z,*scratch,&st).valid:ad->backdoor(q.X,q.Y,q.Z,*scratch,&st).valid;
        if(m=="traversal")return trav.rule(which,q.Y,q.Z,q.X,q.W);
        return m=="zskel"?(which==1?pure->rule1(q.Y,q.Z,q.X,q.W,*scratch,&st):which==2?pure->rule2(q.Y,q.Z,q.X,q.W,*scratch,&st):pure->rule3(q.Y,q.Z,q.X,q.W,*scratch,&st))
                         :(which==1?ad->rule1(q.Y,q.Z,q.X,q.W,*scratch,&st):which==2?ad->rule2(q.Y,q.Z,q.X,q.W,*scratch,&st):ad->rule3(q.Y,q.Z,q.X,q.W,*scratch,&st));};
      // Untimed warm-up passes over the same requests, as in bench.
      for(int w=0;w<o.warmup;++w){ZSkeletonStats ws;for(auto &q:warm)answer(q,ws);}
      phase_times().reset();
      auto q0=Clock::now();
      for(size_t i=0;i<qs.size();++i) {
        auto s0=Clock::now();
        ans[i]=answer(qs[i],st);lat[i]=std::chrono::duration<double,std::micro>(Clock::now()-s0).count();
      }
      double query_ms=elapsed(q0);size_t mism=0,yes=0;
      if(ref.empty())ref=ans;
      for(size_t i=0;i<qs.size();++i){mism+=ans[i]!=ref[i];yes+=ans[i];}
      std::cout<<(first?"":",")<<"{\"method\":\""<<m<<"\",\"query_ms\":"<<query_ms<<",\"mean_us\":"<<(qs.empty()?0:1000*query_ms/double(qs.size()))
        <<",\"p50_us\":"<<pct(lat,.5)<<",\"p90_us\":"<<pct(lat,.9)<<",\"p99_us\":"<<pct(lat,.99)<<",\"max_us\":"<<pct(lat,1)
        <<",\"holds\":"<<yes<<",\"mismatches_vs_first\":"<<mism<<",\"label_build_ms\":"<<(plus_r(m)?reach_labels(g).build_ms:0)<<",\"sep_ms\":"<<phase_times().sep_ns/1e6<<",\"anc_ms\":"<<phase_times().anc_ns/1e6<<",\"adaptive_index\":"<<st.adaptive_index<<",\"adaptive_traversal\":"<<st.adaptive_traversal
        <<",\"h_nodes\":"<<st.h_nodes<<",\"edge_checks\":"<<st.edge_checks<<",\"neg_reach\":"<<st.neg_reach<<",\"cert_free\":"<<st.cert_free<<",\"searches\":"<<st.searches
        <<",\"search_found\":"<<st.search_found<<",\"search_states\":"<<st.search_states<<",\"label_merges\":"<<st.label_merges<<",\"closures\":"<<st.closures<<",\"core_routes\":"<<st.core_routes<<",\"mask_entries\":"<<st.mask_entries<<"}";
      first=false;
      if(mism)std::cerr<<"{\"event\":\"MISMATCH\",\"kind\":\""<<kind<<"\",\"method\":\""<<m<<"\",\"count\":"<<mism<<"}\n";
      std::cerr<<"{\"event\":\"method_done\",\"kind\":\""<<kind<<"\",\"zmode\":\""<<zmode<<"\",\"k\":"<<k<<",\"method\":\""<<m<<"\",\"query_ms\":"<<query_ms<<"}\n";
    }
    std::cout<<"],\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
  }
  return 0;
}

// End-to-end workload: literature-scale back-door screening. Hypotheses
// (x, y) have y a descendant of x (random downward walk, 1..max_dist steps).
// A seeded fraction of vertices is "observed" (cohort covariates). For each
// hypothesis: try Z1 = observed parents of x (<= cap); if invalid, try
// Z2 = Z1 plus observed parents of y minus x, y (<= cap). Report validity,
// per-method end-to-end time including index build, and agreement.
int screen_main(const Options &o,size_t hypotheses,double observed_frac,int max_dist,size_t cap) {
  auto list=[](const std::string &v){std::vector<std::string> r;std::stringstream ss(v);std::string x;while(std::getline(ss,x,','))r.push_back(x);return r;};
  auto methods=list(o.methods);
  auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
  std::cerr<<"{\"event\":\"graph_loaded\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"ms\":"<<load_ms<<"}\n";
  std::mt19937_64 rng(o.seed);
  std::vector<uint8_t> observed(g.n);std::bernoulli_distribution obs(observed_frac);
  for(auto &v:observed)v=obs(rng);
  std::vector<int> sources;for(int v=0;v<g.n;++v)if(!g.ch[v].empty())sources.push_back(v);
  struct H{int x,y;std::vector<int> z1,z2;};std::vector<H> hs;hs.reserve(hypotheses);
  while(hs.size()<hypotheses) {
    H h;h.x=sources[rng()%sources.size()];int v=h.x,len=1+int(rng()%max_dist);
    for(int i=0;i<len&&!g.ch[v].empty();++i)v=g.ch[v][rng()%g.ch[v].size()];
    if(v==h.x)continue;
    h.y=v;
    for(int p:g.pa[h.x])if(observed[p]&&p!=h.y&&h.z1.size()<cap)h.z1.push_back(p);
    h.z2=h.z1;std::unordered_set<int> in(h.z1.begin(),h.z1.end());
    for(int p:g.pa[h.y])if(observed[p]&&p!=h.x&&!in.count(p)&&h.z2.size()<cap){h.z2.push_back(p);in.insert(p);}
    hs.push_back(std::move(h));
  }
  std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"hypotheses\":"<<hs.size()
    <<",\"observed_frac\":"<<observed_frac<<",\"max_dist\":"<<max_dist<<",\"cap\":"<<cap<<",\"seed\":"<<o.seed<<",\"methods\":[";
  std::vector<uint8_t> ref;bool first=true;
  for(auto &m:methods) {
    double build_ms=0;std::unique_ptr<ZSkeleton> idx;std::unique_ptr<ZSkeleton::Scratch> sc;std::unique_ptr<AdaptiveZSkeleton> ad;
    std::unique_ptr<CausalQueries<AdaptiveZSkeleton>> qa;std::unique_ptr<CausalQueries<const ZSkeleton>> qp;std::unique_ptr<TraversalCausal> tr;
    std::unique_ptr<CoreTreeZSkeleton> ct;std::unique_ptr<CoreTreeZSkeleton::Scratch> cts;std::unique_ptr<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>> qc;
    std::unique_ptr<BackdoorBaseline<BidirectionalCIfly>> bb;std::unique_ptr<BackdoorBaseline<AncestralBidirectional>> ba;
    std::unique_ptr<BackdoorBaseline<StreamingCIfly>> bc;std::unique_ptr<BackdoorBaseline<RPReachable>> br;
    std::unique_ptr<BackdoorBaseline<GeigerReachable>> bg;std::unique_ptr<BackdoorBaseline<DCollision>> bd;
    auto b=Clock::now();
    if(m=="traversal")tr=std::make_unique<TraversalCausal>(g);
    else if(m=="geiger")bg=std::make_unique<BackdoorBaseline<GeigerReachable>>(g);
    else if(m=="dcollision")bd=std::make_unique<BackdoorBaseline<DCollision>>(g,false);
    else if(m=="bidirectional")bb=std::make_unique<BackdoorBaseline<BidirectionalCIfly>>(g);
    else if(m=="bidirectional-anc")ba=std::make_unique<BackdoorBaseline<AncestralBidirectional>>(g);
    else if(m=="cifly")bc=std::make_unique<BackdoorBaseline<StreamingCIfly>>(g);
    else if(m=="rp")br=std::make_unique<BackdoorBaseline<RPReachable>>(g,true);
    else if(m=="bidirectional-r"){bb=std::make_unique<BackdoorBaseline<BidirectionalCIfly>>(g);bb->use_labels(reach_labels(g));}
    else if(m=="bidirectional-anc-r"){ba=std::make_unique<BackdoorBaseline<AncestralBidirectional>>(g);ba->use_labels(reach_labels(g));}
    else if(m=="cifly-r"){bc=std::make_unique<BackdoorBaseline<StreamingCIfly>>(g);bc->use_labels(reach_labels(g));}
    else if(m=="rp-r"){br=std::make_unique<BackdoorBaseline<RPReachable>>(g,true);br->use_labels(reach_labels(g));}
    else if(m=="zskel"||m=="zskel-adaptive") {
      idx=std::make_unique<ZSkeleton>(g,o.progress,o.hub_order,o.compress,nullptr,o.bp);configure(*idx,o);sc=std::make_unique<ZSkeleton::Scratch>(g.n);
      if(m=="zskel"){qp=std::make_unique<CausalQueries<const ZSkeleton>>(*idx,*idx);}
      else {ad=std::make_unique<AdaptiveZSkeleton>(*idx,g);qa=std::make_unique<CausalQueries<AdaptiveZSkeleton>>(*idx,*ad);}
    } else if(m=="ctskel") {
      ct=std::make_unique<CoreTreeZSkeleton>(g,o.progress,o.hub_order,o.compress,o.width,o.bp);configure(ct->inner,o);ct->closure_budget=o.closure_budget;ct->use_anc_masks=o.ancmask;ct->surgery_first=o.surgery_first;ct->prune_dominated=o.dominated;ct->parents_blocked_live=o.parent_live;
      cts=std::make_unique<CoreTreeZSkeleton::Scratch>(ct->make_scratch());qc=std::make_unique<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>>(*ct,*ct);
    } else throw std::invalid_argument("screen method "+m);
    build_ms=plus_r(m)?reach_labels(g).build_ms:elapsed(b);
    // outcome per hypothesis: 0 none valid, 1 Z1 valid, 2 Z2 valid
    std::vector<uint8_t> out(hs.size());size_t tests=0;
    phase_times().reset();
    auto q0=Clock::now();
    for(size_t i=0;i<hs.size();++i) {
      auto &h=hs[i];
      auto test=[&](const std::vector<int> &z){++tests;
        if(tr)return tr->backdoor({h.x},{h.y},z).valid;
        if(bb)return bb->valid({h.x},{h.y},z);
        if(ba)return ba->valid({h.x},{h.y},z);
        if(bc)return bc->valid({h.x},{h.y},z);
        if(br)return br->valid({h.x},{h.y},z);
        if(bg)return bg->valid({h.x},{h.y},z);
        if(bd)return bd->valid({h.x},{h.y},z);
        if(qp)return qp->backdoor({h.x},{h.y},z,*sc).valid;
        if(qc)return qc->backdoor({h.x},{h.y},z,*cts).valid;
        return qa->backdoor({h.x},{h.y},z,*sc).valid;};
      out[i]=test(h.z1)?1:(h.z2.size()>h.z1.size()&&test(h.z2))?2:0;
    }
    double query_ms=elapsed(q0);
    size_t mism=0,c1=0,c2=0;if(ref.empty())ref=out;
    for(size_t i=0;i<hs.size();++i){mism+=out[i]!=ref[i];c1+=out[i]==1;c2+=out[i]==2;}
    std::cout<<(first?"":",")<<"{\"method\":\""<<m<<"\",\"build_ms\":"<<build_ms<<",\"query_ms\":"<<query_ms
      <<",\"total_ms\":"<<build_ms+query_ms<<",\"tests\":"<<tests<<",\"us_per_test\":"<<1000*query_ms/double(std::max<size_t>(tests,1))
      <<",\"valid_z1\":"<<c1<<",\"valid_z2\":"<<c2<<",\"mismatches_vs_first\":"<<mism<<",\"sep_ms\":"<<phase_times().sep_ns/1e6<<",\"anc_ms\":"<<phase_times().anc_ns/1e6<<"}";
    first=false;
    std::cerr<<"{\"event\":\"screen_done\",\"method\":\""<<m<<"\",\"build_ms\":"<<build_ms<<",\"query_ms\":"<<query_ms<<"}\n";
  }
  std::cout<<"],\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
  return 0;
}

// Replay a file of real back-door checks: one request per line,
//   B <nx> x.. <ny> y.. <nz> z..
// (back-door validity of Z for X -> Y). Used for workloads taken from
// published studies (literature replay, paper Section 7.4).
int replay_main(const Options &o,const std::string &requests) {
  auto list=[](const std::string &v){std::vector<std::string> r;std::stringstream ss(v);std::string x;while(std::getline(ss,x,','))r.push_back(x);return r;};
  auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
  struct Req{std::vector<int> X,Y,Z;};std::vector<Req> rs;
  {std::ifstream in(requests);std::string kind;
   auto vec=[&](std::vector<int> &v){size_t k;in>>k;v.resize(k);for(auto &x:v)in>>x;};
   while(in>>kind){if(kind!="B")throw std::invalid_argument("replay kind "+kind);Req r;vec(r.X);vec(r.Y);vec(r.Z);rs.push_back(std::move(r));}}
  std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"load_ms\":"<<load_ms
    <<",\"requests\":\""<<requests<<"\",\"count\":"<<rs.size()<<",\"warmup\":"<<o.warmup<<",\"methods\":[";
  std::vector<uint8_t> ref;bool first=true;
  for(auto &m:list(o.methods)) {
    double build_ms=0;std::unique_ptr<ZSkeleton> idx;std::unique_ptr<ZSkeleton::Scratch> sc;std::unique_ptr<AdaptiveZSkeleton> ad;
    std::unique_ptr<CausalQueries<AdaptiveZSkeleton>> qa;std::unique_ptr<CausalQueries<const ZSkeleton>> qp;std::unique_ptr<TraversalCausal> tr;
    std::unique_ptr<CoreTreeZSkeleton> ct;std::unique_ptr<CoreTreeZSkeleton::Scratch> cts;std::unique_ptr<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>> qc;
    std::unique_ptr<BackdoorBaseline<BidirectionalCIfly>> bb;std::unique_ptr<BackdoorBaseline<AncestralBidirectional>> ba;
    std::unique_ptr<BackdoorBaseline<StreamingCIfly>> bc;std::unique_ptr<BackdoorBaseline<RPReachable>> br;
    std::unique_ptr<BackdoorBaseline<GeigerReachable>> bg;std::unique_ptr<BackdoorBaseline<DCollision>> bd;
    auto b=Clock::now();
    if(m=="traversal")tr=std::make_unique<TraversalCausal>(g);
    else if(m=="geiger")bg=std::make_unique<BackdoorBaseline<GeigerReachable>>(g);
    else if(m=="dcollision")bd=std::make_unique<BackdoorBaseline<DCollision>>(g,false);
    else if(m=="bidirectional")bb=std::make_unique<BackdoorBaseline<BidirectionalCIfly>>(g);
    else if(m=="bidirectional-anc")ba=std::make_unique<BackdoorBaseline<AncestralBidirectional>>(g);
    else if(m=="cifly")bc=std::make_unique<BackdoorBaseline<StreamingCIfly>>(g);
    else if(m=="rp")br=std::make_unique<BackdoorBaseline<RPReachable>>(g,true);
    else if(m=="bidirectional-r"){bb=std::make_unique<BackdoorBaseline<BidirectionalCIfly>>(g);bb->use_labels(reach_labels(g));}
    else if(m=="bidirectional-anc-r"){ba=std::make_unique<BackdoorBaseline<AncestralBidirectional>>(g);ba->use_labels(reach_labels(g));}
    else if(m=="cifly-r"){bc=std::make_unique<BackdoorBaseline<StreamingCIfly>>(g);bc->use_labels(reach_labels(g));}
    else if(m=="rp-r"){br=std::make_unique<BackdoorBaseline<RPReachable>>(g,true);br->use_labels(reach_labels(g));}
    else if(m=="zskel"||m=="zskel-adaptive") {
      idx=std::make_unique<ZSkeleton>(g,o.progress,o.hub_order,o.compress,nullptr,o.bp);configure(*idx,o);sc=std::make_unique<ZSkeleton::Scratch>(g.n);
      if(m=="zskel"){qp=std::make_unique<CausalQueries<const ZSkeleton>>(*idx,*idx);}
      else {ad=std::make_unique<AdaptiveZSkeleton>(*idx,g);ad->probe_budget=o.probe;ad->large_h=o.large_h;ad->initial_budget=o.initial_budget;qa=std::make_unique<CausalQueries<AdaptiveZSkeleton>>(*idx,*ad);}
    } else if(m=="ctskel") {
      ct=std::make_unique<CoreTreeZSkeleton>(g,o.progress,o.hub_order,o.compress,o.width,o.bp);configure(ct->inner,o);ct->closure_budget=o.closure_budget;ct->use_anc_masks=o.ancmask;ct->surgery_first=o.surgery_first;ct->prune_dominated=o.dominated;ct->parents_blocked_live=o.parent_live;
      cts=std::make_unique<CoreTreeZSkeleton::Scratch>(ct->make_scratch());qc=std::make_unique<CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton>>(*ct,*ct);
    } else throw std::invalid_argument("replay method "+m);
    build_ms=plus_r(m)?reach_labels(g).build_ms:elapsed(b);
    auto run=[&](const Req &r){return bg?bg->valid(r.X,r.Y,r.Z):bd?bd->valid(r.X,r.Y,r.Z):bb?bb->valid(r.X,r.Y,r.Z):ba?ba->valid(r.X,r.Y,r.Z):bc?bc->valid(r.X,r.Y,r.Z):br?br->valid(r.X,r.Y,r.Z):tr?tr->backdoor(r.X,r.Y,r.Z).valid:qp?qp->backdoor(r.X,r.Y,r.Z,*sc).valid:qc?qc->backdoor(r.X,r.Y,r.Z,*cts).valid:qa->backdoor(r.X,r.Y,r.Z,*sc).valid;};
    for(int w=0;w<o.warmup;++w)for(auto &r:rs)run(r);
    std::vector<uint8_t> ans(rs.size());std::vector<double> lat(rs.size());
    phase_times().reset();
    auto q0=Clock::now();
    for(size_t i=0;i<rs.size();++i){auto s=Clock::now();ans[i]=run(rs[i]);lat[i]=std::chrono::duration<double,std::micro>(Clock::now()-s).count();}
    double query_ms=elapsed(q0);size_t mism=0,valid=0;if(ref.empty())ref=ans;
    for(size_t i=0;i<rs.size();++i){mism+=ans[i]!=ref[i];valid+=ans[i];}
    std::cout<<(first?"":",")<<"{\"method\":\""<<m<<"\",\"build_ms\":"<<build_ms<<",\"query_ms\":"<<query_ms
      <<",\"mean_us\":"<<1000*query_ms/double(std::max<size_t>(rs.size(),1))<<",\"p50_us\":"<<pct(lat,.5)<<",\"p99_us\":"<<pct(lat,.99)
      <<",\"valid\":"<<valid<<",\"mismatches_vs_first\":"<<mism<<",\"sep_ms\":"<<phase_times().sep_ns/1e6<<",\"anc_ms\":"<<phase_times().anc_ns/1e6<<"}";
    first=false;
    std::cerr<<"{\"event\":\"replay_done\",\"method\":\""<<m<<"\",\"query_ms\":"<<query_ms<<"}\n";
  }
  std::cout<<"],\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
  return 0;
}

// Component-decomposed index for very large, fragmented graphs (CausalBank:
// 580M vertices in 268M weakly connected components). Vertices in different
// WCCs are d-separated for every Z (conditioning never creates skeleton
// paths), and Z outside the endpoints' WCC is irrelevant. So: one index per
// WCC of size >= threshold (built single-threaded, one after another); exact
// bidirectional traversal for queries inside smaller WCCs; constant answer
// across WCCs. Baseline: bidirectional traversal on the full graph.
int components_main(const Options &o,uint32_t threshold) {
  auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
  std::cerr<<"{\"event\":\"graph_loaded\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"ms\":"<<load_ms<<"}\n";
  auto t1=Clock::now();
  std::vector<uint32_t> parent(g.n);std::iota(parent.begin(),parent.end(),0u);
  auto find=[&](uint32_t x){while(parent[x]!=x){parent[x]=parent[parent[x]];x=parent[x];}return x;};
  for(int v=0;v<g.n;++v)for(int c:g.ch[v]){auto a=find(v),b=find(c);if(a!=b)parent[a]=b;}
  std::vector<uint32_t> comp(g.n),size_of(g.n,0);
  for(int v=0;v<g.n;++v){comp[v]=find(v);++size_of[comp[v]];}
  std::vector<uint32_t>().swap(parent);
  std::vector<uint32_t> big;for(int v=0;v<g.n;++v)if(comp[v]==uint32_t(v)&&size_of[v]>=threshold)big.push_back(v);
  std::sort(big.begin(),big.end(),[&](uint32_t a,uint32_t b){return size_of[a]>size_of[b];});
  uint64_t wcc=0,nonisolated=0;for(int v=0;v<g.n;++v){if(comp[v]==uint32_t(v))++wcc;if(size_of[comp[v]]>1)++nonisolated;}
  double wcc_ms=elapsed(t1);
  std::cerr<<"{\"event\":\"components\",\"wcc\":"<<wcc<<",\"big\":"<<big.size()<<",\"largest\":"<<(big.empty()?0:size_of[big[0]])<<",\"ms\":"<<wcc_ms<<"}\n";
  // Per big component: local graph, index, adaptive engine.
  struct Part{uint32_t root;std::vector<int> members;std::unique_ptr<Graph> g;std::unique_ptr<ZSkeleton> idx;
    std::unique_ptr<ZSkeleton::Scratch> sc;std::unique_ptr<AdaptiveZSkeleton> ad;double build_ms=0;};
  std::vector<Part> parts(big.size());std::unordered_map<uint32_t,uint32_t> part_of;
  std::vector<int> local(g.n,-1);
  for(size_t i=0;i<big.size();++i)part_of[big[i]]=uint32_t(i);
  for(int v=0;v<g.n;++v){auto it=part_of.find(comp[v]);if(it!=part_of.end()){auto &P=parts[it->second];local[v]=int(P.members.size());P.members.push_back(v);}}
  double build_ms=0;size_t index_bytes=0,entries=0;
  for(size_t i=0;i<parts.size();++i) {
    auto &P=parts[i];P.root=big[i];auto b=Clock::now();
    std::vector<std::pair<int,int>> es;for(int v:P.members)for(int c:g.ch[v])es.emplace_back(local[v],local[c]);
    P.g=std::make_unique<Graph>(int(P.members.size()),es);
    P.idx=std::make_unique<ZSkeleton>(*P.g,o.progress&&i==0,o.hub_order,o.compress,nullptr,o.bp);configure(*P.idx,o);
    P.sc=std::make_unique<ZSkeleton::Scratch>(P.g->n);
    P.ad=std::make_unique<AdaptiveZSkeleton>(*P.idx,*P.g);P.ad->initial_budget=o.initial_budget;P.ad->large_h=o.large_h;
    P.build_ms=elapsed(b);build_ms+=P.build_ms;index_bytes+=P.idx->labels.bytes();entries+=P.idx->labels.entries();
    if(i<5||i+1==parts.size())std::cerr<<"{\"event\":\"part_built\",\"i\":"<<i<<",\"n\":"<<P.g->n<<",\"m\":"<<P.g->m<<",\"ms\":"<<P.build_ms<<"}\n";
  }
  BidirectionalCIfly full(g);   // baseline on the whole graph, and small-WCC path
  // Workloads: x uniform over non-isolated vertices ("all") or over the
  // largest WCC ("largest"); y by a random skeleton walk from x; Z near /
  // ancestor / random within x's WCC.
  std::mt19937_64 rng(o.seed);
  std::vector<int> pool_all;pool_all.reserve(nonisolated);for(int v=0;v<g.n;++v)if(size_of[comp[v]]>1)pool_all.push_back(v);
  auto step=[&](int v,bool up){size_t a=g.pa[v].size(),b=up?0:g.ch[v].size();if(!a&&!b)return v;auto r=rng()%(a+b);return r<a?g.pa[v][r]:g.ch[v][r-a];};
  auto walk=[&](int v,int len,bool up){for(int i=0;i<len;++i)v=step(v,up);return v;};
  std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"load_ms\":"<<load_ms
    <<",\"wcc\":"<<wcc<<",\"wcc_ms\":"<<wcc_ms<<",\"threshold\":"<<threshold<<",\"indexed_components\":"<<parts.size()
    <<",\"largest_component\":"<<(big.empty()?0:size_of[big[0]])<<",\"indexed_vertices\":";
  {uint64_t iv=0;for(auto &P:parts)iv+=P.members.size();std::cout<<iv;}
  std::cout<<",\"build_ms\":"<<build_ms<<",\"index_bytes\":"<<index_bytes<<",\"label_entries\":"<<entries<<",\"runs\":[";
  auto list=[](const std::string &v){std::vector<std::string> r;std::stringstream ss(v);std::string x;while(std::getline(ss,x,','))r.push_back(x);return r;};
  bool firstrun=true;
  for(const std::string pool_name:{"all","largest"})for(auto &zmode:list(o.zmode))for(auto &kz:list(o.zsizes)) {
    int k=std::stoi(kz);
    struct Q{int x,y;std::vector<int> z;};std::vector<Q> qs;
    while(qs.size()<o.queries) {
      int x=pool_name=="all"?pool_all[rng()%pool_all.size()]:(big.empty()?pool_all[0]:parts[0].members[rng()%parts[0].members.size()]);
      int y=x;for(int tr=0;tr<16&&y==x;++tr)y=walk(x,1+int(rng()%6),false);
      if(y==x)continue;
      std::unordered_set<int> used{x,y};std::vector<int> z;int budget=64*k+64;
      const auto *members=part_of.count(comp[x])?&parts[part_of[comp[x]]].members:nullptr;
      while(int(z.size())<k&&budget-->0){
        int v;   // walk length first, then the start vertex (g++ argument order, made explicit)
        if(zmode=="near"){int len=1+int(rng()%4);int from=rng()%2?x:y;v=walk(from,len,false);}
        else if(zmode=="ancestor"){int len=1+int(rng()%8);int from=rng()%2?x:y;v=walk(from,len,true);}
        else v=members?(*members)[rng()%members->size()]:walk(x,1+int(rng()%6),false);
        if(used.insert(v).second)z.push_back(v);
      }
      qs.push_back({x,y,z});
    }
    std::vector<uint8_t> ref(qs.size()),ans(qs.size());std::vector<double> lb(qs.size()),li(qs.size());
    size_t cross=0,small=0,indexed=0;
    for(int w=0;w<=o.warmup;++w) {
      const bool timed=w==o.warmup;
      for(size_t i=0;i<qs.size();++i){auto s=Clock::now();ref[i]=full.separated(qs[i].x,qs[i].y,qs[i].z);if(timed)lb[i]=std::chrono::duration<double,std::micro>(Clock::now()-s).count();}
      for(size_t i=0;i<qs.size();++i) {
        auto &q=qs[i];auto s=Clock::now();
        if(comp[q.x]!=comp[q.y]){ans[i]=1;if(timed)++cross;}
        else{auto it=part_of.find(comp[q.x]);
          if(it==part_of.end()){ans[i]=full.separated(q.x,q.y,q.z);if(timed)++small;}
          else{auto &P=parts[it->second];std::vector<int> zl;for(int v:q.z)if(comp[v]==comp[q.x])zl.push_back(local[v]);
            ans[i]=P.ad->separated(local[q.x],local[q.y],zl,*P.sc);if(timed)++indexed;}}
        if(timed)li[i]=std::chrono::duration<double,std::micro>(Clock::now()-s).count();
      }
    }
    size_t mism=0,sep=0;for(size_t i=0;i<qs.size();++i){mism+=ans[i]!=ref[i];sep+=ans[i];}
    double mb=0,mi=0;for(size_t i=0;i<qs.size();++i){mb+=lb[i];mi+=li[i];}
    mb/=double(qs.size());mi/=double(qs.size());
    std::cout<<(firstrun?"":",")<<"{\"pool\":\""<<pool_name<<"\",\"zmode\":\""<<zmode<<"\",\"k\":"<<k<<",\"queries\":"<<qs.size()
      <<",\"bidirectional_mean_us\":"<<mb<<",\"bidirectional_p50_us\":"<<pct(lb,.5)<<",\"bidirectional_p99_us\":"<<pct(lb,.99)
      <<",\"decomposed_mean_us\":"<<mi<<",\"decomposed_p50_us\":"<<pct(li,.5)<<",\"decomposed_p99_us\":"<<pct(li,.99)
      <<",\"answered_indexed\":"<<indexed<<",\"answered_small_wcc\":"<<small<<",\"answered_cross\":"<<cross
      <<",\"separated\":"<<sep<<",\"mismatches\":"<<mism<<"}";
    firstrun=false;std::cout.flush();
    std::cerr<<"{\"event\":\"run_done\",\"pool\":\""<<pool_name<<"\",\"zmode\":\""<<zmode<<"\",\"k\":"<<k<<",\"bi\":"<<mb<<",\"dec\":"<<mi<<",\"mism\":"<<mism<<"}\n";
  }
  std::cout<<"],\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
  return 0;
}

// Index statistics: build TrekIndex once and print the quantities of the complexity
// analysis in Sections 5.2 and 6.4: fill edges f, shortcut arcs created m_s, core shortcut arcs, the
// largest number of parallel shortcut arcs between two states, delta (shortcut arcs at a vertex), the
// chain length c, the label size l, and the witness sizes. Not a timed run.
int stats_main(const Options &o) {
  auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
  auto b=Clock::now();
  CoreTreeZSkeleton ct(g,o.progress,o.hub_order,o.compress,o.width,o.bp);
  const double build_ms=elapsed(b);
  // parallel multiplicity among all shortcut arcs created, and among the core shortcut arcs
  auto maxmult=[](auto &&keys){std::unordered_map<uint64_t,uint32_t> c;uint32_t mx=0;for(uint64_t k:keys)mx=std::max(mx,++c[k]);return mx;};
  std::vector<uint64_t> ka;ka.reserve(ct.arcs.size());for(auto &a:ct.arcs)ka.push_back(uint64_t(a.from)<<32|a.to);
  std::vector<uint64_t> kc;for(size_t st=0;st<ct.sc.out.size();++st)for(auto &e:ct.sc.out[st])kc.push_back(uint64_t(st)<<32|e.to);
  const uint32_t mult_all=maxmult(ka),mult_core=maxmult(kc);
  // delta: local shortcut arcs of a peeled vertex; shortcut arcs at the four states of a core vertex
  size_t dp=0,dc=0;const size_t np=ct.lsc_off.empty()?0:ct.lsc_off.size()-1;
  for(size_t q=0;q<np;++q)dp=std::max<size_t>(dp,ct.lsc_off[q+1]-ct.lsc_off[q]);
  for(size_t v=0;v*4<ct.sc.out.size();++v){size_t d=0;for(int k=0;k<4;++k)d+=ct.sc.out[4*v+k].size()+ct.sc.in[4*v+k].size();dc=std::max(dc,d);}
  // chain length c
  size_t cmax=0;double csum=0;for(size_t q=0;q<np;++q){size_t c=ct.ch_off[q+1]-ct.ch_off[q];cmax=std::max(cmax,c);csum+=double(c);}
  // label size l on the witness DAG
  const auto &L=ct.inner.labels;const size_t S=L.out_off.empty()?0:L.out_off.size()-1;
  uint64_t lmax=0;double lsum=0;for(size_t st=0;st<S;++st){uint64_t a=L.out_len(uint32_t(st)),c=L.in_len(uint32_t(st));lmax=std::max({lmax,a,c});lsum+=double(a+c);}
  // witness sizes (states) of the core shortcut arcs
  std::function<size_t(int32_t)> wsize=[&](int32_t w)->size_t{size_t r=0;while(w>=0){const auto &x=ct.wits[w];r+=size_t(__builtin_popcount(x.mask));if(x.right>=0)r+=wsize(x.right);w=x.left;}return r;};
  size_t wmax=0;double wsum=0;size_t wn=0;
  for(size_t st=0;st<ct.sc.out.size();++st)for(auto &e:ct.sc.out[st]){size_t z=wsize(ct.arcs[e.arc].wit);wmax=std::max(wmax,z);wsum+=double(z);++wn;}
  std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"load_ms\":"<<load_ms<<",\"build_ms\":"<<build_ms
    <<",\"width\":"<<ct.width<<",\"peeled\":"<<ct.peeled<<",\"core_n\":"<<ct.core.n<<",\"core_m\":"<<ct.core.m<<",\"fill_edges\":"<<ct.fill_edges
    <<",\"shortcuts_created\":"<<ct.arcs.size()<<",\"core_shortcuts\":"<<ct.sc.count<<",\"max_parallel_all\":"<<mult_all<<",\"max_parallel_core\":"<<mult_core
    <<",\"delta_peeled\":"<<dp<<",\"delta_core\":"<<dc<<",\"chain_max\":"<<cmax<<",\"chain_mean\":"<<(np?csum/double(np):0)<<",\"max_depth\":"<<ct.max_depth
    <<",\"label_entries\":"<<L.entries()<<",\"label_max\":"<<lmax<<",\"label_mean\":"<<(S?lsum/double(2*S):0)
    <<",\"witness_max\":"<<wmax<<",\"witness_mean\":"<<(wn?wsum/double(wn):0)<<",\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
  return 0;
}

// Causal-discovery workload (PC tests of Exp-8): the conditional-independence tests that the skeleton phase
// of PC-stable (also the first phase of FCI) issues when it runs with a d-separation oracle over a
// set S of observed variables of the known DAG G. A region S is either An(v) for a random v with
// min_size <= |An(v)| <= region_size ("ancestral", the default: causally sufficient, so PC recovers
// the skeleton of the induced DAG), or the first region_size vertices of a breadth-first search in the
// skeleton ("ball": the other variables are latent; the oracle still answers d-separation in G).
// Defaults: 20 ancestral regions of 20 to 50 vertices, level 1 (the paper: --min-size 800 --region-size 1000).
// Level l tests every edge x - y still present against every subset
// Z of size l of adj(x) \ {y}, then of adj(y) \ {x}, with the adjacencies frozen at the start of the
// level (PC-stable), up to level L or until no edge has enough neighbours; a region stops after
// max_tests tests. Writes the tests as a workload file (read_workload) and prints statistics.
int pcgen_main(const Options &o,int regions,int region_size,int max_level,size_t max_tests,const std::string &region_kind,int min_size) {
  auto g=read_binary(o.graph);
  AncestralBidirectional oracle(g);   // exact and fast on ancestral regions; the tests do not depend on the oracle
  std::mt19937_64 rng(o.seed);
  std::vector<Request> out;std::vector<size_t> per_level(max_level+1,0),sep_level(max_level+1,0);size_t capped=0,edges_left=0;
  std::vector<int> mark(g.n,-1);
  auto t=Clock::now();
  for(int r=0;r<regions;++r) {
    // region: BFS ball of region_size vertices in the skeleton
    std::vector<int> S;
    if(region_kind=="ancestral") {
      // S = An(v) for a random v with min_size <= |An(v)| <= region_size: causally sufficient, so the
      // oracle answers d-separation in the induced DAG and PC recovers its skeleton.
      for(int tries=0;tries<200000&&S.empty();++tries) {
        int v=int(rng()%uint64_t(g.n));const int tag=-(r*200000+tries)-2;
        std::vector<int> q{v};mark[v]=tag;bool big=false;
        for(size_t i=0;i<q.size()&&!big;++i)for(int p:g.pa[q[i]])if(mark[p]!=tag){mark[p]=tag;q.push_back(p);if(int(q.size())>region_size){big=true;break;}}
        if(!big&&int(q.size())>=min_size)S=q;
      }
      if(S.empty())throw std::runtime_error("no ancestral region of the requested size");
    }
    for(int tries=0;tries<1000&&int(S.size())<region_size&&region_kind!="ancestral";++tries) {
      S.clear();int root=int(rng()%uint64_t(g.n));if(g.pa[root].empty()&&g.ch[root].empty())continue;
      std::vector<int> q{root};mark[root]=r*1000+tries;
      for(size_t i=0;i<q.size()&&int(q.size())<region_size;++i){int v=q[i];
        for(auto *L:{&g.pa[v],&g.ch[v]})for(int u:*L)if(mark[u]!=r*1000+tries&&int(q.size())<region_size){mark[u]=r*1000+tries;q.push_back(u);}}
      if(int(q.size())==region_size)S=q;
    }
    if(region_kind!="ancestral"&&int(S.size())<region_size)throw std::runtime_error("no region of the requested size");
    const int k=int(S.size());
    std::vector<std::vector<uint8_t>> adj(k,std::vector<uint8_t>(k,1));for(int i=0;i<k;++i)adj[i][i]=0;
    size_t tests=0;bool stop=false;
    for(int l=0;l<=max_level&&!stop;++l) {
      std::vector<std::vector<int>> nb(k);for(int i=0;i<k;++i)for(int j=0;j<k;++j)if(adj[i][j])nb[i].push_back(j);   // frozen
      bool any=false;
      for(int i=0;i<k&&!stop;++i)for(int j=i+1;j<k&&!stop;++j) {
        if(!adj[i][j])continue;
        bool removed=false;
        for(int side=0;side<(l==0?1:2)&&!removed&&!stop;++side) {   // level 0: one test (Z empty)
          const int a=side?j:i,b=side?i:j;
          std::vector<int> cand;for(int u:nb[a])if(u!=b)cand.push_back(u);
          if(int(cand.size())<l)continue;
          any=true;
          std::vector<int> idx(l);for(int x=0;x<l;++x)idx[x]=x;
          while(true) {
            Request q;q.x=S[i];q.y=S[j];for(int x:idx)q.z.push_back(S[cand[x]]);
            const bool sep=oracle.separated(q.x,q.y,q.z);
            out.push_back(q);++per_level[l];++tests;if(sep)++sep_level[l];
            if(sep){adj[i][j]=adj[j][i]=0;removed=true;break;}
            if(tests>=max_tests){stop=true;++capped;break;}
            int x=l-1;while(x>=0&&idx[x]==int(cand.size())-l+x)--x;   // next subset in lexicographic order
            if(x<0)break;
            ++idx[x];for(int y=x+1;y<l;++y)idx[y]=idx[y-1]+1;
          }
        }
      }
      if(!any)break;
    }
    for(int i=0;i<k;++i)for(int j=i+1;j<k;++j)edges_left+=adj[i][j];
  }
  write_workload(o.workload_out,out);
  size_t sep=0;for(size_t l=0;l<per_level.size();++l)sep+=sep_level[l];
  std::cout<<"{\"graph\":\""<<o.graph<<"\",\"region_kind\":\""<<region_kind<<"\",\"min_size\":"<<min_size<<",\"regions\":"<<regions<<",\"region_size\":"<<region_size<<",\"max_level\":"<<max_level<<",\"max_tests\":"<<max_tests
    <<",\"seed\":"<<o.seed<<",\"tests\":"<<out.size()<<",\"separated\":"<<sep<<",\"capped_regions\":"<<capped<<",\"edges_left\":"<<edges_left<<",\"per_level\":[";
  for(size_t l=0;l<per_level.size();++l)std::cout<<(l?",":"")<<per_level[l];
  std::cout<<"],\"separated_per_level\":[";for(size_t l=0;l<sep_level.size();++l)std::cout<<(l?",":"")<<sep_level[l];
  std::cout<<"],\"ms\":"<<elapsed(t)<<"}"<<std::endl;
  return 0;
}

int main(int argc,char **argv) {
  try {
    const bool causal=argc>=2&&std::string(argv[1])=="causal";
    const bool screen=argc>=2&&std::string(argv[1])=="screen";
    const bool replay=argc>=2&&std::string(argv[1])=="replay";
    const bool components=argc>=2&&std::string(argv[1])=="components";
    const bool stats=argc>=2&&std::string(argv[1])=="stats";
    const bool pcgen=argc>=2&&std::string(argv[1])=="pcgen";
    if(argc<2||(std::string(argv[1])!="bench"&&!causal&&!screen&&!replay&&!components&&!stats&&!pcgen))
      throw std::invalid_argument("usage: zskeleton bench|causal|screen|replay|components [components: --component-threshold N] [replay: --requests FILE] [screen: --hypotheses N --observed F --max-dist L --cap K] [--kinds backdoor,rule1,rule2,rule3,frontdoor,seqbackdoor,mediation,selection; causal zmodes add parents (back-door) and children (front-door)] --graph G.bin [--queries N --zsize K --zmode random,near,ancestor --zsizes 1,4,8 --pairmode random|walk --seed S --methods zskel,bidirectional,cifly --workload-in F --workload-out F --progress 0|1]");
    Options o;std::string kinds_arg="backdoor,rule1,rule2,rule3";
    size_t hypotheses=100000,cap=32;double observed_frac=0.3;int max_dist=4;std::string requests;uint32_t comp_threshold=1000;
    int pc_regions=20,pc_size=50,pc_level=1,pc_min=20;size_t pc_max_tests=1000000;std::string pc_region="ancestral";
    for(int i=2;i<argc;i+=2) {
      if(i+1>=argc)throw std::invalid_argument("option value missing");
      std::string k=argv[i],v=argv[i+1];
      if(k=="--graph")o.graph=v; else if(k=="--methods")o.methods=v; else if(k=="--zmode"||k=="--zmodes")o.zmode=v;
      else if(k=="--pairmode")o.pairmode=v; else if(k=="--workload-in")o.workload_in=v; else if(k=="--workload-out")o.workload_out=v;
      else if(k=="--queries")o.queries=std::stoull(v); else if(k=="--zsize"||k=="--zsizes")o.zsizes=v;
      else if(k=="--seed")o.seed=std::stoull(v); else if(k=="--progress")o.progress=v=="1";
      else if(k=="--kinds")kinds_arg=v;
      else if(k=="--order")o.hub_order=v;
      else if(k=="--hypotheses")hypotheses=std::stoull(v);
      else if(k=="--requests")requests=v;
      else if(k=="--regions")pc_regions=std::stoi(v);
      else if(k=="--region-size")pc_size=std::stoi(v);
      else if(k=="--max-level")pc_level=std::stoi(v);
      else if(k=="--max-tests")pc_max_tests=std::stoull(v);
      else if(k=="--region")pc_region=v;
      else if(k=="--min-size")pc_min=std::stoi(v);
      else if(k=="--component-threshold")comp_threshold=uint32_t(std::stoul(v));
      else if(k=="--observed")observed_frac=std::stod(v);
      else if(k=="--max-dist")max_dist=std::stoi(v);
      else if(k=="--cap")cap=std::stoull(v);
      else if(k=="--emit-answers")o.emit_answers=v=="1";
      else if(k=="--search"){if(v!="bi"&&v!="uni")throw std::invalid_argument("search");o.bidi=v=="bi";}
      else if(k=="--live")o.live=v=="1";
      else if(k=="--probe")o.probe=std::stoull(v);
      else if(k=="--warmup")o.warmup=std::stoi(v);
      else if(k=="--large-h")o.large_h=std::stoull(v);
      else if(k=="--initial-budget")o.initial_budget=std::stoull(v);
      else if(k=="--compress")o.compress=v=="1";
      else if(k=="--closure")o.closure_budget=std::stoull(v);
      else if(k=="--surgery-first")o.surgery_first=v=="1";
      else if(k=="--dominated")o.dominated=v=="1";
      else if(k=="--parent-live")o.parent_live=v=="1";
      else if(k=="--width")o.width=std::stoi(v);
      else if(k=="--probe-threshold")o.probe_threshold=std::stoull(v);
      else if(k=="--ancmask")o.ancmask=v=="1";
      else if(k=="--warmup-disjoint")o.warmup_disjoint=v=="1";
      else if(k=="--bp")o.bp=std::stoi(v);
      else if(k=="--topo")TrekLabels::topo_filter=v=="1";
      else if(k=="--probe-per-z")o.probe_per_z=v=="1";
      else if(k=="--balance"){if(v!="arcs"&&v!="states")throw std::invalid_argument("balance");o.arcs=v=="arcs";}
      else throw std::invalid_argument("unknown option "+k);
    }
    if(stats)return stats_main(o);
    if(pcgen)return pcgen_main(o,pc_regions,pc_size,pc_level,pc_max_tests,pc_region,pc_min);
    if(components)return components_main(o,comp_threshold);
    if(replay){if(o.methods=="bidirectional,zskel")o.methods="traversal,zskel,zskel-adaptive";
      if(requests.empty())throw std::invalid_argument("replay needs --requests FILE");
      return replay_main(o,requests);}
    if(screen){if(o.methods=="bidirectional,zskel")o.methods="traversal,zskel-adaptive";
      return screen_main(o,hypotheses,observed_frac,max_dist,cap);}
    if(causal){if(o.methods=="bidirectional,zskel")o.methods="traversal,zskel,zskel-adaptive";
      std::vector<std::string> kinds;std::stringstream ss(kinds_arg);std::string x;while(std::getline(ss,x,','))kinds.push_back(x);
      for(auto &k:kinds)if(k!="backdoor"&&k!="rule1"&&k!="rule2"&&k!="rule3"&&k!="frontdoor"&&k!="seqbackdoor"&&k!="mediation"&&k!="selection")throw std::invalid_argument("kind");
      return causal_main(o,kinds);}
    auto list=[](const std::string &v){std::vector<std::string> r;std::stringstream ss(v);std::string x;while(std::getline(ss,x,','))r.push_back(x);return r;};
    auto zmodes=list(o.zmode),methods=list(o.methods);std::vector<int> zsizes;for(auto &v:list(o.zsizes))zsizes.push_back(std::stoi(v));
    for(auto &z:zmodes)if(z!="random"&&z!="near"&&z!="ancestor")throw std::invalid_argument("zmode");
    if(o.pairmode!="random"&&o.pairmode!="walk")throw std::invalid_argument("pairmode");
    auto t=Clock::now();auto g=read_binary(o.graph);double load_ms=elapsed(t);
    std::cerr<<"{\"event\":\"graph_loaded\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"ms\":"<<load_ms<<"}\n";
    // Build every method once; workloads below never trigger a rebuild.
    struct Method {std::string name;double build_ms=0;const ZSkeleton *shared=nullptr;std::unique_ptr<ZSkeleton> idx;std::unique_ptr<ZSkeleton::Scratch> zs;
      std::unique_ptr<BidirectionalCIfly> bi;std::unique_ptr<StreamingCIfly> cf;std::unique_ptr<AdaptiveZSkeleton> ad;
      std::unique_ptr<AncestralBidirectional> an;std::unique_ptr<RPReachable> rp;std::unique_ptr<DarwichePruning> dw;std::unique_ptr<TrekSearch> ts;
      std::unique_ptr<GeigerReachable> ge;std::unique_ptr<DCollision> dc;
      std::shared_ptr<CoreTreeZSkeleton> ct;std::unique_ptr<CoreTreeZSkeleton::Scratch> cts;bool ct_built=false;};
    std::vector<Method> ms;
    for(const auto &name:methods) {
      Method m;m.name=name;auto b=Clock::now();
      if(name=="zskel"||name=="zskel-adaptive") {
        // Reuse an already built index: the adaptive engine adds no index state.
        for(auto &prev:ms)if(prev.idx){m.shared=prev.idx.get();break;}
        if(!m.shared){m.idx=std::make_unique<ZSkeleton>(g,o.progress,o.hub_order,o.compress,nullptr,o.bp);configure(*m.idx,o);m.shared=m.idx.get();}
        m.zs=std::make_unique<ZSkeleton::Scratch>(g.n);
        if(name=="zskel-adaptive"){m.ad=std::make_unique<AdaptiveZSkeleton>(*m.shared,g);m.ad->probe_budget=o.probe;m.ad->probe_per_z=o.probe_per_z;m.ad->large_h=o.large_h;m.ad->initial_budget=o.initial_budget;}
      }
      else if(name=="ctskel"||name=="ctskel-noclo") {   // core–tree index: labels on the pendant-peeled core only
        // "ctskel-noclo": the same index with up-closures off (closure budget 0 during its passes);
        // it shares the index of an earlier "ctskel" in the method list.
        for(auto &prev:ms)if(prev.ct){m.ct=prev.ct;m.build_ms=prev.build_ms;break;}
        if(!m.ct){m.ct=std::make_unique<CoreTreeZSkeleton>(g,o.progress,o.hub_order,o.compress,o.width,o.bp);configure(m.ct->inner,o);m.ct->closure_budget=o.closure_budget;m.ct->probe_threshold=o.probe_threshold;m.ct->use_anc_masks=o.ancmask;m.ct->surgery_first=o.surgery_first;m.ct->prune_dominated=o.dominated;m.ct->parents_blocked_live=o.parent_live;m.ct_built=true;}
        m.cts=std::make_unique<CoreTreeZSkeleton::Scratch>(m.ct->make_scratch());
      }
      else if(name=="bidirectional")m.bi=std::make_unique<BidirectionalCIfly>(g);
      else if(name=="bidirectional-anc")m.an=std::make_unique<AncestralBidirectional>(g);
      // Darwiche's pruning test (networkx d_separated): ancestral set, Z out-edges deleted, connectivity.
      else if(name=="darwiche")m.dw=std::make_unique<DarwichePruning>(g);
      // Classical baselines: Geiger et al.'s procedure (An(Z) marking, one-sided search over edge pairs) and the
      // d-collision graph method (Neo4j UDF variant, ported to C++; "-full" runs both searches to completion).
      else if(name=="geiger")m.ge=std::make_unique<GeigerReachable>(g);
      else if(name=="dcollision"||name=="dcollision-full")m.dc=std::make_unique<DCollision>(g,name=="dcollision-full");
      // Ablation of the index: the trek-DAG reduction with plain searches (per edge test, one sweep per query,
      // or one bidirectional search per query).
      else if(name=="trek-edge"||name=="trek-sweep"||name=="trek-bidir")
        m.ts=std::make_unique<TrekSearch>(g,name=="trek-edge"?TrekSearch::Edge:name=="trek-sweep"?TrekSearch::Sweep:TrekSearch::Bidir);
      else if(name=="cifly")m.cf=std::make_unique<StreamingCIfly>(g);
      // RP-separation (Butz et al., RPReachable with early exit). d-Separation has no directed-path
      // tests, so the "+R" baselines equal the plain ones here and are not offered.
      else if(name=="rp")m.rp=std::make_unique<RPReachable>(g,true);
      else throw std::invalid_argument("method "+name);
      if(!m.ct||m.ct_built)m.build_ms=elapsed(b);
      if(m.ct_built)std::cerr<<"{\"event\":\"fringe\",\"mib\":"<<m.ct->fringe_breakdown()<<"}\n";
      if(m.ct_built)std::cerr<<"{\"event\":\"index_built\",\"method\":\"ctskel\",\"ms\":"<<m.build_ms<<",\"entries\":"<<m.ct->inner.labels.entries()<<",\"bytes\":"<<m.ct->inner.labels.bytes()<<",\"core_n\":"<<m.ct->core.n<<",\"peeled\":"<<m.ct->peeled<<",\"peak_rss\":"<<peak_rss()<<"}\n";
      if(m.idx)std::cerr<<"{\"event\":\"index_built\",\"ms\":"<<m.build_ms<<",\"entries\":"<<m.idx->labels.entries()<<",\"bytes\":"<<m.idx->labels.bytes()<<",\"peak_rss\":"<<peak_rss()<<"}\n";
      if(!m.idx&&m.shared)for(auto &prev:ms)if(prev.idx.get()==m.shared)m.build_ms=prev.build_ms;
      ms.push_back(std::move(m));
    }
    std::vector<std::pair<std::string,int>> configs;
    // --workload-in accepts a comma list; each file is its own configuration.
    std::vector<std::string> files=list(o.workload_in);
    for(size_t f=0;f<files.size();++f)configs.push_back({"file:"+files[f],int(f)});
    if(files.empty())for(auto &z:zmodes)for(int k:zsizes)configs.push_back({z,k});
    for(auto &[zmode,zsize]:configs) {
      std::vector<Request> qs;
      std::vector<Request> wq;   // warm-up requests
      if(zmode.rfind("file:",0)==0)qs=read_workload(zmode.substr(5));
      else {Options w=o;w.zmode=zmode;w.zsize=zsize;qs=generate(g,w);
        if(o.warmup_disjoint){w.seed=o.seed+1000;wq=generate(g,w);}}
      const std::vector<Request> &warm=o.warmup_disjoint&&!wq.empty()?wq:qs;
      if(!o.workload_out.empty())write_workload(o.workload_out+"."+zmode+"-"+std::to_string(zsize),qs);
      std::cout<<std::setprecision(10)<<"{\"graph\":\""<<o.graph<<"\",\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"load_ms\":"<<load_ms
        <<",\"queries\":"<<qs.size()<<",\"zsize\":"<<zsize<<",\"zmode\":\""<<zmode<<"\",\"pairmode\":\""<<o.pairmode
        <<"\",\"seed\":"<<o.seed<<",\"search\":\""<<(o.bidi?"bi":"uni")<<"\",\"live\":"<<o.live<<",\"balance\":\""<<(o.arcs?"arcs":"states")<<"\""<<",\"warmup\":"<<o.warmup<<",\"warmup_disjoint\":"<<o.warmup_disjoint<<",\"methods\":[";
      std::vector<uint8_t> reference;std::vector<double> ref_lat;bool first=true;
      for(auto &m:ms) {
        if(m.ct)m.ct->closure_budget=m.name=="ctskel-noclo"?0:o.closure_budget;
        std::vector<uint8_t> ans(qs.size());std::vector<double> lat(qs.size());ZSkeletonStats stats;uint64_t anc_total=0;
        // Untimed warm-up passes over the same requests: every method then
        // starts its timed pass with warm caches (first-pass effects were up
        // to 10x on a fresh process).
        for(int w=0;w<o.warmup;++w)for(auto &q:warm) {
          if(m.ct)m.ct->separated(q.x,q.y,q.z,*m.cts);
          else if(m.ad)m.ad->separated(q.x,q.y,q.z,*m.zs);
          else if(m.shared)m.shared->separated(q.x,q.y,q.z,*m.zs);
          else if(m.bi)m.bi->separated(q.x,q.y,q.z);
          else if(m.an)m.an->separated(q.x,q.y,q.z);
          else if(m.dw)m.dw->separated(q.x,q.y,q.z);
          else if(m.ts)m.ts->separated(q.x,q.y,q.z);
          else if(m.ge)m.ge->separated(q.x,q.y,q.z);
          else if(m.dc)m.dc->separated(q.x,q.y,q.z);
          else if(m.rp)m.rp->separated(q.x,q.y,q.z);
          else m.cf->separated(q.x,q.y,q.z);
        }
        phase_times().reset();
        if(m.ts)m.ts->edge_tests=m.ts->states=0;
        auto q0=Clock::now();
        for(size_t i=0;i<qs.size();++i) {
          auto &q=qs[i];auto s=Clock::now();
          if(m.ct)ans[i]=m.ct->separated(q.x,q.y,q.z,*m.cts,&stats);
          else if(m.ad)ans[i]=m.ad->separated(q.x,q.y,q.z,*m.zs,&stats);
          else if(m.shared)ans[i]=m.shared->separated(q.x,q.y,q.z,*m.zs,&stats);
          else if(m.bi)ans[i]=m.bi->separated(q.x,q.y,q.z);
          else if(m.an){ans[i]=m.an->separated(q.x,q.y,q.z);anc_total+=m.an->last_ancestral;}
          else if(m.dw){ans[i]=m.dw->separated(q.x,q.y,q.z);anc_total+=m.dw->last_ancestral;}
          else if(m.ts)ans[i]=m.ts->separated(q.x,q.y,q.z);
          else if(m.ge)ans[i]=m.ge->separated(q.x,q.y,q.z);
          else if(m.dc)ans[i]=m.dc->separated(q.x,q.y,q.z);
          else if(m.rp)ans[i]=m.rp->separated(q.x,q.y,q.z);
          else ans[i]=m.cf->separated(q.x,q.y,q.z);
          lat[i]=std::chrono::duration<double,std::micro>(Clock::now()-s).count();
        }
        double query_ms=elapsed(q0);
        size_t mismatches=0,separated=0;
        if(reference.empty()){reference=ans;ref_lat=lat;}
        size_t faster=0,ties=0;for(size_t i=0;i<qs.size();++i){faster+=lat[i]<ref_lat[i];ties+=lat[i]==ref_lat[i];}
        for(size_t i=0;i<qs.size();++i){mismatches+=ans[i]!=reference[i];separated+=ans[i];}
        std::vector<double> lat_sep,lat_con;
        for(size_t i=0;i<qs.size();++i)(ans[i]?lat_sep:lat_con).push_back(lat[i]);
        double mean_us=qs.empty()?0:1000.0*query_ms/double(qs.size());
        std::cout<<(first?"":",")<<"{\"method\":\""<<m.name<<"\",\"build_ms\":"<<m.build_ms<<",\"query_ms\":"<<query_ms
          <<",\"mean_us\":"<<mean_us<<",\"p50_us\":"<<pct(lat,.5)<<",\"p90_us\":"<<pct(lat,.9)<<",\"p99_us\":"<<pct(lat,.99)<<",\"max_us\":"<<pct(lat,1)
          <<",\"sep_p50_us\":"<<pct(lat_sep,.5)<<",\"sep_p99_us\":"<<pct(lat_sep,.99)<<",\"con_p50_us\":"<<pct(lat_con,.5)<<",\"con_p99_us\":"<<pct(lat_con,.99)
          <<",\"separated\":"<<separated<<",\"mismatches_vs_first\":"<<mismatches
          <<",\"frac_faster_than_first\":"<<(qs.empty()?0:double(faster)/double(qs.size()))<<",\"frac_ties\":"<<(qs.empty()?0:double(ties)/double(qs.size()))
          <<",\"total_ms\":"<<(m.build_ms+query_ms)<<",\"anc_ms\":"<<phase_times().anc_ns/1e6;
        if(m.shared)std::cout<<",\"compressed\":"<<m.shared->labels.compressed<<",\"hub_order\":\""<<m.shared->labels.order_name<<"\",\"order_ms\":"<<m.shared->labels.order_ms<<",\"index_bytes\":"<<m.shared->labels.bytes()<<",\"label_entries\":"<<m.shared->labels.entries()<<",\"bp_hubs\":"<<m.shared->labels.bp_hubs<<",\"trek_states\":"<<m.shared->labels.states
          <<",\"pll_bfs_visits\":"<<m.shared->labels.bfs_visits<<",\"adaptive_index\":"<<stats.adaptive_index<<",\"adaptive_traversal\":"<<stats.adaptive_traversal<<",\"adaptive_rounds\":"<<stats.adaptive_rounds<<",\"h_nodes\":"<<stats.h_nodes<<",\"edge_checks\":"<<stats.edge_checks<<",\"neg_reach\":"<<stats.neg_reach
          <<",\"cert_free\":"<<stats.cert_free<<",\"searches\":"<<stats.searches<<",\"search_found\":"<<stats.search_found
          <<",\"search_states\":"<<stats.search_states<<",\"label_merges\":"<<stats.label_merges;
        if(m.an||m.dw)std::cout<<",\"mean_ancestral\":"<<(qs.empty()?0:double(anc_total)/double(qs.size()));
        if(m.ts)std::cout<<",\"trek_searches\":"<<m.ts->edge_tests<<",\"trek_states_expanded\":"<<m.ts->states;
        if(m.ct)std::cout<<",\"compressed\":"<<m.ct->inner.labels.compressed<<",\"hub_order\":\""<<m.ct->inner.labels.order_name<<"\",\"index_bytes\":"<<m.ct->inner.labels.bytes()<<",\"label_entries\":"<<m.ct->inner.labels.entries()<<",\"bp_hubs\":"<<m.ct->inner.labels.bp_hubs
          <<",\"trek_states\":"<<m.ct->inner.labels.states<<",\"core_n\":"<<m.ct->core.n<<",\"core_m\":"<<m.ct->core.m<<",\"width\":"<<m.ct->width<<",\"shortcuts\":"<<m.ct->sc.count<<",\"fringe_bytes\":"<<m.ct->fringe_bytes()<<",\"ancmask\":"<<m.ct->use_anc_masks<<",\"dominated\":"<<m.ct->prune_dominated<<",\"parent_live\":"<<m.ct->parents_blocked_live<<",\"closure_budget\":"<<m.ct->closure_budget<<",\"peeled\":"<<m.ct->peeled<<",\"max_depth\":"<<m.ct->max_depth<<",\"peel_ms\":"<<m.ct->peel_ms
          <<",\"h_nodes\":"<<stats.h_nodes<<",\"edge_checks\":"<<stats.edge_checks<<",\"neg_reach\":"<<stats.neg_reach<<",\"cert_free\":"<<stats.cert_free<<",\"searches\":"<<stats.searches<<",\"search_found\":"<<stats.search_found
          <<",\"search_states\":"<<stats.search_states<<",\"label_merges\":"<<stats.label_merges<<",\"fringe_resolved\":"<<stats.fringe_resolved<<",\"core_routes\":"<<stats.core_routes<<",\"closures\":"<<stats.closures<<",\"candidates\":"<<stats.candidates<<",\"probe_scans\":"<<stats.probe_scans<<",\"proj_fast\":"<<stats.proj_fast<<",\"proj_dfs\":"<<stats.proj_dfs<<",\"mask_builds\":"<<stats.mask_builds<<",\"mask_entries\":"<<stats.mask_entries<<",\"et_proj\":"<<stats.et_proj<<",\"et_refuted\":"<<stats.et_refuted<<",\"et_certified\":"<<stats.et_certified<<",\"et_searched\":"<<stats.et_searched<<",\"et_closure\":"<<stats.et_closure<<",\"topo_tests\":"<<m.ct->inner.labels.topo_tests<<",\"topo_pruned\":"<<m.ct->inner.labels.topo_pruned;
        if(o.emit_answers){std::cout<<",\"answers\":\"";for(auto a:ans)std::cout<<int(a);std::cout<<"\"";}
        std::cout<<"}";first=false;
        if(mismatches)std::cerr<<"{\"event\":\"MISMATCH\",\"method\":\""<<m.name<<"\",\"zmode\":\""<<zmode<<"\",\"count\":"<<mismatches<<"}\n";
        std::cerr<<"{\"event\":\"method_done\",\"method\":\""<<m.name<<"\",\"zmode\":\""<<zmode<<"\",\"zsize\":"<<zsize<<",\"query_ms\":"<<query_ms<<"}\n";
      }
      std::cout<<"],\"peak_rss_bytes\":"<<peak_rss()<<"}"<<std::endl;
    }
    return 0;
  } catch(const std::exception &e){std::cerr<<"error: "<<e.what()<<'\n';return 2;}
}
