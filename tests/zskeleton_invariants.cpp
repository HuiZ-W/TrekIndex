#include "../src/zskeleton.hpp"
#include <iostream>
using namespace rs;
uint64_t checks=0;
void check(bool v,const char *what) { ++checks; if(!v) throw std::logic_error(std::string("invariant failed: ")+what); }
template<class F> void rejects(F f) { bool caught=false; try{f();}catch(const std::invalid_argument &){caught=true;} check(caught,"rejects invalid input"); }

Graph random_dag(int n,double p,std::mt19937_64 &rng) {
  std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);
  std::bernoulli_distribution edge(p);std::vector<std::pair<int,int>> es;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(edge(rng))es.emplace_back(order[i],order[j]);
  return Graph(n,es);
}
// Independent: plain BFS on the explicit trek DAG.
bool bfs_reach(const Graph &g,uint32_t s,uint32_t t) {
  std::vector<uint8_t> seen(4*size_t(g.n));std::vector<uint32_t> st{s};seen[s]=1;
  while(!st.empty()){auto v=st.back();st.pop_back();if(v==t)return true;trek_out(g,v,[&](uint32_t w){if(!seen[w]){seen[w]=1;st.push_back(w);}});}
  return false;
}
bool oracle(const Graph &g,int x,int y,const std::vector<int> &z,StreamingCIfly &cf) {
  BayesBall b(g,z);MutableRoot root(g,z);
  bool r=b.separated(x,y);
  {AncestralBidirectional an(g);check(r==an.separated(x,y,z),"bayes-ball vs ancestral bidirectional");}
  {DarwichePruning dw(g);check(r==dw.separated(x,y,z),"bayes-ball vs darwiche pruning");}
  {GeigerReachable ge(g);DCollision dc(g,false),dcf(g,true);check(r==ge.separated(x,y,z),"bayes-ball vs geiger");
   check(r==dc.separated(x,y,z),"bayes-ball vs d-collision");check(r==dcf.separated(x,y,z),"bayes-ball vs d-collision full");}
  {TrekSearch e(g,TrekSearch::Edge),w(g,TrekSearch::Sweep),b(g,TrekSearch::Bidir);check(r==e.separated(x,y,z),"bayes-ball vs trek-edge");
   check(r==w.separated(x,y,z),"bayes-ball vs trek-sweep");check(r==b.separated(x,y,z),"bayes-ball vs trek-bidir");}
  check(r==root.separated(x,y),"bayes-ball vs rootsep");
  check(r==cf.separated(x,y,z),"bayes-ball vs cifly");
  return r;
}
int main() {
  std::mt19937_64 rng(20261001);
  ZSkeletonStats total;
  // Hand cases: chain, fork, collider, collider descendant, parallel chains.
  {
    Graph g(5,{{0,1},{1,2},{3,1},{1,4}});   // 0->1->2, 3->1, 1->4
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);
    check(!idx.separated(0,2,{},s),"chain open");
    check(idx.separated(0,2,{1},s),"chain blocked");
    check(idx.separated(0,3,{},s),"collider closed");
    check(!idx.separated(0,3,{1},s),"collider opened");
    check(!idx.separated(0,3,{4},s),"collider descendant opens");
    check(!idx.separated(0,3,{4,4},s),"duplicate Z");
    check(!idx.separated(2,4,{},s),"fork open");
    check(idx.separated(2,4,{1},s),"fork blocked");
    rejects([&]{idx.separated(0,0,{},s);});
    rejects([&]{idx.separated(0,2,{0},s);});
    rejects([&]{idx.separated(0,2,{7},s);});
    rejects([&]{idx.separated_sets({0},{0},{},s);});
    rejects([&]{idx.separated_sets({0},{2},{2},s);});
  }
  {
    Graph g(4,{{0,1},{0,2},{1,3},{2,3}});   // two parallel chains 0->{1,2}->3
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);
    check(!idx.separated(0,3,{1},s),"one chain still open");
    check(idx.separated(0,3,{1,2},s),"both chains blocked");
    check(!idx.separated(1,2,{0,3},s),"collider 3 reopens after blocking fork 0");
  }
  // Labels equal explicit reachability on T.
  for(int trial=0;trial<300;++trial) {
    int n=1+trial%12;auto g=random_dag(n,std::uniform_real_distribution<>(0.05,0.8)(rng),rng);
    TrekLabels L(g),B(g,false,"butterfly");
    TrekLabels C(g,false,"degree",true);   // compressed rows
    TrekLabels P(g,false,"degree",trial%2==0,nullptr,1+trial%12);   // bit-parallel top hubs
    check(C.entries()==L.entries(),"compression keeps the canonical labels");
    for(uint32_t a=0;a<4*uint32_t(n);++a)for(uint32_t b=0;b<4*uint32_t(n);++b){
      bool r=bfs_reach(g,a,b);check(L.reach(a,b)==r,"label reach (degree)");check(B.reach(a,b)==r,"label reach (butterfly)");
      check(C.reach(a,b)==r,"label reach (compressed)");check(P.reach(a,b)==r,"label reach (bit-parallel)");}
  }
  // Exhaustive (x,y,Z) on small DAGs.
  for(int trial=0;trial<500;++trial) {
    int n=2+trial%6;auto g=random_dag(n,std::uniform_real_distribution<>(0.1,0.9)(rng),rng);
    ZSkeleton idx(g,false,trial%2?"butterfly":"degree",trial%4==1,nullptr,trial%3==0?0:trial%3==1?2:64);ZSkeleton::Scratch s(g.n);StreamingCIfly cf(g);
    for(int x=0;x<n;++x)for(int y=0;y<n;++y)if(x!=y) {
      std::vector<int> rest;for(int v=0;v<n;++v)if(v!=x&&v!=y)rest.push_back(v);
      for(unsigned mask=0;mask<(1u<<rest.size());++mask) {
        std::vector<int> z;for(size_t i=0;i<rest.size();++i)if(mask>>i&1)z.push_back(rest[i]);
        check(idx.separated(x,y,z,s,&total)==oracle(g,x,y,z,cf),"exhaustive answer");
      }
    }
  }
  // Larger random DAGs, sampled Z of every density; plus set endpoints.
  for(int trial=0;trial<1500;++trial) {
    int n=8+trial%40;auto g=random_dag(n,std::uniform_real_distribution<>(0.02,0.4)(rng),rng);
    ZSkeleton idx(g,false,"degree",false,nullptr,trial%2?0:1+trial%64);ZSkeleton::Scratch s(g.n),as(g.n);StreamingCIfly cf(g);
    idx.bidirectional_search=trial%3!=0;idx.balance_by_arcs=trial%5!=0;   // ablation variants stay exact
    AdaptiveZSkeleton adaptive(idx,g);adaptive.initial_budget=1+trial%7;adaptive.probe_budget=trial%3;   // tiny budgets force both arms
    TrekSearch te(g,TrekSearch::Edge),tw(g,TrekSearch::Sweep),tb(g,TrekSearch::Bidir);   // index ablation, scratch reused across queries
    GeigerReachable gr(g);DCollision dcr(g,false),dcfr(g,true);   // classical baselines, workspace reused
    for(int q=0;q<40;++q) {
      std::vector<int> perm(n);std::iota(perm.begin(),perm.end(),0);std::shuffle(perm.begin(),perm.end(),rng);
      int k=std::uniform_int_distribution<>(0,n-2)(rng);
      std::vector<int> z(perm.begin()+2,perm.begin()+2+k);
      bool expect=oracle(g,perm[0],perm[1],z,cf);
      check(idx.separated(perm[0],perm[1],z,s,&total)==expect,"sampled answer");
      check(adaptive.separated(perm[0],perm[1],z,as,&total)==expect,"adaptive answer");
      check(idx.separated(perm[0],perm[1],z,as)==expect,"scratch reusable after aborted rounds");
      {auto zd=z;if(!zd.empty())zd.push_back(zd[0]);   // duplicate Z entries are legal
       check(te.separated(perm[0],perm[1],zd)==expect,"trek-edge, reused");check(tw.separated(perm[0],perm[1],zd)==expect,"trek-sweep, reused");
       check(tb.separated(perm[0],perm[1],zd)==expect,"trek-bidir, reused");
       check(gr.separated(perm[0],perm[1],zd)==expect,"geiger, reused");check(dcr.separated(perm[0],perm[1],zd)==expect,"d-collision, reused");
       check(dcfr.separated(perm[0],perm[1],zd)==expect,"d-collision full, reused");}
      if(n-2-k>=2) {
        int nx=1+int(rng()%2),ny=1+int(rng()%2);if(nx+ny>n-k)ny=1;
        std::vector<int> xs(perm.begin()+2+k,perm.begin()+2+k+nx);
        std::vector<int> ys{perm[0]};if(ny>1)ys.push_back(perm[1]);else xs.push_back(perm[1]);
        if(nx+ny>n-k)continue;
        bool all=true;for(int a:xs)for(int b:ys)all=all&&oracle(g,a,b,z,cf);
        check(idx.separated_sets(xs,ys,z,s,&total)==all,"set endpoints");
      }
    }
  }
  check(total.cert_free>0 && total.searches>0 && total.search_found>0 && total.neg_reach>0,"all query paths exercised");
  check(total.adaptive_index>0 && total.adaptive_traversal>0 && total.adaptive_rounds>total.adaptive_index,"both adaptive arms and restarts exercised");
  std::cout<<"{\"checks\":"<<checks<<",\"queries\":"<<total.queries<<",\"edge_checks\":"<<total.edge_checks
    <<",\"neg_reach\":"<<total.neg_reach<<",\"cert_free\":"<<total.cert_free
    <<",\"searches\":"<<total.searches<<",\"search_found\":"<<total.search_found
    <<",\"adaptive_index\":"<<total.adaptive_index<<",\"adaptive_traversal\":"<<total.adaptive_traversal<<",\"adaptive_rounds\":"<<total.adaptive_rounds<<"}\n";
}
