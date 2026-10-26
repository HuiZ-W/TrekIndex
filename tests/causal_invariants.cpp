#include "../src/causal.hpp"
#include "../src/coretree.hpp"
#include "../src/backdoor_baselines.hpp"
#include "../src/surgery_baselines.hpp"
#include <iostream>
using namespace rs;
// Independent oracles: physically build every surgered graph, then run
// Bayes-ball and RootSep on it; ancestors/descendants by plain BFS.
uint64_t checks=0;
void check(bool v,const char *what){++checks;if(!v)throw std::logic_error(std::string("invariant failed: ")+what);}
template<class F> void rejects(F f){bool caught=false;try{f();}catch(const std::invalid_argument &){caught=true;}check(caught,"rejects invalid input");}

std::vector<std::pair<int,int>> edges_of(const Graph &g){std::vector<std::pair<int,int>> e;for(int v=0;v<g.n;++v)for(int c:g.ch[v])e.emplace_back(v,c);return e;}
Graph surgered(const Graph &g,const std::vector<int> &cut_in,const std::vector<int> &cut_out) {
  std::vector<uint8_t> in(g.n),out(g.n);for(int v:cut_in)in[v]=1;for(int v:cut_out)out[v]=1;
  std::vector<std::pair<int,int>> e;for(auto [p,c]:edges_of(g))if(!in[c]&&!out[p])e.emplace_back(p,c);
  return Graph(g.n,e);
}
bool oracle_sep(const Graph &g,const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &zs) {
  bool all=true;
  for(int x:xs)for(int y:ys){
    BayesBall b(g,zs);MutableRoot r(g,zs);bool s=b.separated(x,y);
    check(s==r.separated(x,y),"oracle bayes-ball vs rootsep");all=all&&s;
  }
  return all;
}
std::vector<uint8_t> reach_down(const Graph &g,const std::vector<int> &from){   // proper descendants
  std::vector<uint8_t> seen(g.n);std::vector<int> st;
  for(int v:from)for(int c:g.ch[v])if(!seen[c]){seen[c]=1;st.push_back(c);}
  while(!st.empty()){int v=st.back();st.pop_back();for(int c:g.ch[v])if(!seen[c]){seen[c]=1;st.push_back(c);}}
  return seen;
}
std::vector<uint8_t> reach_up(const Graph &g,const std::vector<int> &from){     // proper ancestors
  std::vector<uint8_t> seen(g.n);std::vector<int> st;
  for(int v:from)for(int p:g.pa[v])if(!seen[p]){seen[p]=1;st.push_back(p);}
  while(!st.empty()){int v=st.back();st.pop_back();for(int p:g.pa[v])if(!seen[p]){seen[p]=1;st.push_back(p);}}
  return seen;
}
Graph random_dag(int n,double p,std::mt19937_64 &rng) {
  std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);
  std::bernoulli_distribution edge(p);std::vector<std::pair<int,int>> es;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(edge(rng))es.emplace_back(order[i],order[j]);
  return Graph(n,es);
}
std::vector<int> cat(std::vector<int> a,const std::vector<int> &b){a.insert(a.end(),b.begin(),b.end());return a;}

int main() {
  ZSkeletonStats stats;
  // Textbook cases. Vertex names in comments.
  {
    Graph g(3,{{0,1},{1,2},{0,2}});   // Z=0 -> X=1 -> Y=2, Z -> Y
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);CausalQueries<const ZSkeleton> q(idx,idx);
    check(q.backdoor({1},{2},{0},s).valid,"confounder adjusted");
    auto open=q.backdoor({1},{2},{},s);check(!open.valid&&open.open_backdoor_path,"confounding open");
    check(!q.rule2({2},{1},{},{},s),"observational != interventional under confounding");
    check(q.rule2({2},{1},{},{0},s),"exchange after adjusting the confounder");
  }
  {
    Graph g(5,{{1,0},{1,2},{3,2},{3,4},{0,4}});   // M-bias: X=0 <- A=1 -> M=2 <- B=3 -> Y=4, X -> Y
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);CausalQueries<const ZSkeleton> q(idx,idx);
    check(q.backdoor({0},{4},{},s).valid,"M-bias: empty set valid");
    check(!q.backdoor({0},{4},{2},s).valid,"M-bias: conditioning on collider M opens");
    check(q.backdoor({0},{4},{2,1},s).valid,"M-bias: M plus A blocks again");
  }
  {
    Graph g(3,{{0,1},{1,2}});   // X=0 -> M=1 -> Y=2
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);CausalQueries<const ZSkeleton> q(idx,idx);
    auto r=q.backdoor({0},{2},{1},s);check(!r.valid&&r.descendant_in_z&&r.x==0&&r.z==1,"mediator is a descendant");
    check(q.rule3({0},{2},{},{},s),"do(Y) has no effect on upstream X");
    check(!q.rule3({2},{0},{},{},s),"do(X) affects downstream Y");
    check(q.rule1({2},{0},{1},{},s),"rule 1: Y independent of X given do(M)");
    rejects([&]{q.rule1({2},{2},{},{},s);});
    rejects([&]{q.backdoor({0},{2},{0},s);});
  }
  {
    Graph g(4,{{0,1},{0,3},{1,2},{2,3}});   // front-door: U=0 -> X=1 -> M=2 -> Y=3, U -> Y
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);CausalQueries<const ZSkeleton> q(idx,idx);
    check(q.frontdoor({1},{3},{2},s).valid,"front-door: textbook case");
    auto r=q.frontdoor({1},{3},{},s);check(!r.valid&&r.unintercepted,"front-door: empty M does not intercept");
    check(!q.backdoor({1},{3},{},s).valid,"front-door case: back-door fails without U");
    check(q.mediation({1},{2},{3},{0},s),"mediation identified given U");
    check(!q.mediation({1},{2},{3},{},s),"mediation not identified with U unobserved");
    check(q.selection_recoverable({1},{3},{2},s)==false,"selection on a mediator");
    check(q.selection_recoverable({1,2,0},{3},{2},s),"selection on a conditioned variable");
  }
  {
    Graph g(5,{{0,1},{0,4},{1,2},{2,4},{3,2},{3,4}});   // U=0 -> X=1 -> M=2 -> Y=4, U -> Y, V=3 -> M, V -> Y
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n);CausalQueries<const ZSkeleton> q(idx,idx);
    auto r=q.frontdoor({1},{4},{2},s);check(!r.valid&&r.backdoor_m_y,"front-door: M-Y confounding through V");
  }
  // Randomized surgery and criteria against explicit graphs.
  std::mt19937_64 rng(4242);
  for(int trial=0;trial<3000;++trial) {
    int n=3+trial%14;auto g=random_dag(n,std::uniform_real_distribution<>(0.05,0.6)(rng),rng);
    ZSkeleton idx(g,false,"degree",false,nullptr,trial%3==0?0:trial%3==1?4:64);ZSkeleton::Scratch s(g.n),as(g.n);
    idx.live_states=trial%2==0;idx.bidirectional_search=trial%3!=0;   // ablation variants stay exact
    AdaptiveZSkeleton adaptive(idx,g);adaptive.initial_budget=1+trial%9;adaptive.probe_budget=trial%4?trial%17:0;
    CausalQueries<const ZSkeleton> pure(idx,idx);CausalQueries<AdaptiveZSkeleton> ad(idx,adaptive);
    // Core–tree index through the same criteria (random DAGs have few pendants; tests/coretree_invariants.cpp covers deep fringes).
    CoreTreeZSkeleton ct(g,false,"degree",false,1+trial%2,trial%5==0?0:16);ct.surgery_first=(trial/2)%2==1;ct.parents_blocked_live=(trial/4)%2==0;auto cs=ct.make_scratch();CausalQueries<const CoreTreeZSkeleton,CoreTreeZSkeleton> qc(ct,ct);
    BaselineCausal<BidirectionalGeneral> b1(g);BaselineCausal<AncestralGeneral> b2(g);BaselineCausal<CIflyGeneral> b3(g);BaselineCausal<RPGeneral> b4(g);
    BaselineCausal<GeigerGeneral> b5(g);BaselineCausal<DCollisionGeneral> b6(g);   // classical baselines
    for(int q=0;q<12;++q) {
      // Partition: X, Y, Z, W, rest.
      std::vector<int> role(n);for(auto &r:role)r=int(rng()%5);
      std::vector<int> X,Y,Z,W;for(int v=0;v<n;++v){if(role[v]==0)X.push_back(v);if(role[v]==1)Y.push_back(v);if(role[v]==2)Z.push_back(v);if(role[v]==3)W.push_back(v);}
      // Arbitrary surgery, overlapping any role.
      std::vector<int> ci,co;for(int v=0;v<n;++v){if(rng()%4==0)ci.push_back(v);if(rng()%4==0)co.push_back(v);}
      if(!X.empty()&&!Y.empty()) {
        bool expect=oracle_sep(surgered(g,ci,co),X,Y,Z);
        check(idx.separated_general(X,Y,Z,Surgery{ci,co},s,&stats)==expect,"surgered separation (index)");
        check(adaptive.separated_general(X,Y,Z,Surgery{ci,co},as,&stats)==expect,"surgered separation (adaptive)");
        check(b1.sep(X,Y,Z,ci,co)==expect&&b2.sep(X,Y,Z,ci,co)==expect&&b3.sep(X,Y,Z,ci,co)==expect&&b4.sep(X,Y,Z,ci,co)==expect&&b5.sep(X,Y,Z,ci,co)==expect&&b6.sep(X,Y,Z,ci,co)==expect,"surgered separation (six baselines)");
        // Back-door criterion.
        auto de=reach_down(g,X);bool desc=false;for(int z:Z)desc|=bool(de[z]);
        bool bd=!desc&&oracle_sep(surgered(g,{},X),X,Y,Z);
        auto r1=pure.backdoor(X,Y,Z,s,&stats),r2=ad.backdoor(X,Y,Z,as,&stats),r3=qc.backdoor(X,Y,Z,cs,&stats);
        check(r1.valid==bd&&r2.valid==bd&&r3.valid==bd,"back-door criterion");
        check(b1.backdoor(X,Y,Z)==bd&&b2.backdoor(X,Y,Z)==bd&&b3.backdoor(X,Y,Z)==bd&&b4.backdoor(X,Y,Z)==bd&&b5.backdoor(X,Y,Z)==bd&&b6.backdoor(X,Y,Z)==bd,"back-door criterion (six baselines)");
        check(r3.descendant_in_z==desc,"back-door descendant flag (core–tree)");
        check(r1.descendant_in_z==desc,"back-door descendant flag");
        if(desc)check(de[r1.z]&&std::find(X.begin(),X.end(),r1.x)!=X.end(),"back-door witness");
      }
      if(!Y.empty()&&!Z.empty()) {
        auto cond=cat(X,W);
        bool e1=oracle_sep(surgered(g,X,{}),Y,Z,cond);
        bool e2=oracle_sep(surgered(g,X,Z),Y,Z,cond);
        auto anc=reach_up(surgered(g,X,{}),W);std::vector<int> cut=X;for(int z:Z)if(!anc[z])cut.push_back(z);
        bool e3=oracle_sep(surgered(g,cut,{}),Y,Z,cond);
        check(pure.rule1(Y,Z,X,W,s,&stats)==e1&&ad.rule1(Y,Z,X,W,as,&stats)==e1,"do-calculus rule 1");
        check(pure.rule2(Y,Z,X,W,s,&stats)==e2&&ad.rule2(Y,Z,X,W,as,&stats)==e2,"do-calculus rule 2");
        check(pure.rule3(Y,Z,X,W,s,&stats)==e3&&ad.rule3(Y,Z,X,W,as,&stats)==e3,"do-calculus rule 3");
        check(qc.rule1(Y,Z,X,W,cs,&stats)==e1&&qc.rule2(Y,Z,X,W,cs,&stats)==e2&&qc.rule3(Y,Z,X,W,cs,&stats)==e3,"do-calculus rules (core–tree)");
        for(int w=1;w<=3;++w){bool e=w==1?e1:w==2?e2:e3;check(b1.rule(w,Y,Z,X,W)==e&&b2.rule(w,Y,Z,X,W)==e&&b3.rule(w,Y,Z,X,W)==e&&b4.rule(w,Y,Z,X,W)==e&&b5.rule(w,Y,Z,X,W)==e&&b6.rule(w,Y,Z,X,W)==e,"do-calculus rules (six baselines)");}
      }
      // Front-door (M = Z), mediation (M = Z), sequential back-door, selection.
      if(!X.empty()&&!Y.empty()) {
        auto through=reach_down(surgered(g,Z,{}),X);bool open_path=false;for(int y:Y)open_path|=bool(through[y]);
        bool fd=!open_path&&(Z.empty()||(oracle_sep(surgered(g,{},X),X,Z,{})&&oracle_sep(surgered(g,{},Z),Z,Y,X)));
        check(pure.frontdoor(X,Y,Z,s,&stats).valid==fd&&ad.frontdoor(X,Y,Z,as,&stats).valid==fd&&qc.frontdoor(X,Y,Z,cs,&stats).valid==fd,"front-door criterion");
        check(b1.frontdoor(X,Y,Z)==fd&&b2.frontdoor(X,Y,Z)==fd&&b3.frontdoor(X,Y,Z)==fd&&b4.frontdoor(X,Y,Z)==fd&&b5.frontdoor(X,Y,Z)==fd&&b6.frontdoor(X,Y,Z)==fd,"front-door criterion (six baselines)");
        if(!Z.empty()) {
          auto dx=reach_down(g,cat(X,Z));bool desc=false;for(int w:W)desc|=bool(dx[w]);
          bool md=!desc&&oracle_sep(surgered(g,{},Z),Z,Y,cat(X,W))&&oracle_sep(surgered(g,{},X),X,Z,W)&&oracle_sep(surgered(g,Z,X),X,Y,W);
          check(pure.mediation(X,Z,Y,W,s,&stats)==md&&ad.mediation(X,Z,Y,W,as,&stats)==md&&qc.mediation(X,Z,Y,W,cs,&stats)==md,"mediation criterion");
          check(b1.mediation(X,Z,Y,W)==md&&b2.mediation(X,Z,Y,W)==md&&b3.mediation(X,Z,Y,W)==md&&b4.mediation(X,Z,Y,W)==md&&b5.mediation(X,Z,Y,W)==md&&b6.mediation(X,Z,Y,W)==md,"mediation criterion (six baselines)");
        }
        std::vector<std::vector<int>> ws(X.size());for(size_t i=0;i<W.size();++i)ws[i%X.size()].push_back(W[i]);
        bool sb=true;std::vector<int> cond;
        for(size_t k=0;k<X.size()&&sb;++k) {
          auto dk=reach_down(g,std::vector<int>(X.begin()+long(k),X.end()));for(int w:ws[k])if(dk[w])sb=false;
          cond=cat(cond,ws[k]);
          if(sb)sb=oracle_sep(surgered(g,std::vector<int>(X.begin()+long(k)+1,X.end()),{X[k]}),Y,{X[k]},cond);
          cond.push_back(X[k]);
        }
        check(pure.sequential_backdoor(X,ws,Y,s,&stats)==sb&&ad.sequential_backdoor(X,ws,Y,as,&stats)==sb&&qc.sequential_backdoor(X,ws,Y,cs,&stats)==sb,"sequential back-door");
        check(b1.sequential_backdoor(X,ws,Y)==sb&&b2.sequential_backdoor(X,ws,Y)==sb&&b3.sequential_backdoor(X,ws,Y)==sb&&b4.sequential_backdoor(X,ws,Y)==sb&&b5.sequential_backdoor(X,ws,Y)==sb&&b6.sequential_backdoor(X,ws,Y)==sb,"sequential back-door (six baselines)");
      }
      if(!Y.empty()) {
        std::vector<int> ps;for(int v=0;v<n;++v)if(rng()%4==0)ps.push_back(v);
        auto es=edges_of(g);for(int p:ps)es.emplace_back(p,n);Graph gs(n+1,es);
        bool sr=oracle_sep(gs,Y,{n},X);
        check(pure.selection_recoverable(X,Y,ps,s,&stats)==sr&&ad.selection_recoverable(X,Y,ps,as,&stats)==sr&&qc.selection_recoverable(X,Y,ps,cs,&stats)==sr,"selection recoverability");
        check(b1.selection_recoverable(X,Y,ps)==sr&&b2.selection_recoverable(X,Y,ps)==sr&&b3.selection_recoverable(X,Y,ps)==sr&&b4.selection_recoverable(X,Y,ps)==sr&&b5.selection_recoverable(X,Y,ps)==sr&&b6.selection_recoverable(X,Y,ps)==sr,"selection recoverability (six baselines)");
      }
      // Single-pair baselines: plain separation and their back-door mode (out-edges of x removed).
      {
        BackdoorBaseline<BidirectionalCIfly> bb(g);BackdoorBaseline<AncestralBidirectional> ba(g);
        BackdoorBaseline<StreamingCIfly> bc(g);BackdoorBaseline<RPReachable> br(g,true);
        BackdoorBaseline<GeigerReachable> bg(g);BackdoorBaseline<DCollision> bdc(g,false);
        int x=int(rng()%n),y=int(rng()%n);
        if(x!=y) {
          std::vector<int> Zp;for(int w=0;w<n;++w)if(w!=x&&w!=y&&rng()%3==0)Zp.push_back(w);
          bool plain=oracle_sep(g,{x},{y},Zp);
          check(bb.sep.separated(x,y,Zp)==plain&&ba.sep.separated(x,y,Zp)==plain&&bc.sep.separated(x,y,Zp)==plain&&br.sep.separated(x,y,Zp)==plain&&bg.sep.separated(x,y,Zp)==plain&&bdc.sep.separated(x,y,Zp)==plain,"baselines: plain separation");
          bool cut=oracle_sep(surgered(g,{},{x}),{x},{y},Zp);
          check(bb.sep.separated_cut_out(x,y,Zp)==cut&&ba.sep.separated_cut_out(x,y,Zp)==cut&&bc.sep.separated_cut_out(x,y,Zp)==cut&&br.sep.separated_cut_out(x,y,Zp)==cut&&bg.sep.separated_cut_out(x,y,Zp)==cut&&bdc.sep.separated_cut_out(x,y,Zp)==cut,"baselines: separation in the back-door graph");
          auto de=reach_down(g,{x});bool desc=false;for(int z:Zp)desc|=bool(de[z]);
          bool bd=!desc&&cut;
          check(bb.valid({x},{y},Zp)==bd&&ba.valid({x},{y},Zp)==bd&&bc.valid({x},{y},Zp)==bd&&br.valid({x},{y},Zp)==bd&&bg.valid({x},{y},Zp)==bd&&bdc.valid({x},{y},Zp)==bd,"baselines: back-door criterion");
        }
      }
      // Ancestry primitives.
      int u=int(rng()%n),v=int(rng()%n);
      check(idx.is_descendant(u,v)==bool(reach_down(g,{u})[v]),"descendant");
      if(u!=v) {
        std::vector<int> avoid;for(int w=0;w<n;++w)if(w!=u&&w!=v&&rng()%3==0)avoid.push_back(w);
        check(idx.reaches_avoiding(u,{v},avoid,s,&stats)==bool(reach_up(surgered(g,avoid,{}),{v})[u]),"ancestor avoiding");
        std::vector<int> us;for(int w=0;w<n&&us.size()<4;++w)if(w!=v&&std::find(avoid.begin(),avoid.end(),w)==avoid.end())us.push_back(w);
        std::vector<uint8_t> many;idx.reaches_avoiding_many(us,{v},avoid,many,s,&stats);auto up=reach_up(surgered(g,avoid,{}),{v});
        for(size_t i=0;i<us.size();++i)check(bool(many[i])==bool(up[us[i]]),"ancestor avoiding (batched)");
      }
    }
  }
  // Large-Z back-door on larger DAGs: many label-only descendant checks and
  // the h(h-1)/2 adaptive probe.
  for(int trial=0;trial<60;++trial) {
    int n=150+trial*5;auto g=random_dag(n,std::uniform_real_distribution<>(0.03,0.08)(rng),rng);
    ZSkeleton idx(g);ZSkeleton::Scratch s(g.n),as(g.n);AdaptiveZSkeleton adaptive(idx,g);
    CausalQueries<const ZSkeleton> pure(idx,idx);CausalQueries<AdaptiveZSkeleton> ad(idx,adaptive);
    adaptive.large_h=trial%3?8:64;
    for(int q=0;q<20;++q) {
      std::vector<int> perm(n);std::iota(perm.begin(),perm.end(),0);std::shuffle(perm.begin(),perm.end(),rng);
      std::vector<int> X{perm[0]},Y{perm[1]},Z(perm.begin()+2,perm.begin()+2+int(rng()%(n/2)));
      auto de=reach_down(g,X);bool desc=false;for(int z:Z)desc|=bool(de[z]);
      bool bd=!desc&&oracle_sep(surgered(g,{},X),X,Y,Z);
      auto r1=pure.backdoor(X,Y,Z,s,&stats),r2=ad.backdoor(X,Y,Z,as,&stats);
      check(r1.valid==bd&&r2.valid==bd,"large-Z back-door");
      check(r1.descendant_in_z==desc&&r2.descendant_in_z==desc,"large-Z descendant flag");
      if(desc)check(de[r1.z]&&de[r2.z],"large-Z descendant witness");
    }
  }
  check(stats.searches>0&&stats.cert_free>0&&stats.adaptive_index>0&&stats.adaptive_traversal>0,"all paths exercised");
  std::cout<<"{\"checks\":"<<checks<<",\"queries\":"<<stats.queries<<",\"searches\":"<<stats.searches
    <<",\"adaptive_index\":"<<stats.adaptive_index<<",\"adaptive_traversal\":"<<stats.adaptive_traversal<<"}\n";
}
