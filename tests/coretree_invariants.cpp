#include "../src/coretree.hpp"
#include <iostream>
using namespace rs;
uint64_t checks=0;
void check(bool v,const char *what) { ++checks; if(!v) throw std::logic_error(std::string("invariant failed: ")+what); }
// Random DAG with a dense core and pendant trees / chains hanging off it, so
// that peeling removes a large, varied fringe. Vertex ids are shuffled.
Graph fringe_dag(int core_n,int fringe_n,double p,std::mt19937_64 &rng) {
  int n=core_n+fringe_n;std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);
  std::vector<int> topo(n);std::iota(topo.begin(),topo.end(),0);std::shuffle(topo.begin(),topo.end(),rng);   // orientation: all edges follow this order
  std::bernoulli_distribution edge(p);std::vector<std::pair<int,int>> es;
  auto add=[&](int a,int b){int u=order[a],v=order[b];if(topo[u]>topo[v])std::swap(u,v);es.emplace_back(u,v);};
  for(int i=0;i<core_n;++i)for(int j=i+1;j<core_n;++j)if(edge(rng))add(i,j);
  // fringe vertex i attaches to one earlier vertex (pendant trees), or to two (chains, cycles, parallel
  // chains between core vertices); edge directions follow the random topological order
  for(int i=std::max(core_n,1);i<n;++i){int reps=rng()%3==0?2:1;
    for(int r=0;r<reps;++r){int t=int(rng()%uint64_t(i));if(core_n>0&&rng()%3==0)t=int(rng()%uint64_t(core_n));add(t,i);}}
  return Graph(n,es);
}
// Independent oracles on the explicit surgered graph.
Graph surgered(const Graph &g,const Surgery &cut) {
  std::unordered_set<int> in(cut.cut_in.begin(),cut.cut_in.end()),out(cut.cut_out.begin(),cut.cut_out.end());
  std::vector<std::pair<int,int>> es;for(int v=0;v<g.n;++v)for(int c:g.ch[v])if(!in.count(c)&&!out.count(v))es.emplace_back(v,c);
  return Graph(g.n,es);
}
bool oracle_sets(const Graph &g,const std::vector<int> &xs,const std::vector<int> &ys,const std::vector<int> &z) {
  BayesBall b(g,z);for(int x:xs)for(int y:ys)if(!b.separated(x,y))return false;return true;
}
bool bfs_desc(const Graph &g,int u,int v,const std::unordered_set<int> &avoid) {
  std::vector<uint8_t> seen(g.n);std::vector<int> st{u};seen[u]=1;
  while(!st.empty()){int a=st.back();st.pop_back();for(int c:g.ch[a]){if(c==v)return true;if(!seen[c]&&!avoid.count(c)){seen[c]=1;st.push_back(c);}}}
  return false;
}
int main() {
  std::mt19937_64 rng(20261005);
  ZSkeletonStats total;size_t peeled_total=0,shortcuts_total=0,queries=0;
  for(int trial=0;trial<2500;++trial) {
    int core_n=trial%7==0?0:1+int(rng()%12),fringe_n=int(rng()%25);
    auto g=fringe_dag(core_n,fringe_n,std::uniform_real_distribution<>(0.1,0.7)(rng),rng);
    const int width=trial%4;   // 0 = automatic, then 1, 2, 3
    CoreTreeZSkeleton idx(g,false,"degree",false,width,trial%3==0?0:trial%3==1?3:64);ZSkeleton ref(g);
    idx.closure_budget=trial%3==0?0:trial%3==1?1:8;idx.closure_per_node=trial%3==1?0:4;idx.closure_min_h=trial%5==0?2:6;idx.closure_max_down=trial%7==0?64:4;   // off, tiny (forces over-budget fallbacks), default
    idx.probe_threshold=trial%4==0?0:trial%4==1?1:6;   // probing off, always, default
    idx.use_anc_masks=trial%2==0;
    auto s=idx.make_scratch();ZSkeleton::Scratch rs_(g.n);
    peeled_total+=idx.peeled;if(width>=2)shortcuts_total+=idx.sc.count;
    // Lemma A: core labels agree with full labels on core states, except for
    // paths that contain a detour P_c -> U_v -> D_v -> E_c through a pendant
    // parent v of a core vertex c: such a path passes U_c/P_c before and
    // E_c/D_c after the detour, so it is replaceable by U_c -> D_c unless c's
    // own states are its start or end, in which case every query forbids it.
    for(int a=0;a<idx.core.n;++a)for(int b=0;b<idx.core.n;++b)for(uint32_t ka=0;ka<4;++ka)for(uint32_t kb=0;kb<4;++kb) {
      const uint32_t sa=state(a,ka),tb=state(b,kb);
      bool c=idx.inner.labels.reach(sa,tb),f=ref.labels.reach(state(idx.core_vertex[a],ka),state(idx.core_vertex[b],kb));
      if(c==f)continue;
      bool explained=false;
      for(int cc=0;cc<idx.core.n&&!explained;++cc) {
        bool before=idx.inner.labels.reach(sa,state(cc,kUp))||idx.inner.labels.reach(sa,state(cc,kLeave));
        bool after=idx.inner.labels.reach(state(cc,kDown),tb)||idx.inner.labels.reach(state(cc,kEnter),tb);
        bool peeled_parent=false;for(int p:g.pa[idx.core_vertex[cc]])peeled_parent=peeled_parent||!idx.is_core(p);
        explained=f&&!c&&before&&after&&peeled_parent;
      }
      check(explained,"core labels differ from full labels only on detours through peeled parents");
    }
    const int n=g.n;if(n<2)continue;
    for(int q=0;q<30;++q) {
      idx.surgery_first=q%2==1;   // both exploration starts on the same graphs
      idx.prune_dominated=(q/2)%2==0;   // dominance pruning on and off
      idx.parents_blocked_live=(q/4)%2==0;   // parents-blocked endpoints on and off
      std::vector<int> perm(n);std::iota(perm.begin(),perm.end(),0);std::shuffle(perm.begin(),perm.end(),rng);
      int k=std::uniform_int_distribution<>(0,n-2)(rng);
      std::vector<int> z(perm.begin()+2,perm.begin()+2+k);
      // plain
      bool expect=oracle_sets(g,{perm[0]},{perm[1]},z);
      check(idx.separated(perm[0],perm[1],z,s,&total)==expect,"plain answer");++queries;
      check(ref.separated(perm[0],perm[1],z,rs_)==expect,"reference agrees");
      // surgery: random cut sets drawn from all vertices (endpoints and Z included)
      Surgery cut;for(int v=0;v<n;++v){if(rng()%5==0)cut.cut_in.push_back(v);if(rng()%5==0)cut.cut_out.push_back(v);}
      auto gs=surgered(g,cut);
      bool es=oracle_sets(gs,{perm[0]},{perm[1]},z);
      check(idx.separated_general({perm[0]},{perm[1]},z,cut,s,&total)==es,"surgered answer");++queries;
      // sets
      if(n-2-k>=2) {
        std::vector<int> xs{perm[0]},ys{perm[1]};int extra=n-2-k;
        if(extra>=1)xs.push_back(perm[2+k]);if(extra>=2)ys.push_back(perm[3+k]);
        check(idx.separated_sets(xs,ys,z,s,&total)==oracle_sets(g,xs,ys,z),"set endpoints");++queries;
      }
      // descendants and avoiding reachability
      int u=perm[0],v=perm[1];
      check(idx.is_descendant(u,v,s)==bfs_desc(g,u,v,{}),"descendant");
      std::unordered_set<int> av(z.begin(),z.end());std::vector<int> targets{v};if(n-2-k>=1)targets.push_back(perm[2+k]);
      bool ra=false;for(int w:targets)ra=ra||bfs_desc(g,u,w,av);
      check(idx.reaches_avoiding(u,targets,z,s,&total)==ra,"reaches avoiding");
      // Batched form: several sources outside targets and avoid.
      {std::unordered_set<int> used(targets.begin(),targets.end());used.insert(z.begin(),z.end());
       std::vector<int> us;for(int w=0;w<n&&us.size()<4;++w)if(!used.count(w))us.push_back(w);
       std::vector<uint8_t> many;idx.reaches_avoiding_many(us,targets,z,many,s,&total);
       for(size_t i=0;i<us.size();++i){bool r=false;for(int w:targets)r=r||bfs_desc(g,us[i],w,av);check(bool(many[i])==r,"reaches avoiding (batched)");}}
    }
  }
  check(shortcuts_total>0,"width 2 produced shortcut arcs");
  check(peeled_total>0&&total.fringe_resolved>0&&total.core_routes>0&&total.searches>0&&total.cert_free>0&&total.neg_reach>0&&total.closures>0&&total.candidates>0&&total.probe_scans>0,"all paths exercised");
  std::cout<<"{\"checks\":"<<checks<<",\"queries\":"<<queries<<",\"peeled\":"<<peeled_total<<",\"shortcuts\":"<<shortcuts_total<<",\"fringe_resolved\":"<<total.fringe_resolved
    <<",\"closures\":"<<total.closures<<",\"core_routes\":"<<total.core_routes<<",\"searches\":"<<total.searches<<",\"cert_free\":"<<total.cert_free<<"}\n";
}
