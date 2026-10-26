#include "../src/causal.hpp"
#include <iostream>
using namespace rs;
// Reachability gadget (development note; not a result of the paper). For a DAG G0, add
// x_s -> s for every s and t -> y_t for every t. Then for s != t outside F:
//   x_s and y_t are d-connected given F in (G0 + gadget)_{\overline F}
//   iff  t is reachable from s in G0 - F.
// Checked against plain BFS on G0 - F, Bayes-ball on the physically
// surgered graph, and the Z-skeleton index (pure and adaptive).
uint64_t checks=0;
void check(bool v,const char *what){++checks;if(!v)throw std::logic_error(std::string("gadget failed: ")+what);}
int main() {
  std::mt19937_64 rng(1702);
  for(int trial=0;trial<2000;++trial) {
    int n=2+trial%12;
    std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);
    std::bernoulli_distribution e(std::uniform_real_distribution<>(0.05,0.7)(rng));
    std::vector<std::pair<int,int>> es0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(e(rng))es0.emplace_back(order[i],order[j]);
    // Gadget vertices: x_s = n+s, y_t = 2n+t.
    auto es=es0;for(int v=0;v<n;++v){es.emplace_back(n+v,v);es.emplace_back(v,2*n+v);}
    Graph g0(n,es0),g(3*n,es);
    ZSkeleton idx(g);ZSkeleton::Scratch sc(g.n);AdaptiveZSkeleton ad(idx,g);ad.initial_budget=1+trial%5;
    for(int q=0;q<20;++q) {
      int s=int(rng()%n),t=int(rng()%n);if(s==t)continue;
      std::vector<int> F;std::vector<uint8_t> inF(n);
      for(int v=0;v<n;++v)if(v!=s&&v!=t&&rng()%3==0){F.push_back(v);inF[v]=1;}
      // BFS in G0 - F.
      std::vector<uint8_t> seen(n);std::vector<int> st{s};seen[s]=1;bool reach=false;
      while(!st.empty()){int v=st.back();st.pop_back();if(v==t){reach=true;break;}for(int c:g0.ch[v])if(!seen[c]&&!inF[c]){seen[c]=1;st.push_back(c);}}
      // Physically surgered graph, Bayes-ball.
      std::vector<std::pair<int,int>> cut;for(auto [p,c]:es)if(!(c<n&&inF[c]))cut.emplace_back(p,c);
      Graph gs(3*n,cut);BayesBall bb(gs,F);
      bool sep_bb=bb.separated(n+s,2*n+t);
      bool sep_idx=idx.separated_general({n+s},{2*n+t},F,Surgery{F,{}},sc);
      bool sep_ad=ad.separated_general({n+s},{2*n+t},F,Surgery{F,{}},sc);
      check(sep_bb==!reach,"Bayes-ball on surgered gadget equals fault-tolerant reachability");
      check(sep_idx==!reach,"index equals fault-tolerant reachability");
      check(sep_ad==!reach,"adaptive engine equals fault-tolerant reachability");
    }
  }
  // (2) Z = empty: plain d-separation already encodes reachability (the
  // transitive-closure gadget of CIfly, Thm 16): x_s, y_t d-connected given
  // the empty set iff s reaches t in G0. Checked inside the loop above via
  // F = {} cases; here explicitly.
  // (3) OuMv gadget for plain d-separation with large Z: x -> s -> l_i,
  // l_i -> r_j iff M_ij, r_j -> t -> y. Z = {l_i : u_i=0} u {r_j : v_j=0}.
  // x and y are d-connected given Z iff u^T M v = 1.
  uint64_t oumv=0;
  for(int trial=0;trial<3000;++trial) {
    int a=1+int(rng()%7),b=1+int(rng()%7);
    // vertices: x=0, s=1, l_i=2+i, r_j=2+a+j, t=2+a+b, y=3+a+b
    int n=4+a+b,x=0,s0=1,t=2+a+b,y=3+a+b;
    std::vector<std::pair<int,int>> es{{x,s0},{t,y}};
    std::vector<std::vector<int>> M(a,std::vector<int>(b));
    for(int i=0;i<a;++i){es.emplace_back(s0,2+i);for(int j=0;j<b;++j)if(rng()%2){M[i][j]=1;es.emplace_back(2+i,2+a+j);}}
    for(int j=0;j<b;++j)es.emplace_back(2+a+j,t);
    Graph g(n,es);ZSkeleton idx(g);ZSkeleton::Scratch sc(n);
    for(int q=0;q<8;++q) {
      std::vector<int> u(a),v(b),Z;
      for(int i=0;i<a;++i){u[i]=int(rng()%2);if(!u[i])Z.push_back(2+i);}
      for(int j=0;j<b;++j){v[j]=int(rng()%2);if(!v[j])Z.push_back(2+a+j);}
      bool prod=false;for(int i=0;i<a;++i)for(int j=0;j<b;++j)prod|=u[i]&&M[i][j]&&v[j];
      BayesBall bb(g,Z);
      check(bb.separated(x,y)==!prod,"OuMv gadget: Bayes-ball equals u^T M v");
      check(idx.separated(x,y,Z,sc)==!prod,"OuMv gadget: index equals u^T M v");
      ++oumv;
    }
  }
  // Explicit Z = empty reachability check.
  uint64_t empty=0;
  for(int trial=0;trial<1000;++trial) {
    int n=2+trial%12;
    std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);
    std::vector<std::pair<int,int>> es0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(rng()%3==0)es0.emplace_back(order[i],order[j]);
    auto es=es0;for(int v=0;v<n;++v){es.emplace_back(n+v,v);es.emplace_back(v,2*n+v);}
    Graph g0(n,es0),g(3*n,es);BayesBall bb(g,{});
    for(int s=0;s<n;++s)for(int t=0;t<n;++t)if(s!=t){
      std::vector<uint8_t> seen(n);std::vector<int> st{s};seen[s]=1;bool r=false;
      while(!st.empty()){int v=st.back();st.pop_back();if(v==t){r=true;break;}for(int c:g0.ch[v])if(!seen[c]){seen[c]=1;st.push_back(c);}}
      check(bb.separated(n+s,2*n+t)==!r,"Z=empty: plain d-separation equals reachability");++empty;}
  }
  std::cout<<"{\"gadget_checks\":"<<checks<<",\"oumv_checks\":"<<oumv<<",\"empty_z_checks\":"<<empty<<"}\n";
}
