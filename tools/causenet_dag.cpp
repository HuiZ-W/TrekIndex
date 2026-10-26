// Build a DAG from CauseNet relations by evidence-ordered cycle breaking.
// Input:  TSV cause, effect, support, source-types (scripts/extract_causenet.py)
// Output: <out>.bin (FDGBIN01), <out>.names.tsv (id, concept, in/out degree),
//         <out>.stats.json
// Every edge inside a strongly connected component lies on a cycle. For each
// nontrivial SCC, remove its internal edges of the minimum support level, then
// recompute SCCs inside it, until no cycle remains. Edges are never reversed.
// Strong evidence survives; removed edges and their support levels are
// reported. Deterministic.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
struct Edge{uint32_t u,v,w;bool alive=true;};
int main(int argc,char**argv){
  if(argc!=3){std::cerr<<"usage: causenet_dag relations.tsv out_prefix\n";return 2;}
  std::ifstream in(argv[1]);std::string line;
  std::unordered_map<std::string,uint32_t> id;std::vector<std::string> names;std::vector<Edge> es;
  auto nid=[&](const std::string &s){auto it=id.find(s);if(it!=id.end())return it->second;uint32_t k=uint32_t(names.size());id.emplace(s,k);names.push_back(s);return k;};
  size_t lines=0;
  while(std::getline(in,line)){
    ++lines;size_t a=line.find('\t'),b=line.find('\t',a+1),c=line.find('\t',b+1);
    if(a==std::string::npos||b==std::string::npos)continue;
    uint32_t u=nid(line.substr(0,a)),v=nid(line.substr(a+1,b-a-1));
    uint32_t w=uint32_t(std::stoul(line.substr(b+1,(c==std::string::npos?line.size():c)-b-1)));
    if(u!=v)es.push_back({u,v,w});
  }
  // Merge duplicate (u,v) pairs by summing support.
  std::sort(es.begin(),es.end(),[](const Edge&x,const Edge&y){return x.u!=y.u?x.u<y.u:x.v<y.v;});
  {std::vector<Edge> m;for(auto&e:es){if(!m.empty()&&m.back().u==e.u&&m.back().v==e.v)m.back().w+=e.w;else m.push_back(e);}es.swap(m);}
  const uint32_t n=uint32_t(names.size());
  std::cerr<<"{\"event\":\"loaded\",\"lines\":"<<lines<<",\"nodes\":"<<n<<",\"edges\":"<<es.size()<<"}\n";
  // CSR over edge indices (out-edges), rebuilt implicitly via alive flags.
  std::vector<uint64_t> off(n+1);for(auto&e:es)off[e.u+1]++;for(uint32_t i=0;i<n;++i)off[i+1]+=off[i];
  std::vector<uint32_t> adj(es.size());{auto p=off;for(uint32_t k=0;k<es.size();++k)adj[p[es[k].u]++]=k;}
  // Iterative Tarjan restricted to a vertex set (members marked with tag).
  std::vector<uint32_t> idx(n,UINT32_MAX),low(n),comp(n,UINT32_MAX),tag(n,0);std::vector<uint8_t> on(n,0);
  uint32_t tagc=0;
  std::map<uint32_t,uint64_t> removed_by_support;uint64_t removed=0,rounds=0,removed_support_total=0;
  auto sccs=[&](const std::vector<uint32_t>&verts,std::vector<std::vector<uint32_t>>&out){
    ++tagc;for(auto v:verts){tag[v]=tagc;idx[v]=UINT32_MAX;}
    uint32_t counter=0;std::vector<uint32_t> st;std::vector<std::pair<uint32_t,uint64_t>> work;
    for(auto r:verts){if(idx[r]!=UINT32_MAX)continue;work.push_back({r,off[r]});idx[r]=low[r]=counter++;st.push_back(r);on[r]=1;
      while(!work.empty()){auto &[v,i]=work.back();
        if(i<off[v+1]){const Edge&e=es[adj[i++]];if(!e.alive||tag[e.v]!=tagc)continue;uint32_t w=e.v;
          if(idx[w]==UINT32_MAX){idx[w]=low[w]=counter++;st.push_back(w);on[w]=1;work.push_back({w,off[w]});}
          else if(on[w])low[v]=std::min(low[v],idx[w]);}
        else{uint32_t vv=v;work.pop_back();if(!work.empty())low[work.back().first]=std::min(low[work.back().first],low[vv]);
          if(low[vv]==idx[vv]){std::vector<uint32_t> c;uint32_t w;do{w=st.back();st.pop_back();on[w]=0;c.push_back(w);}while(w!=vv);
            if(c.size()>1)out.push_back(std::move(c));}}}}
  };
  std::vector<uint32_t> all(n);for(uint32_t i=0;i<n;++i)all[i]=i;
  std::vector<std::vector<uint32_t>> work;sccs(all,work);
  size_t initial_nontrivial=work.size(),largest=0;for(auto&c:work)largest=std::max(largest,c.size());
  std::cerr<<"{\"event\":\"initial_sccs\",\"nontrivial\":"<<initial_nontrivial<<",\"largest\":"<<largest<<"}\n";
  while(!work.empty()){
    auto c=std::move(work.back());work.pop_back();++rounds;
    ++tagc;for(auto v:c)tag[v]=tagc;const uint32_t ct=tagc;
    uint32_t wmin=UINT32_MAX;
    for(auto v:c)for(uint64_t i=off[v];i<off[v+1];++i){const Edge&e=es[adj[i]];if(e.alive&&tag[e.v]==ct)wmin=std::min(wmin,e.w);}
    for(auto v:c)for(uint64_t i=off[v];i<off[v+1];++i){Edge&e=es[adj[i]];
      if(e.alive&&tag[e.v]==ct&&e.w==wmin){e.alive=false;++removed;removed_support_total+=e.w;++removed_by_support[e.w];}}
    std::vector<std::vector<uint32_t>> sub;sccs(c,sub);for(auto&s:sub)work.push_back(std::move(s));
  }
  // Verify acyclicity (Kahn) and write outputs.
  std::vector<uint32_t> indeg(n,0),outdeg(n,0);uint64_t kept=0,kept_support=0;
  for(auto&e:es)if(e.alive){indeg[e.v]++;outdeg[e.u]++;++kept;kept_support+=e.w;}
  {std::vector<uint32_t> d=indeg,q;for(uint32_t v=0;v<n;++v)if(!d[v])q.push_back(v);size_t seen=0;
   while(!q.empty()){uint32_t v=q.back();q.pop_back();++seen;for(uint64_t i=off[v];i<off[v+1];++i){const Edge&e=es[adj[i]];if(e.alive&&--d[e.v]==0)q.push_back(e.v);}}
   if(seen!=n){std::cerr<<"not acyclic\n";return 3;}}
  std::string pre=argv[2];
  {std::ofstream b(pre+".bin",std::ios::binary);b.write("FDGBIN01",8);uint64_t nn=n,mm=kept;b.write((char*)&nn,8);b.write((char*)&mm,8);
   for(auto&e:es)if(e.alive){uint64_t u=e.u,v=e.v;b.write((char*)&u,8);b.write((char*)&v,8);}}
  {std::ofstream t(pre+".names.tsv");for(uint32_t v=0;v<n;++v)t<<v<<'\t'<<names[v]<<'\t'<<indeg[v]<<'\t'<<outdeg[v]<<'\n';}
  {std::ofstream w(pre+".support.tsv");for(auto&e:es)if(e.alive)w<<e.u<<'\t'<<e.v<<'\t'<<e.w<<'\n';}
  std::ofstream s(pre+".stats.json");
  s<<"{\"method\":\"evidence-ordered-scc-peeling-v1\",\"nodes\":"<<n<<",\"input_edges\":"<<es.size()<<",\"kept_edges\":"<<kept
   <<",\"removed_edges\":"<<removed<<",\"kept_support\":"<<kept_support<<",\"removed_support\":"<<removed_support_total
   <<",\"initial_nontrivial_sccs\":"<<initial_nontrivial<<",\"initial_largest_scc\":"<<largest<<",\"peel_rounds\":"<<rounds
   <<",\"removed_by_support\":{";
  bool f=true;for(auto&[w,c]:removed_by_support){s<<(f?"":",")<<"\""<<w<<"\":"<<c;f=false;}
  s<<"},\"edge_direction_changed\":false}\n";
  std::cerr<<"{\"event\":\"done\",\"kept\":"<<kept<<",\"removed\":"<<removed<<",\"rounds\":"<<rounds<<"}\n";
}
