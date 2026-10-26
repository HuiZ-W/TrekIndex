// Local-Markov implication workloads: x ⊥ y | Pa(x) for y a non-descendant,
// non-parent of x. Every request is separated by construction, so these are
// the queries a graph-vs-data consistency check (dagitty's testable
// implications) issues. Writes one zskeleton workload file per |Pa(x)| bucket
// and per pair mode, plus a JSON line of bucket statistics per file.
//   usage: markov_workload in.bin out_prefix queries seed
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
int main(int argc,char**argv){
  if(argc!=5){fprintf(stderr,"usage: markov_workload in.bin out_prefix queries seed\n");return 2;}
  ifstream in(argv[1],ios::binary);char mg[8];uint64_t n,m;in.read(mg,8);in.read((char*)&n,8);in.read((char*)&m,8);
  if(string(mg,8)!="FDGBIN01"){fprintf(stderr,"bad magic\n");return 2;}
  vector<uint64_t> e(2*m);in.read((char*)e.data(),16*m);
  vector<vector<int>> pa(n),ch(n);
  for(uint64_t i=0;i<m;++i){ch[e[2*i]].push_back(int(e[2*i+1]));pa[e[2*i+1]].push_back(int(e[2*i]));}
  e.clear();e.shrink_to_fit();
  for(auto &p:pa){sort(p.begin(),p.end());p.erase(unique(p.begin(),p.end()),p.end());}
  string prefix=argv[2];size_t queries=stoull(argv[3]);mt19937_64 rng(stoull(argv[4]));
  vector<int> live;for(uint64_t v=0;v<n;++v)if(!pa[v].empty()||!ch[v].empty())live.push_back(int(v));
  // Buckets by |Pa(x)|; "all" samples x uniformly over non-isolated vertices.
  struct B{string name;size_t lo,hi;};
  vector<B> buckets={{"all",0,SIZE_MAX},{"pa0",0,0},{"pa1",1,1},{"pa2-3",2,3},{"pa4-7",4,7},{"pa8-15",8,15},{"pa16-63",16,63},{"pa64+",64,SIZE_MAX}};
  vector<uint32_t> mark(n,0);uint32_t ep=0;vector<int> stack;
  auto descendants=[&](int x){if(++ep==0){fill(mark.begin(),mark.end(),0);ep=1;}
    mark[x]=ep;stack.assign(1,x);
    while(!stack.empty()){int v=stack.back();stack.pop_back();for(int c:ch[v])if(mark[c]!=ep){mark[c]=ep;stack.push_back(c);}}};
  auto step=[&](int v){size_t a=pa[v].size(),b=ch[v].size();if(!a&&!b)return v;auto r=rng()%(a+b);return r<a?pa[v][r]:ch[v][r-a];};
  for(const string mode:{"walk","random"})for(auto &bk:buckets){
    vector<int> pool;for(int v:live)if(pa[v].size()>=bk.lo&&pa[v].size()<=bk.hi)pool.push_back(v);
    if(pool.size()<50){printf("{\"mode\":\"%s\",\"bucket\":\"%s\",\"candidates\":%zu,\"skipped\":1}\n",mode.c_str(),bk.name.c_str(),pool.size());continue;}
    vector<vector<int>> rows;size_t zsum=0,zmax=0,tries=0;
    while(rows.size()<queries&&tries<50*queries){++tries;
      int x=pool[rng()%pool.size()];descendants(x);
      for(int p:pa[x])mark[p]=ep;  // parents are excluded too
      int y=-1;
      for(int t=0;t<16&&y<0;++t){
        int c;
        if(mode=="walk"){c=x;int len=1+int(rng()%6);for(int i=0;i<len;++i)c=step(c);}
        else c=live[rng()%live.size()];
        if(mark[c]!=ep)y=c;
      }
      if(y<0)continue;
      vector<int> r(pa[x].begin(),pa[x].end());r.push_back(x);r.push_back(y);
      zsum+=pa[x].size();zmax=max(zmax,pa[x].size());rows.push_back(std::move(r));
    }
    string path=prefix+"."+mode+"-"+bk.name;ofstream o(path);o<<rows.size()<<'\n';
    for(auto &r:rows){size_t k=r.size()-2;o<<k;for(size_t i=0;i<k;++i)o<<' '<<r[i];o<<' '<<r[k]<<' '<<r[k+1]<<'\n';}
    printf("{\"mode\":\"%s\",\"bucket\":\"%s\",\"file\":\"%s\",\"candidates\":%zu,\"queries\":%zu,\"tries\":%zu,\"z_mean\":%.2f,\"z_max\":%zu}\n",
      mode.c_str(),bk.name.c_str(),path.c_str(),pool.size(),rows.size(),tries,rows.empty()?0.0:double(zsum)/rows.size(),zmax);
    fflush(stdout);
  }
}
