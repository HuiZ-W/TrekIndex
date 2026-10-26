// Standalone driver of RP-separation (Butz et al., PGM 2016; rs::RPReachable in paper_algorithms.hpp)
// for saved d-separation workloads, as used for the RP column of Exp-1 (scripts/replay_paper_baselines.py).
// It was part of a larger driver of earlier baselines that are not in the paper; only RP is kept.
//   usage: paper_baselines bench rp-boolean SEED graph.bin < workload   (rp-full: no early exit)
//          paper_baselines probe graph.bin
#include "paper_algorithms.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/resource.h>
struct Request {std::vector<int> z;int x,y;};
rs::Graph input_graph(const std::string &path=""){
 std::cerr<<"phase: graph_load_start\n";
 uint64_t n=0,m=0;std::ifstream file;
 auto u64=[&](){unsigned char b[8];if(!file.read(reinterpret_cast<char*>(b),8))throw std::invalid_argument("truncated binary");
  uint64_t v=0;for(int i=0;i<8;++i)v|=uint64_t(b[i])<<(8*i);return v;};
 if(path.empty()){if(!(std::cin>>n>>m))throw std::invalid_argument("graph");}
 else {file.open(path,std::ios::binary);char magic[8];if(!file.read(magic,8)||std::string(magic,8)!="FDGBIN01")throw std::invalid_argument("FDGBIN01 magic");n=u64();m=u64();
  file.seekg(0,std::ios::end);auto length=file.tellg();if(m>(UINT64_MAX-24)/16||length<0||uint64_t(length)!=24+16*m)throw std::invalid_argument("binary length");file.seekg(24);}
 // Validate capacity before any graph/edge allocation, including malformed inputs.
 if(n<2||n>uint64_t(std::numeric_limits<int>::max()/2)||m>2147483647)throw std::length_error("native graph ID/size limit");
 std::vector<std::pair<int,int>> es;es.reserve(size_t(m));
 for(uint64_t i=0;i<m;++i){uint64_t u,v;if(path.empty()){if(!(std::cin>>u>>v))throw std::invalid_argument("edge");}else{u=u64();v=u64();}
  if(u>=n||v>=n)throw std::invalid_argument("edge endpoint");es.emplace_back(int(u),int(v));}
 auto graph=rs::Graph(int(n),es);std::cerr<<"phase: graph_load_complete n="<<graph.n<<" m="<<graph.m<<'\n';return graph;
}
std::vector<Request> input_queries(){int n;if(!(std::cin>>n)||n<1)throw std::invalid_argument("queries");
 std::vector<Request> qs(n);for(auto &q:qs){int k;std::cin>>k;if(k<0)throw std::invalid_argument("Z");q.z.resize(k);
 for(int &v:q.z)std::cin>>v;if(!(std::cin>>q.x>>q.y))throw std::invalid_argument("query");}return qs;}
uint64_t rss(){rusage r{};getrusage(RUSAGE_SELF,&r);
#ifdef __APPLE__
 return r.ru_maxrss;
#else
 // /proc reports this executable's address-space high water mark. On Linux,
 // ru_maxrss can retain the Python fork/exec launcher's earlier high water.
 std::ifstream status("/proc/self/status");std::string line;
 while(std::getline(status,line))if(line.rfind("VmHWM:",0)==0){
  std::istringstream value(line.substr(6));uint64_t kib=0;if(value>>kib)return kib*1024;
 }
 return uint64_t(r.ru_maxrss)*1024;
#endif
}
template<class Function>void measure(const std::vector<Request>&qs,double build,Function fn,long long bytes=-1){
 std::cerr<<"phase: build_complete build_ms="<<build<<" query_start\n";
 std::string ans(qs.size(),'0');std::vector<double> lat(qs.size());rusage before{},after{};getrusage(RUSAGE_SELF,&before);auto t=rs::Clock::now();
 for(size_t i=0;i<qs.size();++i){auto s=rs::Clock::now();ans[i]=char('0'+fn(qs[i]));lat[i]=std::chrono::duration<double,std::micro>(rs::Clock::now()-s).count();}
 double query=rs::elapsed(t);getrusage(RUSAGE_SELF,&after);
 std::sort(lat.begin(),lat.end());auto pct=[&](double p){return lat.empty()?0.0:lat[std::min(lat.size()-1,size_t(p*double(lat.size())))];};
 std::cout<<std::setprecision(12)<<"{\"p50_us\":"<<pct(.5)<<",\"p99_us\":"<<pct(.99)<<",";
 std::cerr<<"phase: query_complete query_ms="<<query<<'\n';
 std::cout<<"\"build_ms\":"<<build<<",\"query_ms\":"<<query<<",\"queries\":"<<qs.size()
  <<",\"retained_capacity_bytes\":"<<(bytes<0?"null":std::to_string(bytes))
  <<",\"query_context_switches\":"<<(after.ru_nvcsw-before.ru_nvcsw+after.ru_nivcsw-before.ru_nivcsw)
  <<",\"peak_process_rss_bytes\":"<<rss()<<",\"answers\":\""<<ans<<"\"}\n";
}
void bench(std::string method,int seed,const std::string &graph_path){
 (void)seed;auto g=input_graph(graph_path);auto qs=input_queries();std::cerr<<"phase: build_start method="<<method<<'\n';auto t=rs::Clock::now();
 if(method=="rp-full"||method=="rp-boolean"){
  rs::RPReachable idx(g,method=="rp-boolean");double b=rs::elapsed(t);
  measure(qs,b,[&](const Request&q){return idx.separated(q.x,q.y,q.z);});
 }else throw std::invalid_argument("unknown method (this driver only runs rp-boolean and rp-full)");
}
int main(int argc,char**argv){std::ios::sync_with_stdio(false);std::cin.tie(nullptr);try{
 if(argc==3&&std::string(argv[1])=="probe"){auto t=rs::Clock::now();auto g=input_graph(argv[2]);
  std::cout<<"{\"n\":"<<g.n<<",\"m\":"<<g.m<<",\"load_ms\":"<<rs::elapsed(t)<<",\"peak_process_rss_bytes\":"<<rss()<<"}\n";return 0;}
 if(argc<3||std::string(argv[1])!="bench")throw std::invalid_argument("usage: paper_baselines bench METHOD [seed] [graph.bin]; probe graph.bin");
 bench(argv[2],argc>3?std::stoi(argv[3]):17,argc>4?argv[4]:"");return 0;
 }catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 2;}}
