// Thin construction-only entry point. The build script mechanically extracts
// the unchanged frozen backend before its main(), so the arithmetic is shared.
#include "cm-orbit-kernel.hpp"
int main(int argc,char**argv) {
  if(argc!=5) return 2;
  std::ifstream parameters(argv[2]);CMParameters cp;
  while(parameters>>cp.p>>cp.r>>cp.pair>>cp.lambda>>cp.inv_c>>cp.inv_cu>>cp.x[0]>>cp.x[1]>>cp.t[0]>>cp.t[1]>>cp.y2[0]>>cp.y2[1])
    if(!all_cm_parameters.emplace(cp.pair,cp).second) throw std::runtime_error("duplicate parameters");
  if(!parameters.eof()||all_cm_parameters.empty()) throw std::runtime_error("bad parameters");
  std::ifstream input(argv[1]);std::string line;std::getline(input,line);
  std::vector<Job> jobs;
  while(std::getline(input,line)) {
    std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);Job w;
    s>>w.id>>w.pair>>w.role>>w.fp>>w.order>>w.a>>w.b>>w.gx>>w.gy>>w.qx>>w.qy>>w.seed>>w.repeat;
    if(!s) throw std::runtime_error("bad workload");jobs.push_back(w);
  }
  U pass=std::stoull(argv[3]);std::mt19937_64 rng(20261020+pass);
  std::shuffle(jobs.begin(),jobs.end(),rng);
  // One worker thread, matching the frozen thread-CPU measurement boundary.
  std::thread worker([&]{for(auto w:jobs){w.repeat=pass;run(w,8192,1,-1,argv[4]);}});
  worker.join();
  return 0;
}
