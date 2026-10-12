#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
struct Row { std::string phase,variant,digest; int pair,order; double rate; };
static void need(bool b) { if(!b) throw std::runtime_error("invalid frozen panel"); }
static double median(std::vector<double> v) { need(v.size()==5); std::sort(v.begin(),v.end()); return v[2]; }
static std::vector<std::tuple<std::string,int,int,std::string>> schedule() {
    std::vector<std::tuple<std::string,int,int,std::string>> s={{"warmup",0,1,"control"},{"warmup",0,2,"candidate"}};
    for(const std::string phase : {"aa","ab"}) for(int pair=1;pair<=5;++pair) {
        std::string a=phase=="aa"?"aa1":"control", b=phase=="aa"?"aa2":"candidate";
        if(pair%2==0) std::swap(a,b);
        s.emplace_back(phase,pair,1,a); s.emplace_back(phase,pair,2,b);
    }
    return s;
}
static void emit(const std::vector<Row>& rows,std::ostream& out,const std::string& preflightSha) {
    const auto s=schedule(); need(rows.size()==s.size());
    std::map<std::pair<std::string,int>,std::map<std::string,double>> rates;
    for(size_t i=0;i<rows.size();++i) {
        const auto& r=rows[i]; need(std::make_tuple(r.phase,r.pair,r.order,r.variant)==s[i]);
        need(std::isfinite(r.rate)&&r.rate>0 && r.digest.size()==64);
        need(r.digest.find_first_not_of("0123456789abcdef")==std::string::npos);
        rates[{r.phase,r.pair}][r.variant]=r.rate;
    }
    std::vector<double> control,candidate,ratios; double logs=0,aaMax=0;
    for(int i=1;i<=5;++i) {
        auto& aa=rates.at({"aa",i}); auto& ab=rates.at({"ab",i});
        const double a=aa.at("aa1"),b=aa.at("aa2"); aaMax=std::max(aaMax,2*std::abs(a-b)/(a+b));
        control.push_back(ab.at("control")); candidate.push_back(ab.at("candidate"));
        ratios.push_back(candidate.back()/control.back()); logs+=std::log(ratios.back());
    }
    const double gm=std::exp(logs/5), cm=median(candidate);
    const bool noise=aaMax<0.01;
    const bool promote=noise&&gm>=1.10&&*std::min_element(ratios.begin(),ratios.end())>1;
    out<<std::setprecision(17)<<"{\n  \"schema\":\"ecc2k130_global_hints.v1\",\n  \"timingPanelValid\":true,\n"
       <<"  \"correctnessPreflightPassed\":true,\n  \"preflightSha256\":\""<<preflightSha<<"\",\n"
       <<"  \"updatesPerSample\":201863462912,\n  \"controlMedianMps\":"<<median(control)
       <<",\n  \"candidateMedianMps\":"<<cm<<",\n  \"geometricMeanRatio\":"<<gm
       <<",\n  \"aaMaximumSymmetricDrift\":"<<aaMax<<",\n  \"pairedRatios\":[";
    for(size_t i=0;i<ratios.size();++i) out<<(i?",":"")<<ratios[i];
    out<<"],\n  \"noiseGate\":"<<(noise?"true":"false")
       <<",\n  \"engineeringGate\":"<<(promote?"true":"false")
       <<",\n  \"rateGoalGate\":"<<(noise&&cm>=26000?"true":"false")
       <<",\n  \"decision\":\""<<(!noise?"INCONCLUSIVE_NOISE":promote?"QUALIFIES_ENGINEERING":"DO_NOT_PROMOTE")
       <<"\",\n  \"scope\":\"bounded complete-update benchmark; correctness gates recorded separately\"\n}\n";
}
int main(int argc,char**argv) {
    try {
        if(argc==2&&std::string(argv[1])=="--self-test") {
            std::vector<Row> rows;
            for(auto& [phase,pair,order,variant]:schedule())
                rows.push_back({phase,variant,std::string(64,'0'),pair,order,variant=="candidate"?1200.0:1000.0});
            std::ostringstream out; emit(rows,out,std::string(64,'0')); need(out.str().find("QUALIFIES_ENGINEERING")!=std::string::npos);
            rows[2].rate=0; bool rejected=false; try {emit(rows,out,std::string(64,'0'));} catch(...) {rejected=true;} need(rejected);
            rows[2].rate=1000; rows.pop_back(); rejected=false; try {emit(rows,out,std::string(64,'0'));} catch(...) {rejected=true;} need(rejected);
            std::cout<<"PASS: schedule, positive panel, invalid rate and missing-row rejection\n"; return 0;
        }
        if(argc!=5) {std::cerr<<"usage: summarize samples.tsv preflight.txt preflight-sha256 result.json\n";return 2;}
        std::ifstream preflight(argv[2]); std::string marker,extra;
        need(bool(std::getline(preflight,marker)) &&
             marker=="PASS native/device queue controls, odd300 replay, spread299 nonzero, full/partial corpus identity and bidirectional checkpoints" &&
             !std::getline(preflight,extra));
        const std::string preflightSha=argv[3]; need(preflightSha.size()==64 &&
            preflightSha.find_first_not_of("0123456789abcdef")==std::string::npos);
        std::ifstream in(argv[1]); need(bool(in)); std::string line;
        need(bool(std::getline(in,line))&&line=="phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState");
        std::vector<Row> rows;
        while(std::getline(in,line)) {
            std::istringstream s(line); std::vector<std::string> f; std::string v;
            while(std::getline(s,v,'\t')) f.push_back(v);
            need(f.size()==7); rows.push_back({f[0],f[3],f[5],std::stoi(f[1]),std::stoi(f[2]),std::stod(f[4])});
        }
        std::ofstream out(argv[4]);need(bool(out));emit(rows,out,preflightSha);need(bool(out));
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
