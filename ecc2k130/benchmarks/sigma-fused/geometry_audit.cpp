#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Offline validation of the frozen geometry evidence. Hash validation stays
// in audit-geometry.sh using the platform SHA-256 command.
static std::string read(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("missing " + path);
    return std::string(std::istreambuf_iterator<char>(in), {});
}
static void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
static void contains(const std::string &s, const std::string &v) {
    require(s.find(v) != std::string::npos, "missing marker: " + v);
}
struct Arm { const char *name; int batch, block, minBlocks, workers; };
static const Arm arms[] = {
    {"b16t256",16,256,2,385024}, {"b32t256",32,256,2,192512},
    {"b32t512",32,512,1,192512}, {"b64t256",64,256,2,96256},
    {"b64t512",64,512,1,96256}
};
static void markers(const std::string &log, const Arm &a, bool verify) {
    contains(log, "packed launch bounds: " + std::to_string(a.block) + " threads, " +
                  std::to_string(a.minBlocks) + " min blocks\n");
    contains(log, "packed sigma fused: 1\n");
    contains(log, "packed sigma fused late y: 0\n");
    contains(log, "packed witness: 0\n");
    contains(log, "packed table walk: 0 (0 branches, 0 shared bytes)\n");
    contains(log, "packed native carryless multiply: 1\n");
    contains(log, "packed compact state: 1\n");
    contains(log, "packed shared sigma: 1\n");
    contains(log, "packed kernel: 126 registers/thread, 0 local bytes/thread, 1792 shared bytes/block");
    const int threads = verify ? a.workers / 4 : a.workers;
    contains(log, "backend cuda-packed131: " + std::to_string(threads) + " threads x " +
        std::to_string(a.batch) + " slots x 1 lanes = " +
        std::string(verify ? "1540096 walks, dp weight 48, 95" : "6160384 walks, dp weight 0, 1024") +
        " steps per launch\n");
    require(log.find("MISMATCH") == std::string::npos && log.find("OVERFLOW") == std::string::npos,
            "mismatch/overflow");
}
static void work(const std::string &log, uint64_t expected) {
    const std::regex progress("([0-9]+) iterations[ ]+([0-9]+) dp[ ]+([0-9]+) stored[ ]+([0-9]+) dropped");
    uint64_t last = 0;
    for (auto i = std::sregex_iterator(log.begin(), log.end(), progress); i != std::sregex_iterator(); ++i) {
        last = std::stoull((*i)[1]);
        require(std::stoull((*i)[4]) == 0, "dropped progress reports");
    }
    require(last == expected, "incorrect final work count");
}
int main(int argc, char **argv) {
    if (argc != 2) { std::cerr << "usage: geometry_audit RESULTS\n"; return 2; }
    try {
        const std::string dir = argv[1];
        require(read(dir + "/exit-code") == "0\n", "job exit");
        contains(read(dir + "/host.txt"), "NVIDIA RTX PRO 6000 Blackwell Server Edition");
        contains(read(dir + "/host.txt"), "V13.3.73");
        for (const Arm &a : arms) {
            const std::string log = read(dir + "/verify-" + a.name + ".log");
            markers(log, a, true); work(log, 1024163840ull);
            contains(log, "1710327 distinguished points (300 verified against the reference, 0 dropped)");
        }
        std::ifstream samples(dir + "/samples.tsv");
        std::string line;
        std::getline(samples,line);
        int warmups=0, screens=0;
        std::array<int,5> perArm{};
        while (std::getline(samples,line)) {
            std::istringstream s(line);
            std::array<std::string,7> cols;
            for (int j=0;j<6;++j) require(bool(std::getline(s,cols[j],'\t')), "sample framing");
            std::getline(s,cols[6]);
            int idx=-1;
            for (int j=0;j<5;++j) if(cols[3]==arms[j].name) idx=j;
            require(idx>=0, "unknown arm");
            const bool warm = cols[0]=="warmup";
            require(warm || cols[0]=="screen", "unknown phase");
            const std::string file=dir+"/"+cols[0]+"-"+cols[1]+"-"+cols[2]+"-"+cols[3]+".log";
            const std::string log=read(file);
            markers(log,arms[idx],false);
            work(log,warm ? 100931731456ull : 201863462912ull);
            const std::regex final("finished: ([0-9.]+) M it/s, 0 distinguished points \\(0 verified against the reference, 0 dropped\\)");
            std::smatch m;
            require(std::regex_search(log,m,final), "final rate line");
            require(std::stod(m[1])==std::stod(cols[4]), "rate differs from TSV");
            if(warm) ++warmups; else { ++screens; ++perArm[idx]; }
        }
        require(warmups==5 && screens==15, "incomplete timing panel");
        for(int n:perArm) require(n==3,"incomplete arm");
        std::ifstream corpus(dir+"/corpus-sorted.bin",std::ios::binary);
        require(bool(corpus),"missing sorted corpus");
        std::array<unsigned char,32> record{}, previous{};
        uint64_t records=0;
        while(corpus.read(reinterpret_cast<char*>(record.data()),32)) {
            if(records) require(previous<=record,"corpus not sorted");
            // v1 seed's high 16 bits carry run id 31. Canonical x is 131 bits.
            require(record[6]==31 && record[7]==0,"foreign run id");
            require(record[24]<=7,"noncanonical high limb");
            for(int k=25;k<32;++k) require(record[k]==0,"noncanonical upper bits");
            previous=record; ++records;
        }
        require(corpus.eof() && corpus.gcount()==0 && records==1710327,"corpus framing/count");
        std::cout << "{\"valid\":true,\"verification_arms\":5,\"replayed_per_arm\":300,"
                     "\"warmup_rows\":5,\"screen_rows\":15,\"screen_updates_per_row\":201863462912,"
                     "\"sorted_records\":1710327,\"run_id\":31,\"registers_per_thread\":126,"
                     "\"local_bytes_per_thread\":0,\"shared_bytes_per_block\":1792}\n";
        return 0;
    } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
