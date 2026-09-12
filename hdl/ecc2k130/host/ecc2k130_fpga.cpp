// ecc2k130_fpga.cpp -- host program for the FPGA rho engine (ec2k_axil).
//
// A drop-in replacement for the GPU client under ecc2k130/aws/worker.py: same
// command line, same progress lines, same dp file (32-byte (seed, canonical
// x) records), same checkpoint header, same exit codes.  What differs is
// where the steps happen.  The FPGA holds NENG * 2**ID_W walks and reports
// (walk id, step count, point) whenever one lands on a distinguished point;
// this program owns everything else, exactly as the GPU client's host side
// does: seeds and start points, the corpus file, in-process collision
// detection, verification and checkpoints.
//
//   ecc2k130-fpga --curve 131 --run-id 7 --dp-file dp.bin --checkpoint walk.ck \
//                 --checkpoint-every 600 --verify 0 --device 0
//
//   --device D        FPGA slot (F2: 0..7); ECC_GPU is what worker.py sets
//   --sim             software model of the register block instead of a slot
//                     (built into every binary; the only backend on machines
//                     without the F2 SDK)
//   --dp-weight W     must match the weight the bitstream was built with
//   --launches L      stop after L polls of the device (0 = until solved)
//   --dps N           stop after N distinguished points (0 = no limit)
//   --selftest        check the fast host arithmetic against the reference
//
// Exit codes match the GPU client where worker.py looks at them:
//   0 finished / solved, 1 usage, 3 --verify mismatch, 5 parameter check,
//   6 checkpoint refused, 7 device not found or wrong bitstream.
//
// Seeds.  A walk with flat id g on run r walks seed (r<<48)|(e<<32)|(g<<16)|c:
// c counts restarts of that id, e is an epoch chosen when the run starts
// fresh (so a lost checkpoint does not replay old seeds) and bumped when c
// wraps.  Restarting from a checkpoint continues with every saved seed + 1:
// the walks in flight at the time are gone (the engine's state is not
// readable) and their seeds are simply never reported.
//
// Distinguished start points.  The engine tests the weight after each step,
// never of the point it was loaded with, so a start point that is itself
// distinguished (probability 2**-26 per walk) is reported here directly with
// zero steps, which is what the reference walk returns for that seed.

#include "../../../ecc2k130/include/curveparams.h"
#include "../../../ecc2k130/include/packed131.h"
#include "../../../ecc2k130/include/solver.h"

#include <algorithm>
#include <deque>
#include <errno.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#ifdef ECC_FPGA_PCI
#include <fpga_mgmt.h>
#include <fpga_pci.h>
#endif

typedef Ref<CfgF131> R;
typedef R::Elem Elem;
typedef eccPacked131::P131 P131;
using eccPacked131::add131;
using eccPacked131::inv131;
using eccPacked131::mul131;
using eccPacked131::sigma131;
using eccPacked131::sqr131;

static const int M = 131;

// ---------------------------------------------------------------------------
// register map of ec2k_axil.vhd
// ---------------------------------------------------------------------------
namespace reg {
enum : uint32_t {
    MAGIC = 0x000, CTRL = 0x004, STATUS = 0x008, GEOM = 0x00C,
    STEPS_LO = 0x010, STEPS_HI = 0x014, DPS = 0x018, DROPPED = 0x01C,
    LD_ID = 0x020, LD_X = 0x024, LD_Y = 0x038, LD_GO = 0x04C, CLOCK = 0x050,
    DP_ID = 0x080, DP_STEPS_LO = 0x084, DP_STEPS_HI = 0x088,
    DP_X = 0x090, DP_Y = 0x0A4, DP_POP = 0x0B8,
};
static const uint32_t MAGIC_VALUE = 0x2C130001u;
static const uint32_t CTRL_RUN = 1u, CTRL_CLEAR = 2u;
static const uint32_t ST_LD_BUSY = 1u, ST_DP_AVAIL = 2u;   // bit 2 (was DP_OVERFLOW) reads 0
}  // namespace reg

struct Geometry {
    int idW, logW, logNb, dpWeight, neng;
    unsigned walks() const { return (unsigned)neng << idW; }
    static Geometry decode(uint32_t g) {
        Geometry o;
        o.idW = (int)(g & 0xFF);
        o.logW = (int)((g >> 8) & 0xF);
        o.logNb = (int)((g >> 12) & 0xF);
        o.dpWeight = (int)((g >> 16) & 0xFF);
        o.neng = (int)((g >> 24) & 0xFF);
        return o;
    }
    uint32_t encode() const {
        return (uint32_t)idW | ((uint32_t)logW << 8) | ((uint32_t)logNb << 12) |
               ((uint32_t)dpWeight << 16) | ((uint32_t)neng << 24);
    }
};

// ---------------------------------------------------------------------------
// field helpers on the packed (register image) representation
// ---------------------------------------------------------------------------
static P131 pack(const Elem &a) {
    P131 p;
    for (int i = 0; i < 5; i++) p.v[i] = (uint32_t)(a.v[i / 2] >> (32 * (i & 1)));
    return p;
}
static Elem unpack(const P131 &p) {
    unsigned long long a[3] = {p.v[0] | ((unsigned long long)p.v[1] << 32),
                               p.v[2] | ((unsigned long long)p.v[3] << 32), p.v[4]};
    return R::fromLimbs(a);
}
static bool same(const P131 &a, const P131 &b) {
    for (int i = 0; i < 5; i++)
        if (a.v[i] != b.v[i]) return false;
    return true;
}
static int weight(const P131 &a) {
    int w = 0;
    for (int i = 0; i < 5; i++) w += __builtin_popcount(a.v[i]);
    return w;
}

// Affine addition with no special cases: the formula the engine implements.
static void addRaw(P131 &x, P131 &y, const P131 &sx, const P131 &sy) {
    const P131 d = add131(x, sx);
    const P131 lam = mul131(add131(y, sy), inv131(d));
    const P131 x3 = add131(add131(sqr131(lam), lam), d);
    const P131 y3 = add131(add131(mul131(lam, add131(x, x3)), x3), y);
    x = x3;
    y = y3;
}

static int jOf(int hw) { return 3 + ((hw >> 1) & 7); }

static void stepPacked(P131 &x, P131 &y, int hw) {
    const int j = jOf(hw);
    const P131 sx = sigma131(x, j), sy = sigma131(y, j);
    addRaw(x, y, sx, sy);
}

// sigma^i(P) for i < 131 and Q, so a start point is 128 conditional additions
struct StartTable {
    P131 px[131], py[131], qx, qy;
    void setup(const Solver<CfgF131> &sol) {
        for (int i = 0; i < M; ++i) {
            px[i] = pack(R::sigma(sol.basis.x, i));
            py[i] = pack(R::sigma(sol.basis.y, i));
        }
        qx = pack(sol.target.x);
        qy = pack(sol.target.y);
    }
    void startPoint(u64 seed, P131 &x, P131 &y) const {
        const u64 c0 = R::eccPrfHost(seed, 0), c1 = R::eccPrfHost(seed, 1);
        x = qx;
        y = qy;
        for (int i = 0; i < 128; ++i) {
            const u64 b = (i < 64) ? (c0 >> i) : (c1 >> (i - 64));
            if (b & 1) addRaw(x, y, px[i % M], py[i % M]);
        }
    }
};

struct FastWalk {
    bool ok;
    P131 x, y;
    u64 iters;
    unsigned long long counts[8];
};

static FastWalk rewalkPacked(const StartTable &T, u64 seed, int dpWeight, u64 maxIters) {
    FastWalk w;
    w.ok = false;
    w.iters = 0;
    memset(w.counts, 0, sizeof w.counts);
    T.startPoint(seed, w.x, w.y);
    for (u64 it = 0;; ++it) {
        const int hw = weight(w.x);
        if (hw <= dpWeight) {
            w.ok = true;
            w.iters = it;
            return w;
        }
        if (it >= maxIters) return w;
        w.counts[jOf(hw) - 3]++;
        stepPacked(w.x, w.y, hw);
    }
}

// ---------------------------------------------------------------------------
// bus: the register block, real or modelled
// ---------------------------------------------------------------------------
struct Bus {
    virtual ~Bus() {}
    virtual uint32_t peek(uint32_t off) = 0;
    virtual void poke(uint32_t off, uint32_t v) = 0;
};

#ifdef ECC_FPGA_PCI
struct PciBus : Bus {
    pci_bar_handle_t h = PCI_BAR_HANDLE_INIT;
    bool open(int slot) {
        if (fpga_mgmt_init() || fpga_pci_init()) {
            fprintf(stderr, "fpga_mgmt/fpga_pci init failed (is the F2 SDK loaded?)\n");
            return false;
        }
        if (fpga_pci_attach(slot, FPGA_APP_PF, APP_PF_BAR0, 0, &h)) {
            fprintf(stderr, "cannot attach to FPGA slot %d BAR0\n", slot);
            return false;
        }
        return true;
    }
    ~PciBus() {
        if (h != PCI_BAR_HANDLE_INIT) fpga_pci_detach(h);
    }
    uint32_t peek(uint32_t off) override {
        uint32_t v = 0;
        if (fpga_pci_peek(h, off, &v)) { fprintf(stderr, "pci peek 0x%x failed\n", off); exit(7); }
        return v;
    }
    void poke(uint32_t off, uint32_t v) override {
        if (fpga_pci_poke(h, off, v)) { fprintf(stderr, "pci poke 0x%x failed\n", off); exit(7); }
    }
};
#endif

// Behavioural model of ec2k_axil: same registers, same handshakes, walks
// advanced in software each time STATUS is read.  Loads complete at once;
// the DP queue is 64 deep and, as in hardware, a report that finds it full
// waits in its engine (here: a per-walk pending flag) rather than being
// lost, so DROPPED and STATUS.DP_OVERFLOW are always zero.
struct SimBus : Bus {
    Geometry geo;
    unsigned stepsPerPoll;
    struct Walk { bool live = false, pending = false; P131 x, y; uint32_t steps = 0; };
    struct Dp { uint32_t gid, steps; P131 x, y; };
    std::vector<Walk> walks;
    std::deque<Dp> fifo;
    uint32_t ctrl = 0, ldId = 0, ldx[5] = {0}, ldy[5] = {0};
    uint64_t steps = 0;
    uint32_t stepsHi = 0, dps = 0;

    SimBus(const Geometry &g, unsigned spp) : geo(g), stepsPerPoll(spp), walks(g.walks()) {}

    void report(unsigned g) {
        Walk &w = walks[g];
        if (fifo.size() >= 64) { w.pending = true; return; }
        fifo.push_back(Dp{g, w.steps, w.x, w.y});
        w.pending = false;
        dps++;
    }
    void advance() {
        if (!(ctrl & reg::CTRL_RUN)) return;
        for (unsigned g = 0; g < walks.size(); ++g)
            if (walks[g].pending) report(g);
        for (unsigned g = 0; g < walks.size(); ++g) {
            Walk &w = walks[g];
            if (!w.live) continue;
            for (unsigned s = 0; s < stepsPerPoll && w.live; ++s) {
                stepPacked(w.x, w.y, weight(w.x));
                w.steps++;
                steps++;
                if (weight(w.x) <= geo.dpWeight) {
                    w.live = false;
                    report(g);
                }
            }
        }
    }
    uint32_t peek(uint32_t off) override {
        using namespace reg;
        if (off >= LD_X && off < LD_X + 20) return ldx[(off - LD_X) / 4];
        if (off >= LD_Y && off < LD_Y + 20) return ldy[(off - LD_Y) / 4];
        if (off >= DP_X && off < DP_X + 20) return fifo.empty() ? 0 : fifo.front().x.v[(off - DP_X) / 4];
        if (off >= DP_Y && off < DP_Y + 20) return fifo.empty() ? 0 : fifo.front().y.v[(off - DP_Y) / 4];
        switch (off) {
        case MAGIC: return MAGIC_VALUE;
        case CTRL: return ctrl & CTRL_RUN;
        case STATUS:
            advance();
            return (fifo.empty() ? 0 : ST_DP_AVAIL) |
                   ((uint32_t)std::min<size_t>(fifo.size(), 255) << 16);
        case GEOM: return geo.encode();
        case CLOCK: return 333333;   // kHz, as an image built with the defaults reports
        case STEPS_LO: stepsHi = (uint32_t)(steps >> 32); return (uint32_t)steps;
        case STEPS_HI: return stepsHi;
        case DPS: return dps;
        case DROPPED: return 0;
        case LD_ID: return ldId;
        case DP_ID: return fifo.empty() ? 0 : fifo.front().gid;
        case DP_STEPS_LO: return fifo.empty() ? 0 : fifo.front().steps;
        case DP_STEPS_HI: return 0;
        default: return 0;
        }
    }
    void poke(uint32_t off, uint32_t v) override {
        using namespace reg;
        if (off >= LD_X && off < LD_X + 20) { ldx[(off - LD_X) / 4] = ((off - LD_X) / 4 == 4) ? (v & 7) : v; return; }
        if (off >= LD_Y && off < LD_Y + 20) { ldy[(off - LD_Y) / 4] = ((off - LD_Y) / 4 == 4) ? (v & 7) : v; return; }
        switch (off) {
        case CTRL:
            ctrl = v & CTRL_RUN;
            if (!(v & CTRL_RUN)) for (Walk &w : walks) w.live = w.pending = false;
            if (v & CTRL_CLEAR) { steps = 0; dps = 0; fifo.clear(); }
            break;
        case LD_ID: ldId = v; break;
        case LD_GO:
            if ((ctrl & CTRL_RUN) && ldId < walks.size()) {
                Walk &w = walks[ldId];
                memcpy(w.x.v, ldx, sizeof ldx);
                memcpy(w.y.v, ldy, sizeof ldy);
                w.steps = 0;
                w.live = true;
            }
            break;
        case DP_POP: if (!fifo.empty()) fifo.pop_front(); break;
        default: break;
        }
    }
};

// ---------------------------------------------------------------------------
// files shared with the GPU client
// ---------------------------------------------------------------------------
struct DpFileRecord {
    unsigned long long seed;
    unsigned long long canon[3];
};

struct CkptHeader {
    char magic[8];
    unsigned version;
    unsigned m;
    unsigned threads;   // walks on the device
    unsigned batch;     // 1
    unsigned lanes;     // 1
    unsigned runId;
    unsigned long long iterBase;   // steps completed, all walks
};
static const unsigned CKPT_VERSION = 0xF2u;

static double nowSeconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static volatile sig_atomic_t gStop = 0;
static void onStop(int) { gStop = 1; }

struct Options {
    int curve = 131;
    int dpWeight = -1;
    unsigned runId = 1;
    int verify = 8;
    long launches = 0;
    unsigned long long dps = 0;
    u64 maxIters = 0;
    int device = 0;
    bool sim = false;
    int simEng = 2, simIdW = 4, simDpWeight = -1;   // the model's "bitstream" weight
    unsigned simSteps = 1;
    bool selftest = false;
    double pollMs = 5.0;
    std::string dpFile, ckptFile;
    std::vector<std::string> loadFiles;
    unsigned long long loadMax = 0;
    double ckptSeconds = 300.0;
};

// ---------------------------------------------------------------------------
// seeds
// ---------------------------------------------------------------------------
static u64 seedOf(unsigned runId, unsigned epoch, unsigned gid, unsigned ctr) {
    return ((u64)runId << 48) | ((u64)(epoch & 0xFFFF) << 32) | ((u64)(gid & 0xFFFF) << 16) | (ctr & 0xFFFF);
}
static u64 nextSeed(u64 s) {
    const unsigned ctr = (unsigned)(s & 0xFFFF);
    if (ctr == 0xFFFF) {
        const unsigned epoch = (unsigned)((s >> 32) & 0xFFFF) + 1;
        return (s & 0xFFFF0000FFFF0000ull) | ((u64)(epoch & 0xFFFF) << 32);
    }
    return s + 1;
}

// ---------------------------------------------------------------------------
// the run
// ---------------------------------------------------------------------------
struct Client {
    Options o;
    Bus *bus = 0;
    Geometry geo;
    unsigned nw = 0;
    Solver<CfgF131> sol;
    StartTable T;
    std::vector<u64> seeds;
    u64 totalDp = 0, dpBefore = 0, verified = 0, verifyBudget = 0;
    FILE *dpOut = 0;

    void writeElem(uint32_t base, const P131 &v) {
        for (int k = 0; k < 5; ++k) bus->poke(base + 4 * k, v.v[k]);
    }
    P131 readElem(uint32_t base) {
        P131 v;
        for (int k = 0; k < 5; ++k) v.v[k] = bus->peek(base + 4 * k);
        return v;
    }

    // Load a fresh walk into id g.  A start point that is already
    // distinguished is reported at once and the next seed tried.
    bool startWalk(unsigned g, std::vector<DpRecord> &out) {
        for (;;) {
            P131 x, y;
            T.startPoint(seeds[g], x, y);
            if (weight(x) > geo.dpWeight) {
                for (int spins = 0; bus->peek(reg::STATUS) & reg::ST_LD_BUSY; ++spins) {
                    if (spins > 1000000) { fprintf(stderr, "engine never accepted a load\n"); return false; }
                }
                bus->poke(reg::LD_ID, g);
                writeElem(reg::LD_X, x);
                writeElem(reg::LD_Y, y);
                bus->poke(reg::LD_GO, 1);
                return true;
            }
            DpRecord rec;
            rec.seed = seeds[g];
            rec.iters = 0;
            const Elem ex = unpack(x), ey = unpack(y);
            memcpy(rec.x, ex.v, sizeof rec.x);
            memcpy(rec.y, ey.v, sizeof rec.y);
            out.push_back(rec);
            seeds[g] = nextSeed(seeds[g]);
        }
    }

    // Take every queued report off the device.
    void drain(std::vector<DpRecord> &out, unsigned limit) {
        while (limit-- && (bus->peek(reg::STATUS) & reg::ST_DP_AVAIL)) {
            const uint32_t gid = bus->peek(reg::DP_ID);
            const u64 steps = (u64)bus->peek(reg::DP_STEPS_LO) | ((u64)bus->peek(reg::DP_STEPS_HI) << 32);
            const P131 x = readElem(reg::DP_X), y = readElem(reg::DP_Y);
            bus->poke(reg::DP_POP, 1);
            if (gid >= nw) {
                fprintf(stderr, "device reported walk %u of %u; bitstream and host disagree\n", gid, nw);
                exit(7);
            }
            DpRecord rec;
            rec.seed = seeds[gid];
            rec.iters = steps;
            const Elem ex = unpack(x), ey = unpack(y);
            memcpy(rec.x, ex.v, sizeof rec.x);
            memcpy(rec.y, ey.v, sizeof rec.y);
            out.push_back(rec);
            seeds[gid] = nextSeed(seeds[gid]);
            startWalk(gid, out);
        }
    }

    bool saveCheckpoint(u64 steps) const {
        if (o.ckptFile.empty()) return true;
        const std::string tmp = o.ckptFile + ".tmp";
        FILE *f = fopen(tmp.c_str(), "wb");
        if (!f) return false;
        CkptHeader h;
        memset(&h, 0, sizeof h);
        memcpy(h.magic, "ECC2K130", 8);
        h.version = CKPT_VERSION;
        h.m = M;
        h.threads = nw;
        h.batch = 1;
        h.lanes = 1;
        h.runId = o.runId;
        h.iterBase = steps;
        bool ok = fwrite(&h, sizeof h, 1, f) == 1;
        ok = ok && fwrite(seeds.data(), sizeof(u64), seeds.size(), f) == seeds.size();
        const u64 lifetime = dpBefore + totalDp;
        ok = ok && fwrite(&lifetime, sizeof lifetime, 1, f) == 1;
        ok = ok && fflush(f) == 0;
        ok = (fclose(f) == 0) && ok;
        if (ok) ok = rename(tmp.c_str(), o.ckptFile.c_str()) == 0;
        return ok;
    }

    // 1 restored, 0 no checkpoint, -1 refused
    int restoreCheckpoint(u64 *steps) {
        if (o.ckptFile.empty()) return 0;
        FILE *f = fopen(o.ckptFile.c_str(), "rb");
        if (!f) return errno == ENOENT ? 0 : -1;
        CkptHeader h;
        bool ok = fread(&h, sizeof h, 1, f) == 1 && memcmp(h.magic, "ECC2K130", 8) == 0 &&
                  h.version == CKPT_VERSION && h.m == (unsigned)M && h.threads == nw &&
                  h.runId == o.runId;
        if (ok) {
            fseek(f, 0, SEEK_END);
            ok = (unsigned long long)ftell(f) == sizeof h + sizeof(u64) * nw + sizeof(u64);
            fseek(f, sizeof h, SEEK_SET);
        }
        ok = ok && fread(seeds.data(), sizeof(u64), nw, f) == nw;
        ok = ok && fread(&dpBefore, sizeof dpBefore, 1, f) == 1;
        fclose(f);
        if (!ok) return -1;
        for (unsigned g = 0; g < nw; ++g) seeds[g] = nextSeed(seeds[g]);   // those walks are gone
        *steps = h.iterBase;
        return 1;
    }

    // Reload corpus files so a collision against earlier work is found here
    // and not only by the offline merge.  Returns 0 to keep going, otherwise
    // the exit code (a solution found while reloading ends the run).
    int reloadCorpus() {
        std::vector<std::string> corpus = o.loadFiles;
        if (!o.dpFile.empty()) corpus.push_back(o.dpFile);
        for (size_t i = 0; i < corpus.size(); ++i) {
            char *real = realpath(corpus[i].c_str(), 0);
            if (!real) continue;
            corpus[i] = real;
            free(real);
        }
        std::sort(corpus.begin(), corpus.end());
        corpus.erase(std::unique(corpus.begin(), corpus.end()), corpus.end());
        std::sort(corpus.begin(), corpus.end(), [](const std::string &a, const std::string &b) {
            struct stat sa, sb;
            const bool oa = stat(a.c_str(), &sa) == 0, ob = stat(b.c_str(), &sb) == 0;
            if (!oa || !ob) return oa > ob;
            return sa.st_mtime > sb.st_mtime;
        });
        size_t reloaded = 0, skipped = 0;
        for (size_t ci = 0; ci < corpus.size(); ++ci) {
            if (o.loadMax && reloaded >= o.loadMax) { skipped = corpus.size() - ci; break; }
            FILE *in = fopen(corpus[ci].c_str(), "rb");
            if (!in) continue;
            if (o.loadMax && fseek(in, 0, SEEK_END) == 0) {
                const long sz = ftell(in);
                unsigned long long skip = 0;
                if (sz > 0) {
                    const unsigned long long nrec = (unsigned long long)sz / sizeof(DpFileRecord);
                    const unsigned long long remain = o.loadMax - reloaded;
                    if (nrec > remain) skip = nrec - remain;
                }
                if (fseek(in, (long)(skip * sizeof(DpFileRecord)), SEEK_SET) != 0) rewind(in);
            }
            DpFileRecord fr;
            while (fread(&fr, sizeof fr, 1, in) == 1) {
                if (o.loadMax && reloaded >= o.loadMax) break;
                Solver<CfgF131>::Key key;
                key.v[0] = fr.canon[0];
                key.v[1] = fr.canon[1];
                key.v[2] = fr.canon[2];
                Solver<CfgF131>::Entry other;
                ++reloaded;
                if (!sol.insertKey(key, fr.seed, 0, &other)) continue;
                printf("collision found while reloading %s: seeds %016llx and %016llx\n",
                       corpus[ci].c_str(), (unsigned long long)fr.seed, (unsigned long long)other.seed);
                if (solveCollision(fr.seed, other.seed)) {
                    fclose(in);
                    return -1;
                }
            }
            fclose(in);
        }
        if (reloaded) {
            printf("reloaded %zu points from %zu file(s), %zu distinct orbits\n", reloaded,
                   corpus.size() - skipped, sol.inserted);
            if (o.loadMax && reloaded >= o.loadMax)
                printf("  stopped at the --load-max %llu cap; %zu file(s) not read, "
                       "their collisions are left to the offline merge\n",
                       o.loadMax, skipped);
        }
        return 0;
    }

    // Recompute both walks (fast arithmetic for the steps, the reference for
    // the tracked scalar of the start point) and solve.  Prints "k = " on
    // success, which is the line worker.py watches for.
    bool solveCollision(u64 seedA, u64 seedB) {
        const double tr = nowSeconds();
        Solver<CfgF131>::WalkResult W[2];
        const u64 sd[2] = {seedA, seedB};
        for (int i = 0; i < 2; ++i) {
            const FastWalk f = rewalkPacked(T, sd[i], geo.dpWeight, sol.maxIters);
            W[i].ok = f.ok;
            W[i].seed = sd[i];
            W[i].iters = f.iters;
            memcpy(W[i].counts, f.counts, sizeof f.counts);
            W[i].endPoint = R::make(unpack(f.x), unpack(f.y));
            R::startPoint(sd[i], sol.basis, sol.target, &W[i].alpha0, sol.ell, sol.spow);
        }
        U192 k;
        std::string why;
        if (!sol.solve(W[0], W[1], &k, &why)) {
            printf("  unusable (%s), continuing\n", why.c_str());
            return false;
        }
        printf("  recomputed both walks in %.2f s\n", nowSeconds() - tr);
        printf("  k = %s\n", u192_to_dec(k).c_str());
        printf("  verified [k]P == Q\n");
        fflush(stdout);
        return true;
    }

    // 0 keep going, else exit code
    int handleRecord(const DpRecord &rec) {
        if (verifyBudget) {
            --verifyBudget;
            const FastWalk w = rewalkPacked(T, rec.seed, geo.dpWeight, sol.maxIters);
            const bool ok = w.ok && w.iters == rec.iters && same(w.x, pack(R::fromLimbs(rec.x))) &&
                            same(w.y, pack(R::fromLimbs(rec.y)));
            if (!ok) {
                printf("MISMATCH: seed %016llx was not reproduced by the reference walk\n",
                       (unsigned long long)rec.seed);
                return 3;
            }
            ++verified;
        }
        if (dpOut) {
            const Elem cx = R::canonical(R::fromLimbs(rec.x));
            DpFileRecord fr;
            fr.seed = rec.seed;
            fr.canon[0] = cx.v[0];
            fr.canon[1] = cx.v[1];
            fr.canon[2] = cx.v[2];
            fwrite(&fr, sizeof fr, 1, dpOut);
        }
        Solver<CfgF131>::Entry other;
        if (!sol.insert(rec, &other)) return 0;
        printf("collision: seeds %016llx and %016llx meet after %llu and %llu steps\n",
               (unsigned long long)rec.seed, (unsigned long long)other.seed,
               (unsigned long long)rec.iters, (unsigned long long)other.iters);
        return solveCollision(rec.seed, other.seed) ? -1 : 0;
    }

    int run() {
        sol.setup(eccF131::PX, eccF131::PY, eccF131::QX, eccF131::QY, eccF131::ELL_DEC, eccF131::S_DEC,
                  o.dpWeight < 0 ? eccF131::DP_WEIGHT : o.dpWeight,
                  o.maxIters ? o.maxIters : (u64)1 << 40);
        std::string why;
        if (!sol.checkSetup(&why)) {
            printf("parameter check failed: %s\n", why.c_str());
            return 5;
        }
        T.setup(sol);
        if (o.selftest) return selftest();

        // the device
        if (o.sim) {
            Geometry g;
            g.idW = o.simIdW; g.logW = 4; g.logNb = 3; g.neng = o.simEng;
            g.dpWeight = o.simDpWeight < 0 ? sol.dpWeight : o.simDpWeight;
            bus = new SimBus(g, o.simSteps);
        } else {
#ifdef ECC_FPGA_PCI
            PciBus *p = new PciBus;
            if (!p->open(o.device)) return 7;
            bus = p;
#else
            fprintf(stderr, "built without the F2 SDK: only --sim is available\n");
            return 7;
#endif
        }
        if (bus->peek(reg::MAGIC) != reg::MAGIC_VALUE) {
            fprintf(stderr, "slot %d does not hold the ECC2K-130 image (MAGIC 0x%08x)\n", o.device,
                    bus->peek(reg::MAGIC));
            return 7;
        }
        geo = Geometry::decode(bus->peek(reg::GEOM));
        nw = geo.walks();
        if (geo.dpWeight != sol.dpWeight) {
            fprintf(stderr, "bitstream was built for dp weight %d, run asked for %d\n", geo.dpWeight,
                    sol.dpWeight);
            return 7;
        }
        if (nw == 0 || nw > 0x10000) {
            fprintf(stderr, "unsupported geometry: %u walks\n", nw);
            return 7;
        }
        if (o.runId > 0xFFFF) {
            fprintf(stderr, "--run-id must fit 16 bits\n");
            return 1;
        }
        const uint32_t clkKhz = bus->peek(reg::CLOCK);   // 0 from images older than the register
        printf("backend %s: %d engine(s) x %u walks = %u walks, dp weight %d, batches of %d, %d in flight, "
               "engine clock %s\n",
               o.sim ? "fpga-sim" : "fpga", geo.neng, 1u << geo.idW, nw, geo.dpWeight, 1 << geo.logW,
               1 << geo.logNb,
               clkKhz ? (std::to_string(clkKhz / 1000) + "." + std::to_string(clkKhz / 100 % 10) + " MHz").c_str()
                      : "not reported (250 MHz shell clock)");

        // corpus, then walk state
        if (int rc = reloadCorpus()) return rc < 0 ? 0 : rc;
        seeds.resize(nw);
        u64 steps0 = 0;
        const int rs = restoreCheckpoint(&steps0);
        if (rs < 0) {
            fprintf(stderr, "checkpoint %s is incompatible or incomplete; use the matching backend/settings or a new checkpoint path\n",
                    o.ckptFile.c_str());
            return 6;
        }
        if (rs > 0) {
            printf("resumed from %s at iteration %llu\n", o.ckptFile.c_str(), (unsigned long long)steps0);
        } else {
            if (!o.ckptFile.empty()) printf("no checkpoint at %s, starting fresh\n", o.ckptFile.c_str());
            const unsigned epoch = (unsigned)((time(0) / 60) & 0xFFFF);
            for (unsigned g = 0; g < nw; ++g) seeds[g] = seedOf(o.runId, epoch, g, 0);
        }
        verifyBudget = (u64)(o.verify < 0 ? 0 : o.verify);
        dpOut = o.dpFile.empty() ? 0 : fopen(o.dpFile.c_str(), "ab");

        // reset the engines and seed every walk
        bus->poke(reg::CTRL, reg::CTRL_CLEAR);
        bus->poke(reg::CTRL, reg::CTRL_RUN | reg::CTRL_CLEAR);
        std::vector<DpRecord> recs;
        const double tl = nowSeconds();
        for (unsigned g = 0; g < nw; ++g)
            if (!startWalk(g, recs)) return 7;
        printf("loaded %u walks in %.2f s\n", nw, nowSeconds() - tl);
        fflush(stdout);

        signal(SIGINT, onStop);
        signal(SIGTERM, onStop);

        const double t0 = nowSeconds();
        double lastPrint = t0, lastCkpt = t0;
        u64 lost = 0, devSteps = 0;
        bool warnedLost = false;
        int rc = 0;
        for (long poll = 0; o.launches == 0 || poll < o.launches; ++poll) {
            drain(recs, 4096);
            for (size_t i = 0; i < recs.size() && rc == 0; ++i) {
                ++totalDp;
                rc = handleRecord(recs[i]);
            }
            recs.clear();
            if (rc) break;
            if (o.dps && totalDp >= o.dps) break;

            const double now = nowSeconds();
            const bool leaving = gStop || (o.launches && poll + 1 == o.launches);
            if (leaving || now - lastPrint > 2.0 || now - lastCkpt > o.ckptSeconds) {
                devSteps = steps0 + ((u64)bus->peek(reg::STEPS_LO) | ((u64)bus->peek(reg::STEPS_HI) << 32));
                lost = bus->peek(reg::DROPPED);
            }
            if (dpOut && (leaving || now - lastCkpt > o.ckptSeconds)) fflush(dpOut);
            if (!o.ckptFile.empty() && (leaving || now - lastCkpt > o.ckptSeconds)) {
                if (!saveCheckpoint(devSteps)) printf("warning: could not write checkpoint %s\n", o.ckptFile.c_str());
                lastCkpt = now;
            }
            if (gStop) {
                printf("stopping: %llu iterations of %u walks, %llu points reported\n",
                       (unsigned long long)(devSteps - steps0), nw, (unsigned long long)totalDp);
                break;
            }
            if (now - lastPrint > 2.0 || leaving) {
                const double el = now - t0;
                const double it = (double)(devSteps - steps0);
                printf("  %8.1f s  %10.3f M it/s  %10llu iterations  %8llu dp  %8llu stored  %8llu dropped\n",
                       el, el > 0 ? it / el / 1e6 : 0.0, (unsigned long long)it, (unsigned long long)totalDp,
                       (unsigned long long)sol.inserted, (unsigned long long)lost);
                if (lost && !warnedLost && lost * 10 > totalDp + lost) {
                    printf("  warning: dropping %.1f%% of reports; the host is not polling the queue fast enough\n",
                           100.0 * (double)lost / (double)(totalDp + lost));
                    warnedLost = true;
                }
                fflush(stdout);
                lastPrint = now;
            }
            if (!o.sim && o.pollMs > 0) usleep((useconds_t)(o.pollMs * 1000));
        }

        devSteps = steps0 + ((u64)bus->peek(reg::STEPS_LO) | ((u64)bus->peek(reg::STEPS_HI) << 32));
        lost = bus->peek(reg::DROPPED);
        if (dpOut) { fflush(dpOut); fclose(dpOut); dpOut = 0; }
        if (!o.ckptFile.empty()) saveCheckpoint(devSteps);
        if (rc == -1) {
            printf("  solved after %llu iterations of %u walks in %.2f s (%llu distinguished points)\n",
                   (unsigned long long)(devSteps - steps0), nw, nowSeconds() - t0, (unsigned long long)totalDp);
            rc = 0;
        } else if (rc == 0) {
            const double el = nowSeconds() - t0;
            printf("  finished: %.3f M it/s, %llu distinguished points (%llu verified against the reference, %llu dropped)\n",
                   el > 0 ? (double)(devSteps - steps0) / el / 1e6 : 0.0, (unsigned long long)totalDp,
                   (unsigned long long)verified, (unsigned long long)lost);
        }
        fflush(stdout);
        return rc;
    }

    // The fast host arithmetic against the client's independent reference.
    int selftest() {
        int bad = 0;
        u64 rng = 0x2C13F2u;
        auto next = [&]() { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; };
        for (int t = 0; t < 8; ++t) {
            const u64 seed = seedOf(7, 1, t, t);
            P131 x, y;
            T.startPoint(seed, x, y);
            U192 alpha;
            const R::Point want = R::startPoint(seed, sol.basis, sol.target, &alpha, sol.ell, sol.spow);
            if (!same(x, pack(want.x)) || !same(y, pack(want.y))) { printf("FAIL: start point %d\n", t); ++bad; }
            if (!R::onCurve(want)) { printf("FAIL: start point %d off curve\n", t); ++bad; }
            // a few steps of the walk
            R::Point p = want;
            for (int s = 0; s < 6; ++s) {
                const int hw = R::weight(p.x);
                if (hw != weight(x)) { printf("FAIL: weight at step %d\n", s); ++bad; }
                p = R::step(p, hw);
                stepPacked(x, y, hw);
                if (!same(x, pack(p.x)) || !same(y, pack(p.y))) { printf("FAIL: step %d of seed %d\n", s, t); ++bad; break; }
            }
        }
        // rewalk at a loose cutoff against the reference walk
        for (int t = 0; t < 4; ++t) {
            const u64 seed = seedOf(9, 2, t, next() & 0xFF);
            const int w = 56;
            const FastWalk f = rewalkPacked(T, seed, w, 1 << 20);
            Solver<CfgF131> s2 = sol;
            s2.dpWeight = w;
            const Solver<CfgF131>::WalkResult r = s2.rewalk(seed);
            if (!f.ok || !r.ok || f.iters != r.iters || !same(f.x, pack(r.endPoint.x)) ||
                memcmp(f.counts, r.counts, sizeof f.counts) != 0) {
                printf("FAIL: rewalk of seed %016llx\n", (unsigned long long)seed);
                ++bad;
            }
        }
        // seed successor never leaves the id
        for (int t = 0; t < 3; ++t) {
            const u64 s = seedOf(3, 0xFFFF, 5, 0xFFFE);
            const u64 s1 = nextSeed(s), s2 = nextSeed(s1);
            if (((s1 >> 16) & 0xFFFF) != 5 || ((s2 >> 16) & 0xFFFF) != 5 || (s2 & 0xFFFF) != 0 ||
                ((s2 >> 32) & 0xFFFF) != 0 || ((s2 >> 48) != 3)) {
                printf("FAIL: seed successor\n");
                ++bad;
                break;
            }
        }
        printf(bad ? "selftest: FAIL (%d)\n" : "selftest: PASS\n", bad);
        return bad ? 1 : 0;
    }
};

static void usage() {
    printf(
        "ecc2k130-fpga - host program for the ECC2K-130 FPGA rho engine\n"
        "  --curve 131          the only curve the engine implements\n"
        "  --device D           FPGA slot (ECC_GPU from worker.py); default 0\n"
        "  --sim                model the device in software\n"
        "                       (--sim-engines N --sim-idw W --sim-steps S --sim-dp-weight W)\n"
        "  --dp-weight W        must equal the bitstream's (default: the challenge's 34)\n"
        "  --run-id R           16-bit salt making seeds unique across processes\n"
        "  --verify N           recompute the first N reported points on the host\n"
        "  --dp-file F          append distinguished points to F (binary, 32 bytes each)\n"
        "  --load F             preload a corpus file so collisions with earlier runs count\n"
        "  --load-max N         stop reloading after N points (0 = no limit), newest file first\n"
        "  --checkpoint F       save and resume seeds through F\n"
        "  --checkpoint-every S seconds between checkpoints (default 300)\n"
        "  --launches L         stop after L polls of the device (0 = until solved)\n"
        "  --dps N              stop after N distinguished points\n"
        "  --poll-ms T          sleep between polls of a real device (default 5)\n"
        "  --selftest           check the host arithmetic against the reference and exit\n"
        "  --threads --steps --packed --dp-cap --max-iters   accepted for worker.py, ignored\n");
}

int main(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const bool nx = i + 1 < argc;
        if (a == "--curve" && nx) o.curve = atoi(argv[++i]);
        else if (a == "--device" && nx) o.device = atoi(argv[++i]);
        else if (a == "--dp-weight" && nx) o.dpWeight = atoi(argv[++i]);
        else if (a == "--run-id" && nx) o.runId = (unsigned)atoi(argv[++i]);
        else if (a == "--verify" && nx) o.verify = atoi(argv[++i]);
        else if (a == "--dp-file" && nx) o.dpFile = argv[++i];
        else if (a == "--load" && nx) o.loadFiles.push_back(argv[++i]);
        else if (a == "--load-max" && nx) o.loadMax = strtoull(argv[++i], 0, 10);
        else if (a == "--checkpoint" && nx) o.ckptFile = argv[++i];
        else if (a == "--checkpoint-every" && nx) o.ckptSeconds = atof(argv[++i]);
        else if (a == "--launches" && nx) o.launches = atol(argv[++i]);
        else if (a == "--dps" && nx) o.dps = strtoull(argv[++i], 0, 10);
        else if (a == "--max-iters" && nx) o.maxIters = strtoull(argv[++i], 0, 10);
        else if (a == "--poll-ms" && nx) o.pollMs = atof(argv[++i]);
        else if (a == "--sim") o.sim = true;
        else if (a == "--sim-engines" && nx) o.simEng = atoi(argv[++i]);
        else if (a == "--sim-idw" && nx) o.simIdW = atoi(argv[++i]);
        else if (a == "--sim-steps" && nx) o.simSteps = (unsigned)atoi(argv[++i]);
        else if (a == "--sim-dp-weight" && nx) o.simDpWeight = atoi(argv[++i]);
        else if (a == "--selftest") o.selftest = true;
        else if ((a == "--threads" || a == "--steps" || a == "--dp-cap") && nx) ++i;
        else if (a == "--packed") {}
        else if (a == "--help" || a == "-h") { usage(); return 0; }
        else { printf("unknown option %s\n", a.c_str()); usage(); return 1; }
    }
    if (o.curve != 131) {
        fprintf(stderr, "the FPGA engine implements only --curve 131\n");
        return 1;
    }
    if (!o.sim) {
        const char *env = getenv("ECC_GPU");
        if (env && *env && o.device == 0) o.device = atoi(env);
    }
    Client c;
    c.o = o;
    const int rc = c.run();
    delete c.bus;
    return rc;
}
