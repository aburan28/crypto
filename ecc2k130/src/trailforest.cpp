// The walks behind the client's distinguished points, one orbit at a time.
//
// The client records only (seed, canonical endpoint) per distinguished point,
// which is all a collision needs.  This tool walks seeds with the reference
// implementation -- the same start-point derivation, iteration function and
// distinguishing test the client verifies its own reports against under
// --verify -- and prints every orbit a walk passes through on its way to its
// distinguished point.  scripts/site/walk_forest.py reads that and draws the
// forest the trails form, for the campaign status page.
//
// Two modes, one output format:
//
//   --corpus F        replay the records in the client's --dp-file F.  A trail
//                     whose last orbit is not the one the record names is an
//                     error, not a drawing.
//   --generate        walk the client's own seed schedule for --run-id R,
//                     lanes 0 .. --walks-1, exactly as the client would, but
//                     without stopping at the first collision the way a
//                     search does.  On a curve small enough to draw, a search
//                     solves within its first launch, so its corpus holds only
//                     the handful of walks reported before that; this mode
//                     keeps every lane's first walk.  --check-corpus F then
//                     holds the trails to the client's real records: every
//                     record in F for a generated seed must end on the orbit
//                     the replay reached, and the count of records checked is
//                     printed in the header.  --corpus-out F writes the
//                     generated endpoints in the client's 32-byte record
//                     format, so the client can --load them.
//
// Build:  make trailforest            (plain C++, no GPU)
// Usage:  build/trailforest --curve 23 --instance 0 --corpus dps.bin [--max N] > trails.txt
//         build/trailforest --curve 23 --instance 0 --generate --run-id 1 --walks 40 \
//                           --check-corpus dps.bin --corpus-out forest.bin > trails.txt
//
// Output, one record per line:
//
//   # curve 23 instance 0 dp-weight 7 mode generate run-id 1 walks 40 checked 3
//   walk <seed hex> <steps> <orbit0> <orbit1> ... <orbitN>
//
// where each orbit is the canonical representative of the x-coordinate's
// orbit under Frobenius (the collision key the store uses), as hex, and the
// last orbit on the line is the distinguished one.  <steps> is the number of
// iterations, so a line carries <steps> + 1 orbits.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <string>
#include <vector>

#include "../include/curveparams.h"
#include "../include/solver.h"

struct DpFileRecord {
    unsigned long long seed;
    unsigned long long canon[3];
};
static_assert(sizeof(DpFileRecord) == 32, "corpus records must remain 32 bytes");

struct Options {
    int curve = 23;
    int instance = -1;
    int dpWeight = -1;
    bool polyBasis = false;
    bool generate = false;
    unsigned runId = 1;
    unsigned long long walks = 40;
    unsigned long long max = 0;
    unsigned long long skip = 0;
    unsigned long long maxIters = 1ull << 32;
    std::string corpus;
    std::string checkCorpus;
    std::string corpusOut;
};

// The client's seed for lane `walkIndex` of run `runId` (walk.h, eccSeedFor):
// a finished lane continues at seed + 1, so the low 16 bits count restarts.
static unsigned long long seedFor(unsigned runId, unsigned long long walkIndex) {
    return ((unsigned long long)runId << 48) | ((walkIndex & 0xFFFFFFFFull) << 16);
}

static bool readCorpus(const std::string &path, std::vector<DpFileRecord> *recs) {
    FILE *in = fopen(path.c_str(), "rb");
    struct stat st;
    if (!in || fstat(fileno(in), &st) != 0 || st.st_size % sizeof(DpFileRecord) != 0) {
        if (in) fclose(in);
        fprintf(stderr, "cannot read corpus %s (missing, or not a whole number of 32-byte records)\n",
                path.c_str());
        return false;
    }
    recs->resize((size_t)(st.st_size / sizeof(DpFileRecord)));
    const bool ok = recs->empty() || fread(recs->data(), sizeof(DpFileRecord), recs->size(), in) == recs->size();
    fclose(in);
    if (!ok) fprintf(stderr, "short read on %s\n", path.c_str());
    return ok;
}

static void printHex(const unsigned long long *v) {
    // 192 bits, most significant limb first, leading zeros dropped past the
    // first digit so a small field does not print as 48 zeros.
    char buf[64];
    int n = snprintf(buf, sizeof buf, "%llx%016llx%016llx", v[2], v[1], v[0]);
    int i = 0;
    while (i < n - 1 && buf[i] == '0') ++i;
    fputs(buf + i, stdout);
}

template <class Cfg>
static int run(const Options &o, const unsigned long long *px, const unsigned long long *py,
               const unsigned long long *qx, const unsigned long long *qy, const char *ellDec,
               const char *sDec, int defaultW) {
    typedef Ref<Cfg> R;
    const int w = o.dpWeight < 0 ? defaultW : o.dpWeight;
    Solver<Cfg> sol;
    sol.setup(px, py, qx, qy, ellDec, sDec, w, o.maxIters);
    std::string why;
    if (!sol.checkSetup(&why)) {
        fprintf(stderr, "parameter check failed: %s\n", why.c_str());
        return 5;
    }

    // Walk one seed to its distinguished point, collecting every orbit.
    auto walk = [&](unsigned long long seed, std::vector<typename R::Elem> *orbits) -> bool {
        U192 alpha;
        typename R::Point p = R::startPoint(seed, sol.basis, sol.target, &alpha, sol.ell, sol.spow);
        orbits->clear();
        for (unsigned long long it = 0;; ++it) {
            orbits->push_back(R::canonical(p.x));
            const int hw = R::weight(p.x);
            if (hw <= w) return true;
            if (it >= o.maxIters) {
                fprintf(stderr, "seed %016llx did not reach a distinguished point within %llu steps\n",
                        seed, o.maxIters);
                return false;
            }
            p = R::step(p, hw);
        }
    };
    auto emit = [](unsigned long long seed, const std::vector<typename R::Elem> &orbits) {
        printf("walk %016llx %zu", seed, orbits.size() - 1);
        for (size_t k = 0; k < orbits.size(); ++k) {
            putchar(' ');
            printHex(orbits[k].v);
        }
        putchar('\n');
    };
    auto sameKey = [](const typename R::Elem &e, const DpFileRecord &rec) {
        return e.v[0] == rec.canon[0] && e.v[1] == rec.canon[1] && e.v[2] == rec.canon[2];
    };

    std::vector<typename R::Elem> orbits;

    if (o.generate) {
        // The header comes first, so the record count it reports has to be
        // known before any trail is printed; the walks are cheap enough to
        // keep.
        std::vector<unsigned long long> seeds;
        std::vector<std::vector<typename R::Elem> > trails;
        for (unsigned long long i = 0; i < o.walks; ++i) {
            const unsigned long long seed = seedFor(o.runId, i);
            if (!walk(seed, &orbits)) return 3;
            seeds.push_back(seed);
            trails.push_back(orbits);
        }
        size_t checked = 0;
        if (!o.checkCorpus.empty()) {
            std::vector<DpFileRecord> recs;
            if (!readCorpus(o.checkCorpus, &recs)) return 8;
            for (size_t r = 0; r < recs.size(); ++r) {
                for (size_t i = 0; i < seeds.size(); ++i) {
                    if (recs[r].seed != seeds[i]) continue;
                    if (!sameKey(trails[i].back(), recs[r])) {
                        fprintf(stderr, "seed %016llx: the client's record ends on an orbit the replay did not reach\n",
                                recs[r].seed);
                        return 3;
                    }
                    ++checked;
                }
            }
            if (!checked) {
                fprintf(stderr, "no record in %s carries a generated seed; nothing was checked\n",
                        o.checkCorpus.c_str());
                return 3;
            }
        }
        if (!o.corpusOut.empty()) {
            FILE *out = fopen(o.corpusOut.c_str(), "wb");
            if (!out) { fprintf(stderr, "cannot write %s\n", o.corpusOut.c_str()); return 8; }
            for (size_t i = 0; i < seeds.size(); ++i) {
                DpFileRecord fr;
                fr.seed = seeds[i];
                fr.canon[0] = trails[i].back().v[0];
                fr.canon[1] = trails[i].back().v[1];
                fr.canon[2] = trails[i].back().v[2];
                if (fwrite(&fr, sizeof fr, 1, out) != 1) { fprintf(stderr, "short write on %s\n", o.corpusOut.c_str()); return 8; }
            }
            fclose(out);
        }
        printf("# curve %d instance %d dp-weight %d mode generate run-id %u walks %llu checked %zu\n",
               o.curve, o.instance, w, o.runId, o.walks, checked);
        for (size_t i = 0; i < seeds.size(); ++i) emit(seeds[i], trails[i]);
        return 0;
    }

    std::vector<DpFileRecord> recs;
    if (!readCorpus(o.corpus, &recs)) return 8;
    const size_t first = (size_t)(o.skip < recs.size() ? o.skip : recs.size());
    size_t last = recs.size();
    if (o.max && first + o.max < last) last = first + (size_t)o.max;

    printf("# curve %d instance %d dp-weight %d mode replay records %zu skip %zu drawn %zu\n",
           o.curve, o.instance, w, recs.size(), first, last - first);
    for (size_t i = first; i < last; ++i) {
        const DpFileRecord &rec = recs[i];
        if (!walk(rec.seed, &orbits)) return 3;
        if (!sameKey(orbits.back(), rec)) {
            fprintf(stderr, "seed %016llx: the replayed walk ends at an orbit the corpus does not record\n",
                    rec.seed);
            return 3;
        }
        emit(rec.seed, orbits);
    }
    return 0;
}

static void usage() {
    fprintf(stderr,
            "trailforest - replay the walks behind a corpus file, orbit by orbit\n"
            "\n"
            "  --curve M        131, 83, 41, 23 (normal basis) or 97, 41 --poly-basis, 19, 13\n"
            "  --poly-basis     for curve 41, the polynomial-basis parameters\n"
            "  --instance I     planted test instance I (small curves only)\n"
            "  --dp-weight W    the cutoff the corpus was collected at (default: the curve's)\n"
            "  --corpus F       replay the client's --dp-file F (32-byte records)\n"
            "  --skip K         skip the first K records\n"
            "  --max N          replay at most N records (0 = all)\n"
            "  --generate       walk the client's seed schedule instead of a corpus\n"
            "  --run-id R       16-bit run id whose seeds to walk (default 1)\n"
            "  --walks N        lanes 0 .. N-1 of that run, first walk each (default 40)\n"
            "  --check-corpus F hold generated trails to the client's records in F\n"
            "  --corpus-out F   write the generated endpoints as a corpus file\n"
            "  --max-iters N    give up on a walk after N steps (default 2^32)\n");
}

int main(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const bool nx = i + 1 < argc;
        if (a == "--curve" && nx) o.curve = atoi(argv[++i]);
        else if (a == "--poly-basis") o.polyBasis = true;
        else if (a == "--instance" && nx) o.instance = atoi(argv[++i]);
        else if (a == "--dp-weight" && nx) o.dpWeight = atoi(argv[++i]);
        else if (a == "--corpus" && nx) o.corpus = argv[++i];
        else if (a == "--skip" && nx) o.skip = strtoull(argv[++i], 0, 10);
        else if (a == "--max" && nx) o.max = strtoull(argv[++i], 0, 10);
        else if (a == "--max-iters" && nx) o.maxIters = strtoull(argv[++i], 0, 10);
        else if (a == "--generate") o.generate = true;
        else if (a == "--run-id" && nx) o.runId = (unsigned)strtoul(argv[++i], 0, 10) & 0xFFFFu;
        else if (a == "--walks" && nx) o.walks = strtoull(argv[++i], 0, 10);
        else if (a == "--check-corpus" && nx) o.checkCorpus = argv[++i];
        else if (a == "--corpus-out" && nx) o.corpusOut = argv[++i];
        else { usage(); return 1; }
    }
    if (o.corpus.empty() == !o.generate) { usage(); return 1; }

#define TF_DISPATCH(NS, CFG)                                                                  \
    do {                                                                                      \
        const unsigned long long *px = NS::PX, *py = NS::PY, *qx = NS::QX, *qy = NS::QY;       \
        if (o.instance >= 0 && o.instance < NS::NUM_INSTANCES) {                              \
            px = NS::INSTANCE_PX[o.instance];                                                  \
            py = NS::INSTANCE_PY[o.instance];                                                  \
            qx = NS::INSTANCE_QX[o.instance];                                                  \
            qy = NS::INSTANCE_QY[o.instance];                                                  \
        } else if (o.instance >= 0) {                                                          \
            fprintf(stderr, "curve %d has %d planted instances\n", o.curve, NS::NUM_INSTANCES); \
            return 1;                                                                          \
        }                                                                                      \
        return run<CFG>(o, px, py, qx, qy, NS::ELL_DEC, NS::S_DEC, NS::DP_WEIGHT);            \
    } while (0)

    if (o.curve == 131) TF_DISPATCH(eccF131, CfgF131);
    if (o.curve == 83) TF_DISPATCH(eccF83, CfgF83);
    if (o.curve == 41 && !o.polyBasis) TF_DISPATCH(eccF41, CfgF41);
    if (o.curve == 23) TF_DISPATCH(eccF23, CfgF23);
    if (o.curve == 97) TF_DISPATCH(eccP97, CfgP97);
    if (o.curve == 41 && o.polyBasis) TF_DISPATCH(eccP41, CfgP41);
    if (o.curve == 19) TF_DISPATCH(eccP19, CfgP19);
    if (o.curve == 13) TF_DISPATCH(eccP13, CfgP13);
    fprintf(stderr, "unsupported curve %d\n", o.curve);
    return 1;
}
