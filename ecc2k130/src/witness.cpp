// Orbit witnesses for a corpus: the claim artifact a cairn piecework node pays.
//
// The client records only (seed, canonical endpoint), which is all a collision
// needs and is deliberately not checkable: reproducing one costs the steps
// that produced it -- 2^28.41 at the campaign's weight 32, 2^25.27 at the
// weight 34 a cairn job allows -- and a canonical name is free to invent --
// any low-weight bit string rotated to its least rotation is a syntactically
// perfect orbit name.  So a corpus record is worth nothing to anybody who did
// not walk it.
//
// A *witness* makes it worth something, and it comes out of the walk's own
// algebra.  One step is R -> R + sigma^j(R) = [1 + s^j]R and the endomorphism
// ring is commutative, so a trail of any length is [mu]R_0 with
//
//     mu = prod_j (1 + s^j)^{n_j},    n_j = steps that took branch j
//
// in any order.  Eight counters.  A claim carrying them is checked by one
// double scalar multiplication -- about 227 group operations against the
// 4e7 the trail cost -- which is the asymmetry that turns a trail into a
// payable artifact.
//
// There are two ways to get the n_j, and this tool reads both.
//
// A **v2 corpus** carries them: the walk counted its own branches on the
// device and wrote them beside the point, so emitting a claim costs one
// double scalar multiplication -- the same one the payer will do -- and the
// binding check is that mu lands on the orbit the record names.
//
// A **v1 corpus** does not, so the trail has to be replayed on the CPU to
// recover them.  Solver::rewalk already counts the n_j, because collision
// resolution has always needed them, but a replay costs what the trail cost:
// 2^25.27 steps for one ECC2K-130 record.  v1 is supported because every
// corpus written before the counters existed is one, not because it is a
// reasonable way to emit at scale.
//
// Either way the witness is verified against the record it claims before
// anything is printed.
//
// Two coordinate systems meet here and they are not the same one.  This
// client works in the *permuted* type-II ONB, where sigma is the coordinate
// permutation i -> fold(2i) on 1..m, not a rotation.  A cairn job pins a
// plain normal basis by its `nb_generator`, where sigma is a rotation and the
// orbit's name is the least rotation of the coordinate string.  Both bases
// are the same set of conjugates, so the map between them is a permutation:
// the job's generator is this basis's element T, and cairn's coordinate k is
// this one's fold(T * 2^k).  --job finds T by looking the generator up in
// GAMMA_TO_PB, so a job that names a different normal element still works and
// one that names a non-basis element fails loudly instead of emitting names
// nobody can check.
//
// Usage:
//   build/witness --curve 131 --job ecc2k130.json --corpus dps.bin > claims.jsonl
//   build/witness --curve 23 --instance 0 --job small.json --corpus dps.bin
//
// One JSON object per line, each a complete batch artifact:
//
//   {"dps":[{"j":[n3,...,n10],"seed":"<hex>","x":"<canonical orbit, hex>"},...]}
//
// Build:  make witness            (plain C++, no GPU)
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <string>
#include <vector>

#include "../include/curveparams.h"
#include "../include/kernel.h"
#include "../include/solver.h"

struct DpFileRecord {
    unsigned long long seed;
    unsigned long long canon[3];
};
static_assert(sizeof(DpFileRecord) == 32, "corpus records must remain 32 bytes");

// Corpus v2 carries the witness the walk already computed, so emitting a claim
// costs one double scalar multiplication instead of replaying the trail.  A v1
// corpus has no counts and still has to be replayed; both are read here.
struct DpFileRecordV2 {
    unsigned long long seed;
    unsigned long long iters;
    unsigned long long canon[3];
    unsigned counts[8];
};

struct DpFileHeader {
    char magic[8];
    unsigned version;
    unsigned recordBytes;
};

static const char DP_MAGIC_V2[8] = {'E', 'C', 'C', '2', 'K', 'D', 'P', '2'};

// What the rest of this tool works with, whichever format produced it.
struct CorpusRecord {
    unsigned long long seed;
    unsigned long long canon[3];
    unsigned long long iters;
    unsigned long long counts[8];
    bool hasWitness;
};

struct Options {
    int curve = 131;
    int instance = -1;
    bool polyBasis = false;
    int dpWeight = -1;
    unsigned long long maxIters = ECC_REPLAY_MAX_ITERS;
    std::string corpus;
    std::string job;
    std::string nbGenerator;
    unsigned long long skip = 0;
    unsigned long long max = 0;      // 0 = all
    int batch = 64;
    std::string outDir;
    bool quiet = false;
};

static void usage() {
    fprintf(stderr,
            "usage: witness --corpus F [--job J | --nb-generator HEX] [options]\n"
            "  --curve C          23, 41, 83 or 131 (default 131)\n"
            "  --instance I       planted test instance on a small curve\n"
            "  --dp-weight W      distinguishing weight (default: the curve's)\n"
            "  --max-iters N      refuse a trail longer than N steps (default %llu:\n"
            "                     the campaign guard plus its overshoot; never lower)\n"
            "  --job J            cairn job document; pins the normal basis and\n"
            "                     is checked against this binary's constants\n"
            "  --nb-generator H   the job's normal element as polynomial-basis\n"
            "                     hex, when no job file is at hand\n"
            "  --skip N/--max N   a slice of the corpus\n"
            "  --batch N          points per artifact (default 64, cairn's cap)\n"
            "  --out-dir D        one artifact per file, batch-00000.json ... , which\n"
            "                     is what a node's checker reads; default is a stream\n"
            "                     of one batch per line on stdout\n"
            "  --quiet            suppress the progress line on stderr\n",
            (unsigned long long)ECC_REPLAY_MAX_ITERS);
}

// ---- a very small reader for the flat job document ------------------------
// Only the fields this tool has to agree with; anything else in the file is
// the node's business.  A field that is present and disagrees is fatal, and a
// field that is absent is simply not checked -- the point is to refuse a job
// this binary cannot actually walk, not to re-validate cairn's schema.
static bool jobField(const std::string &src, const char *key, std::string *out) {
    const std::string pat = std::string("\"") + key + "\"";
    size_t p = src.find(pat);
    if (p == std::string::npos) return false;
    p = src.find(':', p + pat.size());
    if (p == std::string::npos) return false;
    ++p;
    while (p < src.size() && (src[p] == ' ' || src[p] == '\t' || src[p] == '\n')) ++p;
    if (p < src.size() && src[p] == '"') {
        const size_t e = src.find('"', ++p);
        if (e == std::string::npos) return false;
        *out = src.substr(p, e - p);
        return true;
    }
    size_t e = p;
    while (e < src.size() && src[e] != ',' && src[e] != '}' && src[e] != '\n') ++e;
    *out = src.substr(p, e - p);
    while (!out->empty() && (out->back() == ' ' || out->back() == '\r')) out->pop_back();
    return true;
}

static bool jobInt(const std::string &src, const char *key, long long *out) {
    std::string s;
    if (!jobField(src, key, &s) || s.empty()) return false;
    *out = strtoll(s.c_str(), 0, 10);
    return true;
}

// Expect `key` to equal `want` when it is present.  Returns false on conflict.
static bool jobExpect(const std::string &src, const char *key, long long want, const char *what) {
    long long got;
    if (!jobInt(src, key, &got)) return true;
    if (got == want) return true;
    fprintf(stderr, "job %s is %lld, this binary walks %lld (%s)\n", key, got, want, what);
    return false;
}

static bool parseHex192(const std::string &hex, unsigned long long *v) {
    v[0] = v[1] = v[2] = 0;
    if (hex.empty() || hex.size() > 48) return false;
    for (char c : hex) {
        int d;
        if (c >= '0' && c <= '9') d = c - '0';
        else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
        else return false;
        // v <<= 4
        v[2] = (v[2] << 4) | (v[1] >> 60);
        v[1] = (v[1] << 4) | (v[0] >> 60);
        v[0] = (v[0] << 4) | (unsigned long long)d;
    }
    return true;
}

static bool readCorpus(const std::string &path, std::vector<CorpusRecord> *recs) {
    struct stat st;
    FILE *in = fopen(path.c_str(), "rb");
    if (!in || fstat(fileno(in), &st) != 0) {
        fprintf(stderr, "corpus %s is missing or unreadable\n", path.c_str());
        if (in) fclose(in);
        return false;
    }
    // v2 announces itself with a magic rather than with its size: a truncated
    // v1 file and a v2 file would otherwise each look like a valid file of the
    // other format and mis-frame every record after the first.
    char magic[8];
    const bool v2 = fread(magic, 1, sizeof magic, in) == sizeof magic &&
                    memcmp(magic, DP_MAGIC_V2, sizeof magic) == 0;
    const long base = v2 ? (long)sizeof(DpFileHeader) : 0;
    const long recBytes = v2 ? (long)sizeof(DpFileRecordV2) : (long)sizeof(DpFileRecord);
    if (st.st_size < base || (st.st_size - base) % recBytes != 0) {
        fprintf(stderr, "corpus %s is not a whole number of %ld-byte %s records\n",
                path.c_str(), recBytes, v2 ? "v2" : "v1");
        fclose(in);
        return false;
    }
    if (fseek(in, base, SEEK_SET) != 0) { fclose(in); return false; }
    recs->resize((size_t)((st.st_size - base) / recBytes));
    bool ok = true;
    for (size_t i = 0; i < recs->size() && ok; ++i) {
        CorpusRecord &out = (*recs)[i];
        if (v2) {
            DpFileRecordV2 fr;
            ok = fread(&fr, sizeof fr, 1, in) == 1;
            if (!ok) break;
            out.seed = fr.seed;
            out.iters = fr.iters;
            for (int k = 0; k < 3; ++k) out.canon[k] = fr.canon[k];
            for (int k = 0; k < 8; ++k) out.counts[k] = fr.counts[k];
            out.hasWitness = true;
        } else {
            DpFileRecord fr;
            ok = fread(&fr, sizeof fr, 1, in) == 1;
            if (!ok) break;
            out.seed = fr.seed;
            out.iters = 0;
            for (int k = 0; k < 3; ++k) out.canon[k] = fr.canon[k];
            for (int k = 0; k < 8; ++k) out.counts[k] = 0;
            out.hasWitness = false;
        }
    }
    fclose(in);
    if (!ok) fprintf(stderr, "short read on %s\n", path.c_str());
    return ok;
}

// ---- the coordinate map ---------------------------------------------------
// fold(t) is the index of sigma^k applied to basis element t: the type-II ONB
// is indexed by 1..m standing for zeta^i + zeta^-i, so an index past m folds
// back through the ring of size 2m+1.
static int foldIndex(long long t, int nring, int m) {
    t %= nring;
    if (t < 0) t += nring;
    return (int)(t <= m ? t : nring - t);
}

template <class Cfg>
static int run(const Options &o, const unsigned long long *px, const unsigned long long *py,
               const unsigned long long *qx, const unsigned long long *qy, const char *ellDec,
               const char *sDec, int defaultW, const unsigned long long (*gammaToPb)[3]) {
    typedef Ref<Cfg> R;
    static const int M = Cfg::M;
    static const int NRING = Cfg::NRING;

    const int w = o.dpWeight < 0 ? defaultW : o.dpWeight;
    unsigned long long maxIters = o.maxIters;

    // The job, if given: agreement with this binary, then the normal element.
    std::string gammaHex = o.nbGenerator;
    if (!o.job.empty()) {
        FILE *jf = fopen(o.job.c_str(), "rb");
        if (!jf) { fprintf(stderr, "cannot read job %s\n", o.job.c_str()); return 4; }
        std::string src;
        char buf[4096];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, jf)) > 0) src.append(buf, n);
        fclose(jf);
        if (!jobExpect(src, "m", M, "field degree") ||
            !jobExpect(src, "dp_max_weight", w, "distinguishing weight") ||
            !jobExpect(src, "j_base", 3, "first branch") ||
            !jobExpect(src, "j_count", 8, "branch count") ||
            !jobExpect(src, "start_terms", 128, "start-point terms")) return 4;
        long long cap;
        if (jobInt(src, "max_steps_per_walker", &cap) && cap > 0 &&
            (unsigned long long)cap < maxIters) {
            fprintf(stderr, "note: job caps a trail at %lld steps; using that\n", cap);
            maxIters = (unsigned long long)cap;
        }
        std::string wit;
        if (jobField(src, "witness", &wit) && wit != "j-counts") {
            fprintf(stderr, "job asks for witness \"%s\"; this tool emits j-counts\n", wit.c_str());
            return 4;
        }
        std::string g;
        if (gammaHex.empty() && jobField(src, "nb_generator", &g)) gammaHex = g;
        long long mb;
        if (jobInt(src, "max_batch", &mb) && mb > 0 && o.batch > mb) {
            fprintf(stderr, "job caps a batch at %lld points, --batch is %d\n", mb, o.batch);
            return 4;
        }
    }
    if (gammaHex.empty()) {
        fprintf(stderr, "need --job or --nb-generator: the orbit name is basis-dependent\n");
        return 4;
    }

    // Find the job's normal element among this basis's conjugates.  T is
    // 1-based, matching the ONB index convention (element i is zeta^i+zeta^-i).
    unsigned long long gamma[3];
    if (!parseHex192(gammaHex, gamma)) {
        fprintf(stderr, "nb_generator is not hex\n");
        return 4;
    }
    int T = 0;
    for (int i = 0; i < M; ++i) {
        if (gammaToPb[i][0] == gamma[0] && gammaToPb[i][1] == gamma[1] &&
            gammaToPb[i][2] == gamma[2]) { T = i + 1; break; }
    }
    if (!T) {
        fprintf(stderr,
                "the job's nb_generator is not a conjugate of this basis: no permutation "
                "between them, so no orbit name this tool emits would check\n");
        return 4;
    }
    // cairn coordinate k reads this basis's coordinate perm[k] (0-based).
    std::vector<int> perm(M);
    {
        std::vector<char> seen(M, 0);
        long long e = 1;
        for (int k = 0; k < M; ++k) {
            const int idx = foldIndex((long long)T * e, NRING, M) - 1;
            if (idx < 0 || idx >= M || seen[idx]) {
                fprintf(stderr, "basis map is not a permutation; refusing to emit\n");
                return 4;
            }
            seen[idx] = 1;
            perm[k] = idx;
            e = (2 * e) % NRING;
        }
    }

    Solver<Cfg> sol;
    sol.setup(px, py, qx, qy, ellDec, sDec, w, maxIters);
    std::string why;
    if (!sol.checkSetup(&why)) {
        fprintf(stderr, "parameter check failed: %s\n", why.c_str());
        return 5;
    }

    std::vector<CorpusRecord> recs;
    if (!readCorpus(o.corpus, &recs)) return 8;

    // Settle --out-dir before the first walk: a batch that turns out to be
    // unwritable at flush time is a batch of replays thrown away.
    if (!o.outDir.empty()) {
        struct stat st;
        if ((mkdir(o.outDir.c_str(), 0755) != 0 && errno != EEXIST) ||
            stat(o.outDir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) {
            fprintf(stderr, "cannot use --out-dir %s\n", o.outDir.c_str());
            return 9;
        }
    }

    // ---- the canonical orbit name, in the job's basis ----------------------
    // Read this basis's coordinates through perm to get the job's coordinate
    // string, then take its least rotation, which is what sigma does there.
    auto canonicalName = [&](const typename R::Elem &x, unsigned long long *out) {
        unsigned long long c[3] = {0, 0, 0};
        for (int k = 0; k < M; ++k) {
            const int i = perm[k];
            if ((x.v[i >> 6] >> (i & 63)) & 1ull) c[k >> 6] |= 1ull << (k & 63);
        }
        auto less = [](const unsigned long long *a, const unsigned long long *b) {
            for (int i = 2; i >= 0; --i) {
                if (a[i] != b[i]) return a[i] < b[i];
            }
            return false;
        };
        unsigned long long best[3] = {c[0], c[1], c[2]}, cur[3] = {c[0], c[1], c[2]};
        for (int r = 1; r < M; ++r) {
            // rotate left by one inside M bits
            const int top = M - 1;
            const unsigned long long carry = (cur[top >> 6] >> (top & 63)) & 1ull;
            unsigned long long nx[3];
            nx[0] = cur[0] << 1;
            nx[1] = (cur[1] << 1) | (cur[0] >> 63);
            nx[2] = (cur[2] << 1) | (cur[1] >> 63);
            nx[0] |= carry;
            // mask to M bits
            for (int i = 0; i < 3; ++i) {
                const int lo = i * 64;
                if (lo >= M) nx[i] = 0;
                else if (lo + 64 > M) nx[i] &= (M - lo >= 64) ? ~0ull : ((1ull << (M - lo)) - 1);
            }
            cur[0] = nx[0]; cur[1] = nx[1]; cur[2] = nx[2];
            if (less(cur, best)) { best[0] = cur[0]; best[1] = cur[1]; best[2] = cur[2]; }
        }
        out[0] = best[0]; out[1] = best[1]; out[2] = best[2];
    };

    auto hexOf = [](const unsigned long long *v, std::string *out) {
        char buf[64];
        snprintf(buf, sizeof buf, "%llx%016llx%016llx", v[2], v[1], v[0]);
        int i = 0;
        while (buf[i] == '0' && buf[i + 1]) ++i;   // no leading zero; "0" survives
        *out = buf + i;
    };

    // Clamp so a --skip past the corpus is an empty slice, as it was for the
    // serial loop; unclamped, limit - first below would wrap.
    const unsigned long long first = o.skip < recs.size() ? o.skip : recs.size();
    unsigned long long limit = recs.size();
    if (o.max && first + o.max < limit) limit = first + o.max;

    std::vector<std::string> batch;
    unsigned long long emitted = 0, steps = 0, carried = 0;
    unsigned long long batchNo = 0;
    // False when the batch could not be written; the caller stops there rather
    // than walk more records whose witnesses would go the same way.
    auto flush = [&]() -> bool {
        if (batch.empty()) return true;
        // A node's checker reads one artifact per file, so --out-dir is what a
        // submitter wants; the stream on stdout is for looking at.
        FILE *out = stdout;
        std::string path;
        if (!o.outDir.empty()) {
            char name[64];
            snprintf(name, sizeof name, "/batch-%05llu.json", batchNo);
            path = o.outDir + name;
            out = fopen(path.c_str(), "wb");
            if (!out) {
                fprintf(stderr, "cannot write %s\n", path.c_str());
                return false;
            }
        }
        fputs("{\"dps\":[", out);
        for (size_t i = 0; i < batch.size(); ++i) {
            if (i) fputc(',', out);
            fputs(batch[i].c_str(), out);
        }
        fputs("]}\n", out);
        if (out != stdout) fclose(out);
        ++batchNo;
        emitted += batch.size();
        batch.clear();
        return true;
    };

#if ECC_WALK_TABLE
    fprintf(stderr, "this binary is built with ECC_WALK_TABLE: the replay tracks the "
                    "table's coefficients directly and leaves the j-counts empty, so it "
                    "cannot produce a j-counts witness\n");
    return 7;
#else
    // The replay is the expensive half of this tool and every record is
    // independent of every other, so it runs in parallel while the emission
    // stays serial.  The output does not change: results are written back by
    // index and the batching below walks them in corpus order, so the bytes
    // are the same whatever the thread count -- which is what lets the tests
    // compare a threaded run against a corpus walked one record at a time.
    struct Replayed {
        std::string el;              // the artifact element, ready to emit
        unsigned long long iters;
        bool carried;
        int err;                     // 0, or the exit code this record earns
        std::string msg;             // what to print for that error
    };
    const size_t span = (size_t)(limit - first);
    std::vector<Replayed> done_(span);
    unsigned long long progress = 0;

#pragma omp parallel for schedule(dynamic, 1)
    for (long long t = 0; t < (long long)span; ++t) {
        Replayed &slot = done_[(size_t)t];
        slot.err = 0;
        slot.iters = 0;
        slot.carried = false;
        char note[256];
        const CorpusRecord &rec = recs[(size_t)first + (size_t)t];
        // A v2 corpus already carries the counts, so the claim costs one
        // scalar multiplication rather than a replay of the trail -- the
        // asymmetry the witness exists for.  A v1 corpus has to be walked.
        const typename Solver<Cfg>::WalkResult wr =
            rec.hasWitness ? sol.fromCounts(rec.seed, rec.counts, rec.iters)
                           : sol.rewalk(rec.seed);
        if (!wr.ok) {
            if (rec.hasWitness)
                snprintf(note, sizeof note,
                         "seed %016llx carries a witness that does not sum to its %llu "
                         "steps or does not reach a distinguished point", rec.seed, rec.iters);
            else
                snprintf(note, sizeof note,
                         "seed %016llx did not reach a distinguished point in %llu steps",
                         rec.seed, maxIters);
            slot.err = 6;
            slot.msg = note;
            continue;
        }
        // The replay must land on the orbit the record names, or the corpus
        // and this binary disagree about the walk and nothing below is worth
        // claiming.
        const typename R::Elem canon = R::canonical(wr.endPoint.x);
        if (canon.v[0] != rec.canon[0] || canon.v[1] != rec.canon[1] || canon.v[2] != rec.canon[2]) {
            snprintf(note, sizeof note, "seed %016llx %s a different orbit than its record names",
                     rec.seed, rec.hasWitness ? "carries a witness landing on" : "replays to");
            slot.err = 6;
            slot.msg = note;
            continue;
        }
        if (R::weight(wr.endPoint.x) > w) {
            snprintf(note, sizeof note, "seed %016llx ends on weight %d, past the job's %d",
                     rec.seed, R::weight(wr.endPoint.x), w);
            slot.err = 6;
            slot.msg = note;
            continue;
        }
        // The witness, checked here rather than trusted.
        //
        // For a carried witness this comparison is a tautology -- fromCounts
        // built the endpoint out of mu -- and the check that binds is the
        // orbit comparison above, which is the same thing the payer checks.
        // For a replay it is independent: the endpoint came from stepping and
        // the counts came from counting those steps.
        if (!rec.hasWitness) {
            // [mu]R_0 rather than [mu*alpha0]P + [mu]Q.  R_0 IS [alpha0]P + Q,
            // so this is the payer's statement reached with one scalar
            // multiplication instead of two -- 196 point operations saved on
            // every replayed record, which is most of what a short trail costs.
            if (!R::eq(R::scalarMul(wr.startPt, sol.multiplier(wr.counts)), wr.endPoint)) {
                snprintf(note, sizeof note,
                         "seed %016llx: the j-counts do not reproduce the endpoint", rec.seed);
                slot.err = 6;
                slot.msg = note;
                continue;
            }
        }
        unsigned long long name[3];
        canonicalName(wr.endPoint.x, name);
        std::string xhex, shex;
        hexOf(name, &xhex);
        {
            const unsigned long long s[3] = {rec.seed, 0, 0};
            hexOf(s, &shex);
        }
        std::string el = "{\"j\":[";
        for (int u = 0; u < 8; ++u) {
            char b[32];
            snprintf(b, sizeof b, "%s%llu", u ? "," : "", wr.counts[u]);
            el += b;
        }
        el += "],\"seed\":\"" + shex + "\",\"x\":\"" + xhex + "\"}";
        slot.el = el;
        slot.iters = wr.iters;
        slot.carried = rec.hasWitness;
        if (!o.quiet) {
            unsigned long long seen;
#pragma omp atomic capture
            seen = ++progress;
            if (seen % 64 == 0) fprintf(stderr, "\rwitnessed %llu of %llu", seen, (unsigned long long)span);
        }
    }

    // Serial from here: the first failure IN CORPUS ORDER is the one reported,
    // so a threaded run fails exactly where a single-threaded one would.
    for (size_t t = 0; t < span; ++t) {
        const Replayed &slot = done_[t];
        if (slot.err) {
            fprintf(stderr, "%s\n", slot.msg.c_str());
            return slot.err;
        }
        batch.push_back(slot.el);
        steps += slot.iters;
        if (slot.carried) ++carried;
        if ((int)batch.size() >= o.batch && !flush()) return 9;
    }
#endif
    if (!flush()) return 9;
    if (!o.quiet) {
        // Separate the two numbers: steps replayed is the cost this run paid,
        // steps carried is the cost the walk had already paid for it.
        fprintf(stderr, "\r%llu witnesses from %llu records, %llu steps carried, "
                        "%llu steps replayed\n",
                emitted, limit - first, carried ? steps : 0ull, carried ? 0ull : steps);
    }
    return 0;
}

int main(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const bool nx = i + 1 < argc;
        if (a == "--curve" && nx) o.curve = atoi(argv[++i]);
        else if (a == "--instance" && nx) o.instance = atoi(argv[++i]);
        else if (a == "--poly-basis") o.polyBasis = true;
        else if (a == "--dp-weight" && nx) o.dpWeight = atoi(argv[++i]);
        else if (a == "--max-iters" && nx) o.maxIters = strtoull(argv[++i], 0, 10);
        else if (a == "--corpus" && nx) o.corpus = argv[++i];
        else if (a == "--job" && nx) o.job = argv[++i];
        else if (a == "--nb-generator" && nx) o.nbGenerator = argv[++i];
        else if (a == "--skip" && nx) o.skip = strtoull(argv[++i], 0, 10);
        else if (a == "--max" && nx) o.max = strtoull(argv[++i], 0, 10);
        else if (a == "--batch" && nx) o.batch = atoi(argv[++i]);
        else if (a == "--out-dir" && nx) o.outDir = argv[++i];
        else if (a == "--quiet") o.quiet = true;
        else { usage(); return 1; }
    }
    if (o.corpus.empty() || o.batch < 1) { usage(); return 1; }

#define W_DISPATCH(NS, CFG)                                                                    \
    do {                                                                                       \
        const unsigned long long *px = NS::PX, *py = NS::PY, *qx = NS::QX, *qy = NS::QY;        \
        if (o.instance >= 0 && o.instance < NS::NUM_INSTANCES) {                               \
            px = NS::INSTANCE_PX[o.instance];                                                   \
            py = NS::INSTANCE_PY[o.instance];                                                   \
            qx = NS::INSTANCE_QX[o.instance];                                                   \
            qy = NS::INSTANCE_QY[o.instance];                                                   \
        } else if (o.instance >= 0) {                                                           \
            fprintf(stderr, "curve %d has %d planted instances\n", o.curve, NS::NUM_INSTANCES);  \
            return 1;                                                                           \
        }                                                                                       \
        return run<CFG>(o, px, py, qx, qy, NS::ELL_DEC, NS::S_DEC, NS::DP_WEIGHT,               \
                        NS::GAMMA_TO_PB);                                                       \
    } while (0)

    if (o.curve == 131) W_DISPATCH(eccF131, CfgF131);
    if (o.curve == 83) W_DISPATCH(eccF83, CfgF83);
    if (o.curve == 41 && !o.polyBasis) W_DISPATCH(eccF41, CfgF41);
    if (o.curve == 23) W_DISPATCH(eccF23, CfgF23);
    fprintf(stderr,
            "unsupported curve %d: an orbit witness needs the normal-basis walk "
            "(23, 41, 83 or 131)\n",
            o.curve);
    return 1;
}
