// The iteration constant of the two walks, on the walks' own reference code.
//
//   walk-constant-host --n 23|41 --walk sigma|table [--dist native|uniform|ecc2k130]
//                      [--walks W] [--trials T] [--seed S] [--threads K]
//
// examples/ecc2k130_walk_constant.rs measures c, the factor by which a walk
// needs more iterations than a random mapping on the classes of <sigma, -1>,
// on an emulation: the table walk there orients T_h by a canonical
// representative instead of the device's phase and pivot bit, and avoids
// only the step back to the previous class instead of applying the device's
// tag-history rule.  This runs the same statistic on the device walks
// themselves -- Ref::step for the sigma walk, TableWalk::step for the table
// walk, both exactly as Solver::rewalk calls them, branch from the device's
// normal-basis weight -- so the emulation's substitutions can be checked
// rather than argued.  WALK-CONSTANT.md has the method and the numbers.
//
// `--dist` picks the branch as the emulation's option of the same name does:
// `native` is the device's (HW(x_n)/2) mod H and runs Ref::step /
// TableWalk::step unmodified; `uniform` hashes the class key with a salt
// drawn per trial; `ecc2k130` maps that hash onto the probabilities
// (HW/2) mod H has at n = 131, the distribution the campaign walks see.
// With a hashed branch the table walk still takes its phase, sign bit,
// cycle rule and addend from TableWalk -- only h is replaced -- and the
// sigma walk still adds sigma^(3+b)(R).  The table walk's H is the build's
// ECC_TABLE_BRANCHES, 8 in every campaign build.
//
// A trial: W walks from random starts [r]P, stepped in turn, until a step
// lands on a class (Ref::canonical of x) any of them has visited.  The
// count is the number of distinct classes visited by then, whose mean for a
// random mapping on N = (l - 1)/2m classes is 1 + Q(N).  The table walk
// draws a fresh table T_h = [r_h]P per trial, so trials average over
// mappings as well as starts; the sigma walk has one mapping per curve.
//
// Fruitless cycles: the table walk is additive, so steps whose addends cancel
// in pairs return a walk to a point it has left with no collision.  The rule
// removes the 2- and 4-step pairwise ones, but not longer ones, nor steps
// whose addends sum to zero through Frobenius itself: s^2 + s + 2 = 0 makes
// sigma^(k+2) T + sigma^(k+1) T + sigma^k T + sigma^k T = O in four steps with
// no pair.  Each table step compares the new point with the walk's last
// RECENT points; a formal return (pairwise or tau-relation, counted apart) is
// not a collision, and the walk restarts.  Per-step rates are printed beside
// their leading-order predictions.
//
// One JSON line on stdout, the human line on stderr.
#ifndef ECC_NO_CUDA
#define ECC_NO_CUDA
#endif
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>
#include "../include/curveparams.h"
#include "../include/kernel.h"
#include "../include/solver.h"

#if !ECC_WALK_TABLE
#error "build with -DECC_WALK_TABLE=1"
#endif

// How many of a walk's own recent points a table step is compared with.
static const int RECENT = 16;

struct Tally {
    unsigned long long trials = 0, restarts = 0, steps = 0, ruleFired = 0, shortReturns = 0;
    double visited = 0, visitedSq = 0;
    std::vector<unsigned long long> branchCounts = std::vector<unsigned long long>(16, 0);
    // fruitless[L]: returns to the point L steps back whose tags cancel in
    // pairs; relation[L]: formal returns that are not pairwise (tau-relations).
    std::vector<unsigned long long> fruitless = std::vector<unsigned long long>(RECENT + 1, 0);
    std::vector<unsigned long long> relation = std::vector<unsigned long long>(RECENT + 1, 0);
    void add(const Tally &t) {
        trials += t.trials;
        restarts += t.restarts;
        steps += t.steps;
        ruleFired += t.ruleFired;
        shortReturns += t.shortReturns;
        for (int i = 0; i <= RECENT; ++i) fruitless[i] += t.fruitless[i];
        for (int i = 0; i <= RECENT; ++i) relation[i] += t.relation[i];
        visited += t.visited;
        visitedSq += t.visitedSq;
        for (size_t i = 0; i < branchCounts.size(); ++i) branchCounts[i] += t.branchCounts[i];
    }
};

static unsigned long long splitmix(unsigned long long z) {
    z += 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

// P[(HW/2) mod h] for HW binomial on n coordinates, conditioned even.
static std::vector<double> weightBranchProbabilities(int n, int h) {
    std::vector<double> p(h, 0.0);
    double total = 0;
    for (int k = 0; k <= n; k += 2) {
        const double w = std::exp(std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0));
        p[(k / 2) % h] += w;
        total += w;
    }
    for (double &x : p) x /= total;
    return p;
}

enum Dist { NATIVE, UNIFORM, ECC2K130 };

// Whether steps with these tags return a walk to where it started for every
// table: per branch h, the signed sum of tau^k over its steps is 0 in Z[tau],
// tau^2 + tau + 2 = 0 (Frobenius here: s^2 + s + 2 = 0).  Exponents start
// after the largest cyclic gap in k mod m, exact for relations spanning < m.
static bool returnsFormally(const std::vector<unsigned> &tags, int m) {
    for (unsigned h0 : tags) {
        const int h = eccTagH(h0);
        std::vector<int> ks;
        for (unsigned t : tags) if (eccTagH(t) == h) ks.push_back(eccTagK(t));
        std::sort(ks.begin(), ks.end());
        ks.erase(std::unique(ks.begin(), ks.end()), ks.end());
        int base = ks[0], gap = -1;
        for (size_t i = 0; i < ks.size(); ++i) {
            const int next = ks[(i + 1) % ks.size()], g = (next - ks[i] - 1 + m) % m;
            if (g > gap) { gap = g; base = next; }
        }
        long long a = 0, b = 0;
        for (unsigned t : tags) {
            if (eccTagH(t) != h) continue;
            const int d = (eccTagK(t) - base + m) % m;
            long long x = 1, y = 0;   // tau^d = x + y tau
            for (int j = 0; j < d; ++j) { const long long nx = -2 * y; y = x - y; x = nx; }
            const long long sign = eccTagEps(t) ? -1 : 1;
            a += sign * x;
            b += sign * y;
        }
        if (a || b) return false;
    }
    return true;
}

// Whether step tags cancel in pairs: every tag's opposite is there as often.
static bool cancels(std::vector<unsigned> tags) {
    while (!tags.empty()) {
        const unsigned t = tags.back();
        tags.pop_back();
        size_t i = 0;
        while (i < tags.size() && tags[i] != (t ^ ECC_TAG_EPS)) ++i;
        if (i == tags.size()) return false;
        tags[i] = tags.back();
        tags.pop_back();
    }
    return true;
}
static const char *distName[] = {"native", "uniform", "ecc2k130"};

// Exact 1 + Q(N): the mean number of distinct points a random mapping on N
// points visits before its first repeat.
static double randomMappingRho(double n) {
    double term = 1, sum = 1;
    for (double k = 1; term > 1e-15 && k < n; k += 1) {
        term *= 1 - k / n;
        sum += term;
    }
    return sum;
}

template <class Cfg>
static int run(bool table, Dist dist, int walks, unsigned long long trials,
               unsigned long long seed, int threads) {
    typedef Ref<Cfg> R;
    typedef typename R::Point Point;
    static_assert(Cfg::M <= 64, "the class key is one word");
    if (!TableWalk<Cfg>::applicable()) {
        std::fprintf(stderr, "n = %d has no type-II normal basis\n", Cfg::M);
        return 1;
    }
    Solver<Cfg> *sol = new Solver<Cfg>;
    const char *ell, *s;
    const unsigned long long *PX, *PY, *QX, *QY;
    if (Cfg::M == 23) {
        PX = eccF23::PX; PY = eccF23::PY; QX = eccF23::QX; QY = eccF23::QY;
        ell = eccF23::ELL_DEC; s = eccF23::S_DEC;
    } else {
        PX = eccF41::PX; PY = eccF41::PY; QX = eccF41::QX; QY = eccF41::QY;
        ell = eccF41::ELL_DEC; s = eccF41::S_DEC;
    }
    sol->setup(PX, PY, QX, QY, ell, s, 0, 0);
    std::string why;
    if (!sol->checkSetup(&why)) {
        std::fprintf(stderr, "setup: %s\n", why.c_str());
        return 1;
    }
    const unsigned long long ellU = sol->ell.v[0];
    const int H = table ? ECC_TABLE_BRANCHES : 8;
    // Cumulative n = 131 branch probabilities, scaled to 64 bits.
    std::vector<unsigned long long> cdf;
    {
        double acc = 0;
        for (double p : weightBranchProbabilities(131, H)) {
            acc += p;
            cdf.push_back(acc >= 1.0 ? ~0ull : (unsigned long long)std::ldexp(acc, 64));
        }
    }
    auto hashedBranch = [&](unsigned long long key, unsigned long long salt) {
        const unsigned long long u = splitmix(key ^ salt);
        if (dist == UNIFORM) return int(u % (unsigned long long)H);
        int b = 0;
        while (b < H - 1 && u >= cdf[b]) ++b;
        return b;
    };

    std::atomic<unsigned long long> next(0);
    std::mutex mu;
    Tally total;
    auto worker = [&](int tid) {
        std::mt19937_64 rng(seed * 1000003ull + (unsigned long long)tid);
        std::uniform_int_distribution<unsigned long long> scalar(1, ellU - 1);
        TableWalk<Cfg> *tw = new TableWalk<Cfg>(sol->walk);
        Tally t;
        std::unordered_set<unsigned long long> seen;
        std::vector<Point> p(walks);
        std::vector<unsigned long long> hist(walks), key(walks), len(walks);
        // Each walk's last RECENT points and the tags that left them.
        std::vector<Point> recentPt((size_t)walks * RECENT);
        std::vector<unsigned> recentTag((size_t)walks * RECENT);
        while (next.fetch_add(1) < trials) {
            if (table) {
                for (int h = 0; h < H; ++h) {
                    const Point th = R::scalarMul(sol->basis, u192_from(scalar(rng)));
                    for (int k = 0; k < Cfg::M; ++k) tw->table[h][k] = R::frob(th, k);
                }
            }
            const unsigned long long salt = rng();
            seen.clear();
            bool done = false, restart = false;
            for (int w = 0; w < walks && !done; ++w) {
                p[w] = R::scalarMul(sol->basis, u192_from(scalar(rng)));
                hist[w] = ECC_HIST_EMPTY;
                len[w] = 0;
                key[w] = R::canonical(p[w].x).v[0];
                if (!seen.insert(key[w]).second) done = true;
            }
            while (!done) {
                for (int w = 0; w < walks; ++w) {
                    const int hw = R::weight(p[w].x);
                    Point q;
                    unsigned tag = ECC_TAG_NONE;
                    if (table) {
                        unsigned raw = tw->rawTag(p[w], hw);
                        if (dist != NATIVE)
                            raw = eccTag(hashedBranch(key[w], salt), eccTagK(raw), eccTagEps(raw));
                        t.branchCounts[eccTagH(raw)]++;
                        tag = TableWalk<Cfg>::resolveTag(raw, hist[w]);
                        t.ruleFired += tag != raw;
                        if (dist == NATIVE) {
                            q = tw->step(p[w], hw, &hist[w], nullptr, nullptr, sol->ell, sol->spow);
                        } else {
                            hist[w] = eccHistPush(hist[w], tag);
                            q = R::addPtRaw(p[w], tw->addend(tag));
                        }
                    } else if (dist == NATIVE) {
                        t.branchCounts[R::jOf(hw) - 3]++;
                        q = R::step(p[w], hw);
                    } else {
                        const int b = hashedBranch(key[w], salt);
                        t.branchCounts[b]++;
                        q = R::addPt(p[w], R::frob(p[w], 3 + b));
                    }
                    t.steps++;
                    if (q.inf || !R::onCurve(q)) { restart = true; break; }
                    if (table) {
                        // A return to one of the walk's own recent points whose
                        // tags cancel is a fruitless cycle, not a collision: count
                        // it by length and restart the walk.
                        const size_t base = (size_t)w * RECENT;
                        recentPt[base + len[w] % RECENT] = p[w];
                        recentTag[base + len[w] % RECENT] = tag;
                        len[w]++;
                        int fruitlessLen = 0;
                        for (int l = 1; l <= RECENT && (unsigned long long)l <= len[w]; ++l) {
                            if (!R::eq(recentPt[base + (len[w] - l) % RECENT], q)) continue;
                            std::vector<unsigned> tags;
                            for (int i = 1; i <= l; ++i) tags.push_back(recentTag[base + (len[w] - i) % RECENT]);
                            if (cancels(tags)) { t.fruitless[l]++; fruitlessLen = l; }
                            else if (returnsFormally(tags, Cfg::M)) { t.relation[l]++; fruitlessLen = l; }
                            else t.shortReturns++;
                            break;
                        }
                        if (fruitlessLen) {
                            p[w] = R::scalarMul(sol->basis, u192_from(scalar(rng)));
                            hist[w] = ECC_HIST_EMPTY;
                            len[w] = 0;
                            key[w] = R::canonical(p[w].x).v[0];
                            if (!seen.insert(key[w]).second) { done = true; break; }
                            continue;
                        }
                    }
                    p[w] = q;
                    key[w] = R::canonical(q.x).v[0];
                    if (!seen.insert(key[w]).second) { done = true; break; }
                }
                if (restart) break;
            }
            if (restart) { t.restarts++; continue; }
            const double v = (double)seen.size();
            t.trials++;
            t.visited += v;
            t.visitedSq += v * v;
        }
        delete tw;
        std::lock_guard<std::mutex> g(mu);
        total.add(t);
    };
    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i) pool.emplace_back(worker, i);
    for (auto &th : pool) th.join();

    const double classes = (double)(ellU - 1) / (2.0 * Cfg::M);
    const double expected = randomMappingRho(classes);
    const double mean = total.visited / total.trials;
    const double var = total.visitedSq / total.trials - mean * mean;
    const double se = std::sqrt(var / total.trials);
    double s2 = 0;
    unsigned long long counted = 0;
    for (int b = 0; b < H; ++b) counted += total.branchCounts[b];
    for (int b = 0; b < H; ++b) s2 += std::pow((double)total.branchCounts[b] / counted, 2);
    // First-order models (WALK-CONSTANT.md §2): a walk injective on each
    // branch pays 1/sqrt(1 - sum p^2); the table walk, whose same-branch
    // class collisions are suppressed only when the frames agree, pays
    // 1/sqrt(1 - sum p^2 / 2m).
    const double c = mean / expected, cse = se / expected;
    const double injective = 1 / std::sqrt(1 - s2), classFrame = 1 / std::sqrt(1 - s2 / (2.0 * Cfg::M));
    const double model = table ? classFrame : injective;
    // Fruitless cycles the rule lets through, to leading order: 4 pairwise
    // 6-step patterns, three branches each drawn twice, (sum p^2 / 2m)^3; and
    // 24 4-step tau-relations, one branch drawn four times, sum p^4 / (2m)^3.
    double s4 = 0;
    for (int b = 0; b < H; ++b) s4 += std::pow((double)total.branchCounts[b] / counted, 4);
    const double twoM = 2.0 * Cfg::M;
    const double pairwisePredicted = 4 * std::pow(s2 / twoM, 3), relationPredicted = 24 * s4 / std::pow(twoM, 3);
    auto byLength = [](const std::vector<unsigned long long> &v, unsigned long long *sum) {
        std::string out;
        *sum = 0;
        for (int l = 1; l <= RECENT; ++l) {
            *sum += v[l];
            if (!v[l]) continue;
            if (!out.empty()) out += ",";
            out += "\"" + std::to_string(l) + "\":" + std::to_string(v[l]);
        }
        return out;
    };
    unsigned long long fruitlessTotal = 0, relationTotal = 0;
    const std::string fruitlessJson = byLength(total.fruitless, &fruitlessTotal);
    const std::string relationJson = byLength(total.relation, &relationTotal);
    const double fruitlessRate = (double)fruitlessTotal / (double)total.steps;
    const double relationRate = (double)relationTotal / (double)total.steps;
    const char *walk = table ? "table" : "sigma";
    std::fprintf(stderr,
                 "n = %d, device %s walk, H = %d, %s branches, W = %d: %llu trials, mean %.1f "
                 "classes visited (random mapping %.1f), c = %.4f +- %.4f; sum p^2 = %.5f, "
                 "model %.4f, c / model = %.4f; cycle rule %llu times in %llu steps; fruitless: "
                 "pairwise %.3e/step (predicted %.3e) by length {%s}, tau-relation %.3e/step "
                 "(predicted %.3e) by length {%s}; short returns %llu; restarts %llu\n",
                 Cfg::M, walk, H, distName[dist], walks, total.trials, mean, expected, c, cse, s2,
                 model, c / model, total.ruleFired, total.steps, fruitlessRate, pairwisePredicted,
                 fruitlessJson.c_str(), relationRate, relationPredicted, relationJson.c_str(),
                 total.shortReturns, total.restarts);
    std::printf("{\"n\":%d,\"walk\":\"device-%s\",\"branches\":%d,\"dist\":\"%s\",\"walks\":%d,"
                "\"trials\":%llu,\"seed\":%llu,\"classes\":%.1f,\"random_mapping_rho\":%.4f,"
                "\"mean_visited\":%.4f,\"se_visited\":%.4f,\"c\":%.5f,\"c_se\":%.5f,\"sum_p2\":%.5f,"
                "\"injective_model\":%.5f,\"class_frame_model\":%.5f,\"cycle_rule\":%llu,"
                "\"restarts\":%llu,\"steps\":%llu,\"sum_p4\":%.6f,\"fruitless\":{%s},"
                "\"fruitless_per_step\":%.6e,\"fruitless_predicted\":%.6e,\"relation\":{%s},"
                "\"relation_per_step\":%.6e,\"relation_predicted\":%.6e,\"short_returns\":%llu}\n",
                Cfg::M, walk, H, distName[dist], walks, total.trials, seed, classes, expected, mean,
                se, c, cse, s2, injective, classFrame, total.ruleFired, total.restarts, total.steps,
                s4, fruitlessJson.c_str(), fruitlessRate, pairwisePredicted, relationJson.c_str(),
                relationRate, relationPredicted, total.shortReturns);
    std::fflush(stdout);
    delete sol;
    return 0;
}

int main(int argc, char **argv) {
    int n = 23, walks = 8, threads = 4;
    unsigned long long trials = 1000, seed = 1;
    bool table = false;
    Dist dist = NATIVE;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "--n")) n = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--walk")) {
            if (!std::strcmp(argv[i + 1], "table")) table = true;
            else if (std::strcmp(argv[i + 1], "sigma")) { std::fprintf(stderr, "unknown walk %s\n", argv[i + 1]); return 2; }
        } else if (!std::strcmp(argv[i], "--dist")) {
            if (!std::strcmp(argv[i + 1], "native")) dist = NATIVE;
            else if (!std::strcmp(argv[i + 1], "uniform")) dist = UNIFORM;
            else if (!std::strcmp(argv[i + 1], "ecc2k130")) dist = ECC2K130;
            else { std::fprintf(stderr, "unknown dist %s\n", argv[i + 1]); return 2; }
        } else if (!std::strcmp(argv[i], "--walks")) walks = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--trials")) trials = std::strtoull(argv[i + 1], nullptr, 10);
        else if (!std::strcmp(argv[i], "--seed")) seed = std::strtoull(argv[i + 1], nullptr, 10);
        else if (!std::strcmp(argv[i], "--threads")) threads = std::atoi(argv[i + 1]);
        else { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (walks < 1 || threads < 1) { std::fprintf(stderr, "walks and threads must be positive\n"); return 2; }
    if (n == 23) return run<CfgF23>(table, dist, walks, trials, seed, threads);
    if (n == 41) return run<CfgF41>(table, dist, walks, trials, seed, threads);
    std::fprintf(stderr, "n must be 23 or 41\n");
    return 2;
}
