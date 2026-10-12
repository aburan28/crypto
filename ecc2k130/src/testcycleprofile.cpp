// Exact v3 cold-probe accounting on synthetic graphs. No GPU and no ECDLP.
#define ECC_NO_CUDA
#define ECC_WALK_TABLE 1
#ifndef ECC_CYCLE_PROFILE
#define ECC_CYCLE_PROFILE 1
#endif
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/tablewalk.h"

static void require(bool ok, const char *why) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", why);
        std::exit(1);
    }
}

struct GraphOps {
    const unsigned *tags;
    int n;
    bool closes;
    int dp;
    int exceptional;

    bool distinguished(int p) const { return p == dp; }
    unsigned tag(int p) const { return tags[p % n]; }
    bool next(int p, unsigned t, int *q) const {
        require(t == tag(p), "profiled probe uses raw steps");
        if (p == exceptional) return false;
        *q = closes ? (p + 1) % n : p + 1;
        return true;
    }
    bool oppositeCloses(int start, int after, unsigned t) const {
        require(t == tag(after), "fast2 checks the next raw tag");
        return closes && (after + 1) % n == start;
    }
    bool equal(int a, int b) const { return a == b; }
    bool less(int a, int b) const { return a < b; }
};

static unsigned profile(const GraphOps &ops, EccCycleProfile *totals) {
    EccCycleProbeProfile probe = {};
    const unsigned raw = ops.tag(0);
    const unsigned result = eccCycleAnchorTagProfile(0, raw, ops, 131, 8, &probe);
    require(probe.outcome != ECC_CYCLE_OUTCOME_NONE, "every probe records a terminal outcome");
    eccCycleProfileAccumulate(totals, probe);
    return result;
}

static std::vector<unsigned> distinctTags(int n) {
    std::vector<unsigned> tags;
    for (int i = 0; i < n; ++i) tags.push_back(eccTag(i & 7, 11 * i, 0));
    return tags;
}

int main() {
    EccCycleProfile totals = {};

    // One exact general cycle of each length. Length four uses a tau relation
    // whose start is the least eligible vertex, so it also exercises an exit.
    for (int n = 1; n <= 8; ++n) {
        std::vector<unsigned> tags = distinctTags(n);
        if (n == 2) {
#if ECC_CYCLE_FAST2
            // Keep this in the general path; fast2 gets its own control below.
            tags = {eccTag(0, 5, 0), eccTag(1, 5, 0)};
#else
            tags = {eccTag(0, 5, 0), eccTag(0, 5, 1)};
#endif
        }
        if (n == 4) {
            const unsigned r = eccTag(0, 37, 0);
            tags = {r, r, eccTagAdvanceK(r, 1, 131), eccTagAdvanceK(r, 2, 131)};
        }
        const unsigned raw = tags[0];
        const unsigned result = profile(GraphOps{tags.data(), n, true, -1, -1}, &totals);
        if (n == 4 || (n == 2 && !ECC_CYCLE_FAST2))
            require(result != raw, "least eligible general-cycle anchor exits");
        else
            require(result == raw, "ineligible general cycle keeps the raw edge");
    }

#if ECC_CYCLE_FAST2
    {
        const unsigned raw = eccTag(0, 19, 0);
        const unsigned tags[2] = {raw, raw ^ ECC_TAG_EPS};
        require(profile(GraphOps{tags, 2, true, -1, -1}, &totals) != raw,
                "fast2 cycle exits at its least eligible anchor");
    }
#endif

    // The three non-cycle terminal classes, with exact affine-call budgets.
    std::vector<unsigned> openTags = distinctTags(9);
    profile(GraphOps{openTags.data(), 9, false, -1, -1}, &totals);
    profile(GraphOps{openTags.data(), 9, false, 3, -1}, &totals);
    profile(GraphOps{openTags.data(), 9, false, -1, 4}, &totals);

    const unsigned long long fast2 = ECC_CYCLE_FAST2 ? 1 : 0;
    const unsigned long long expectedHints = 11 + fast2;
    const unsigned long long expectedNextCalls = 52 + fast2;
    require(totals.hints == expectedHints, "total hints match the synthetic corpus");
    require(totals.fast2Hits == fast2, "fast2 hits match the build and corpus");
    for (int i = 0; i < 8; ++i)
        require(totals.generalCycles[i] == 1, "one general cycle of every length 1..8");
    require(totals.openEightStepProbes == 1, "one open eight-step probe");
    require(totals.dpAborts == 1, "one distinguished-point abort");
    require(totals.exceptionalDenominatorAborts == 1,
            "one exceptional-denominator abort");
    require(totals.affineNextCalls == expectedNextCalls,
            "affine next calls reconcile with all probe paths");
    require(totals.anchorExits == 2, "general and fast2/control anchor exits reconcile");
    require(eccCycleProfileTerminalTotal(totals) == totals.hints,
            "terminal outcomes partition total hints");

    std::printf("PASS: hints %llu; fast2 %llu; cycles 1,1,1,1,1,1,1,1; "
                "open 1; dp 1; exceptional 1; next %llu; exits %llu\n",
                totals.hints, totals.fast2Hits, totals.affineNextCalls,
                totals.anchorExits);
}
