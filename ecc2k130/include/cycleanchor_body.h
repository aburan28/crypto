// Shared function body; deliberately included by separate host and device
// templates. No include guard: this is not a standalone translation unit.
// Keeping execution spaces explicit avoids compiling host-only Ref operations
// into an unused device specialization (or device intrinsics for the host).
    Point points[8];
    unsigned tags[8];
    Point cur = start;
    int length = 0;
    for (int i = 0; i < 8; ++i) {
        if (ops.distinguished(cur)) return raw; // never escape past a report
        points[i] = cur;
        tags[i] = i == 0 ? raw : ops.tag(cur);
        Point next;
        if (!ops.next(cur, tags[i], &next)) return raw;
        cur = next;
        if (ops.equal(cur, start)) { length = i + 1; break; }
    }
    if (!length) return raw;
    int anchor = -1;
    bool startEligible = false;
    for (int i = 0; i < length; ++i) {
        unsigned long long cyclic = ECC_HIST_EMPTY;
        for (int back = 4; back >= 1; --back)
            cyclic = eccHistPush(cyclic, tags[(i + 4 * length - back) % length]);
        if (!eccTagFruitless(tags[i], cyclic, m)) continue;
        if (i == 0) startEligible = true;
        if (anchor < 0 || ops.less(points[i], points[anchor])) anchor = i;
    }
    // Equivalent orbit keys may tie: covariance gives equivalent exit edges.
    if (!startEligible || anchor < 0 || ops.less(points[anchor], start)) return raw;
    return eccTag((eccTagH(raw) + 1) & (branches - 1), eccTagK(raw), eccTagEps(raw));
