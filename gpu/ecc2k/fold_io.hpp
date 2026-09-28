/* fold_io.hpp -- the files a folded build reads and writes.  Host only.
 *
 * A **plan** is what a build needs and cannot compute on its own: the
 * factor base sorted by signed Frobenius orbit, one representative per
 * orbit, where each orbit's points begin, and the normal basis the keys
 * are named in.  `examples/dump_fold_plan.rs` writes it from
 * `PairSumTable::folded_rows`, so a device is handed the CPU's plan
 * rather than a second derivation of it.
 *
 * A **table** is what a build produces, as `PairSumTable` stores it:
 * bucket offsets, tagged words and the presence filter.
 * `examples/load_fold_table.rs` loads it with
 * `PairSumTable::from_folded_parts` and checks it against the table the
 * CPU builds itself.
 *
 * Both carry the degree and the base's selection (`seed`, the point
 * count asked for), which is all the CPU needs to rebuild the base.
 * Little-endian, the layout an x86 host writes natively; the Rust side
 * reads it that way wherever it runs, and the writers here refuse to run
 * on a big-endian host rather than write something else.
 *
 * `test_pairtable_emu.cpp` uses both, which is what keeps the formats
 * and these functions under test: the launcher in `fold2k.cu` adds only
 * the device calls around them.
 *
 *   plan   0   "PTPLAN1\0"
 *          8   u32 degree, u32 base_request, u64 seed
 *          24  u32 n_points, u32 n_orbits, u32 n_reps, u32 canon_bytes
 *          40  u64 canon_tables[canon_bytes * 256]
 *              u64 x[n_points], u64 y[n_points]   (sorted by orbit)
 *              u32 suffix[n_orbits + 1]
 *              u32 rep_orbit[n_reps]
 *              u64 rep_x[n_reps], u64 rep_y[n_reps]
 *
 *   table  0   "PTFOLD1\0"
 *          8   u32 degree, u32 base_request, u64 seed
 *          24  u32 bucket_shift, u32 buckets, u32 words, u32 present_words
 *          40  u64 present_mask
 *          48  u32 canon_bytes, u32 0
 *          56  u64 canon_tables[canon_bytes * 256]
 *              u32 bucket_start[buckets + 1]
 *              u32 words[words]
 *              u64 present[present_words]
 */
#ifndef GPU_ECC2K_FOLD_IO_HPP
#define GPU_ECC2K_FOLD_IO_HPP

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "koblitz.cuh"

/* A point from the low 64 bits of its coordinates, which is all a field
 * of degree at most 62 has. */
inline pt2k pt_from_u64(uint64_t x, uint64_t y) {
    pt2k p;
    memset(&p, 0, sizeof p);
    p.x.v[0] = (uint32_t)x;
    p.x.v[1] = (uint32_t)(x >> 32);
    p.y.v[0] = (uint32_t)y;
    p.y.v[1] = (uint32_t)(y >> 32);
    p.inf = 0;
    return p;
}

struct PtFoldPlan {
    uint32_t degree = 0, base_request = 0;
    uint64_t seed = 0;
    int n_orbits = 0;
    int canon_bytes = 0;
    std::vector<uint64_t> canon_tables;
    std::vector<pt2k> by_orbit, rep_pts;
    std::vector<uint32_t> suffix, rep_orbit;
};

inline bool pt_little_endian() {
    const uint16_t probe = 1;
    return *(const uint8_t *)&probe == 1;
}

/* Reads `count` values of `size` bytes, or fails. */
inline bool pt_read_n(FILE *f, void *out, size_t size, size_t count) {
    return count == 0 || fread(out, size, count, f) == count;
}

/* **Read a plan**, and refuse one a kernel could not index safely: a
 * suffix that is not a non-decreasing run from 0 to `n_points`, or a
 * representative naming an orbit the plan does not have. */
inline bool pt_read_plan(const char *path, PtFoldPlan &plan, std::string &err) {
    if (!pt_little_endian()) {
        err = "big-endian host";
        return false;
    }
    FILE *f = fopen(path, "rb");
    if (!f) {
        err = std::string("cannot open ") + path;
        return false;
    }
    char magic[8];
    uint32_t head[2], counts[4];
    bool ok = pt_read_n(f, magic, 1, 8) && memcmp(magic, "PTPLAN1", 8) == 0 &&
              pt_read_n(f, head, 4, 2) && pt_read_n(f, &plan.seed, 8, 1) &&
              pt_read_n(f, counts, 4, 4);
    if (!ok) {
        fclose(f);
        err = std::string(path) + ": not a plan file";
        return false;
    }
    plan.degree = head[0];
    plan.base_request = head[1];
    const uint32_t n_points = counts[0], n_orbits = counts[1], n_reps = counts[2];
    plan.n_orbits = (int)n_orbits;
    plan.canon_bytes = (int)counts[3];
    plan.canon_tables.resize((size_t)plan.canon_bytes * 256);
    std::vector<uint64_t> x(n_points), y(n_points), rx(n_reps), ry(n_reps);
    plan.suffix.resize((size_t)n_orbits + 1);
    plan.rep_orbit.resize(n_reps);
    ok = pt_read_n(f, plan.canon_tables.data(), 8, plan.canon_tables.size()) &&
         pt_read_n(f, x.data(), 8, n_points) && pt_read_n(f, y.data(), 8, n_points) &&
         pt_read_n(f, plan.suffix.data(), 4, plan.suffix.size()) &&
         pt_read_n(f, plan.rep_orbit.data(), 4, n_reps) && pt_read_n(f, rx.data(), 8, n_reps) &&
         pt_read_n(f, ry.data(), 8, n_reps);
    const bool at_end = ok && fgetc(f) == EOF;
    fclose(f);
    if (!ok || !at_end) {
        err = std::string(path) + (ok ? ": trailing bytes" : ": ends early");
        return false;
    }
    bool sane = plan.suffix[0] == 0 && plan.suffix[n_orbits] == n_points;
    for (uint32_t o = 0; sane && o < n_orbits; o++) sane = plan.suffix[o] <= plan.suffix[o + 1];
    for (uint32_t r = 0; sane && r < n_reps; r++) sane = plan.rep_orbit[r] < n_orbits;
    if (!sane) {
        err = std::string(path) + ": suffix or representatives out of range";
        return false;
    }
    plan.by_orbit.resize(n_points);
    for (uint32_t i = 0; i < n_points; i++) plan.by_orbit[i] = pt_from_u64(x[i], y[i]);
    plan.rep_pts.resize(n_reps);
    for (uint32_t r = 0; r < n_reps; r++) plan.rep_pts[r] = pt_from_u64(rx[r], ry[r]);
    return true;
}

/* A built table, as `PairSumTable::from_folded_parts` takes it. */
struct PtFoldTable {
    int bucket_shift = 0;
    std::vector<uint32_t> bucket_start, words;
    std::vector<uint64_t> present;
    uint64_t present_mask = 0;
};

/* **Write a table**, with the plan's identity and basis beside it. */
inline bool pt_write_table(const char *path, const PtFoldPlan &plan, const PtFoldTable &t,
                           std::string &err) {
    if (!pt_little_endian()) {
        err = "big-endian host";
        return false;
    }
    FILE *f = fopen(path, "wb");
    if (!f) {
        err = std::string("cannot open ") + path + " for writing";
        return false;
    }
    const uint32_t head[2] = {plan.degree, plan.base_request};
    const uint32_t sizes[4] = {(uint32_t)t.bucket_shift, (uint32_t)(t.bucket_start.size() - 1),
                               (uint32_t)t.words.size(), (uint32_t)t.present.size()};
    const uint32_t canon[2] = {(uint32_t)plan.canon_bytes, 0};
    const size_t tables = (size_t)plan.canon_bytes * 256;
    bool ok = fwrite("PTFOLD1", 1, 8, f) == 8 && fwrite(head, 4, 2, f) == 2 &&
              fwrite(&plan.seed, 8, 1, f) == 1 && fwrite(sizes, 4, 4, f) == 4 &&
              fwrite(&t.present_mask, 8, 1, f) == 1 && fwrite(canon, 4, 2, f) == 2 &&
              fwrite(plan.canon_tables.data(), 8, tables, f) == tables &&
              fwrite(t.bucket_start.data(), 4, t.bucket_start.size(), f) == t.bucket_start.size() &&
              fwrite(t.words.data(), 4, t.words.size(), f) == t.words.size() &&
              fwrite(t.present.data(), 8, t.present.size(), f) == t.present.size();
    ok = (fclose(f) == 0) && ok;
    if (!ok) err = std::string("writing ") + path + " failed";
    return ok;
}

#endif /* GPU_ECC2K_FOLD_IO_HPP */
