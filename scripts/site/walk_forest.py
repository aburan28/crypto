#!/usr/bin/env python3
"""Draw the forest that the campaign's walks build, from their trails.

Input is the output of ecc2k130/src/trailforest.cpp: one line per walk,
listing every orbit the walk passed through on its way to its distinguished
point.  Walks are deterministic in the point they stand on, so two walks that
reach the same orbit share every step after it; the trails therefore form a
forest whose roots are the distinguished points, and this script draws that
forest as an SVG for the campaign status page.

    python3 scripts/site/walk_forest.py \
        --trails docs/ecc2k130-status/walk-forest/trails.txt \
        --corpus docs/ecc2k130-status/walk-forest/forest.bin \
        --out docs/ecc2k130-status/walk-forest.svg

What is drawn is exactly what is in the trails file: one hollow node per
orbit visited, one edge per iteration, a filled node at each distinguished
point.  The first few places where two walks from different starts meet are
picked out in colour, one colour per arriving walk, with the tail they then
share drawn in the distinguished-point colour.

`--corpus` binds the figure to the client's record format: every trail must
end on the orbit the corresponding 32-byte record names, and the script
refuses to draw otherwise.  Nothing here is random; the layout is a
deterministic function of the trails, so the SVG regenerates identically.

Layout: a circle.  Every distinguished point sits on the rim and every walk
comes in to it from inside, as far from the rim as it had iterations to go,
so a walk's seed sits as deep as the walk was long; each walk has its own
slice of the circle and its own colour.  Pure Python, no dependencies.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import struct
import sys

RECORD = struct.Struct("<4Q")  # seed, then the canonical orbit as three limbs

# Canvas units.  Edge length sets the scale of everything else.
EDGE = 10.0
NODE_R = 2.6
DP_R = 3.6
STROKE = 1.1
HIGHLIGHT_STROKE = 2.4
HIGHLIGHTS = 4

# The dashboard palette, baked in: an <img> cannot read the page's CSS
# variables, and the page is dark only.
COLORS = {
    "node_fill": "#141a2f",
    "dp": "#4fe8ae",
    "walk_a": "#f8cd78",
    "walk_b": "#6f92ff",
    "shared": "#4fe8ae",
    "ring": "#2a3354",
    "tick": "#7d89aa",
    "bg": "#141a2f",
}


# ---------------------------------------------------------------------------
# trails
# ---------------------------------------------------------------------------

class Forest:
    def __init__(self):
        self.header = ""
        self.params = {}
        self.walks = []       # [(seed_hex, [orbit, ...])]
        self.succ = {}        # orbit -> next orbit
        self.pred = {}        # orbit -> [orbit, ...] in walk order
        self.roots = []       # distinguished orbits, in walk order
        self.order = []       # every orbit, first-seen order
        self.every = 1        # iterations between consecutive drawn orbits
        self.steps = []       # iterations per walk, in walk order


def expected_nodes(steps, every):
    """A trail sampled every `every` steps from step 0, plus its endpoint
    when that did not fall on a stride."""
    return steps // every + 1 + (1 if steps % every else 0)


def parse_header(line):
    # "# curve 23 instance 0 dp-weight 7 mode generate run-id 1 walks 96 checked 6"
    words = line[1:].split()
    return {words[i]: words[i + 1] for i in range(0, len(words) - 1, 2)}


def read_trails(path):
    forest = Forest()
    seen = set()
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if not forest.header:
                    forest.header = line
                    forest.params = parse_header(line)
                    forest.every = int(forest.params.get("every", 1))
                continue
            words = line.split()
            if words[0] != "walk" or len(words) < 4:
                raise SystemExit("unrecognised trail line: %r" % line[:80])
            seed, steps, orbits = words[1], int(words[2]), words[3:]
            if len(orbits) != expected_nodes(steps, forest.every):
                raise SystemExit("walk %s: %d steps but %d orbits at a stride of %d"
                                 % (seed, steps, len(orbits), forest.every))
            for a, b in zip(orbits, orbits[1:]):
                if a in forest.succ and forest.succ[a] != b:
                    raise SystemExit("orbit %s steps to two different orbits; the walk is not deterministic" % a)
                if a not in forest.succ:
                    forest.succ[a] = b
                    forest.pred.setdefault(b, []).append(a)
            for orbit in orbits:
                if orbit not in seen:
                    seen.add(orbit)
                    forest.order.append(orbit)
            root = orbits[-1]
            if root in forest.succ:
                raise SystemExit("seed %s ends on %s, which another walk stepped through" % (seed, root))
            if root not in forest.roots:
                forest.roots.append(root)
            forest.walks.append((seed, orbits))
            forest.steps.append(steps)
    if not forest.walks:
        raise SystemExit("no walks in %s" % path)
    return forest


def read_corpus(path):
    """The endpoints trails must land on, keyed the way the trails name
    walks: a client corpus (32-byte records, keyed by seed) or, for a
    sampled forest whose orbits are hashed, trailforest's --hashes-out file
    of `<walk> <name>` lines."""
    if path.endswith(".hashes"):
        out = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                words = line.split()
                if len(words) == 2:
                    out[words[0]] = words[1]
        return out
    with open(path, "rb") as fh:
        data = fh.read()
    if len(data) % RECORD.size:
        raise SystemExit("%s is not a whole number of %d-byte records" % (path, RECORD.size))
    out = {}
    for i in range(0, len(data), RECORD.size):
        seed, k0, k1, k2 = RECORD.unpack_from(data, i)
        out["%016x" % seed] = "%x%016x%016x" % (k2, k1, k0)
    return out


def orbit_name(orbit_hex):
    """trailforest's public name for an orbit: the first 16 hex digits of
    SHA-256 over the canonical representative's 24 little-endian bytes."""
    return hashlib.sha256(int(orbit_hex, 16).to_bytes(24, "little")).hexdigest()[:16]


def canon_hex(text):
    # trailforest prints the 192-bit value with leading zeros dropped past the
    # first digit; the corpus reader above prints it in full.
    return text.lstrip("0") or "0"


def bind_to_corpus(forest, corpus):
    """Every trail must end on the orbit its record names."""
    hashed = "hash" in forest.params
    checked = 0
    for seed, orbits in forest.walks:
        if seed not in corpus:
            raise SystemExit("walk %s has a trail but no corpus record" % seed)
        recorded, reached = corpus[seed], orbits[-1]
        if not hashed:
            recorded, reached = canon_hex(recorded), canon_hex(reached)
        if recorded != reached:
            raise SystemExit("walk %s: the trail ends on %s but the record names %s"
                             % (seed, orbits[-1], corpus[seed]))
        checked += 1
    return checked


# ---------------------------------------------------------------------------
# highlights: where walks from different starts meet
# ---------------------------------------------------------------------------

def find_meetings(forest, limit):
    """The first `limit` orbits, in walk order, that two walks arrived at from
    different predecessors, each with the two arriving trails and the shared
    tail to the distinguished point."""
    arrived = {}   # orbit -> (walk index, position) of the first walk through it
    meetings = []
    for w, (seed, orbits) in enumerate(forest.walks):
        for pos, orbit in enumerate(orbits):
            if orbit not in arrived:
                arrived[orbit] = (w, pos)
                continue
            if pos == 0:
                break   # a start orbit another walk already visited: no arrival to show
            first_w, first_pos = arrived[orbit]
            if first_pos == 0 or orbits[pos - 1] == forest.walks[first_w][1][first_pos - 1]:
                break   # same predecessor: the walks had already merged, or the first started here
            meetings.append({
                "orbit": orbit,
                "a": forest.walks[first_w][1][:first_pos + 1],
                "b": orbits[:pos + 1],
                "shared": orbits[pos:],
            })
            break   # from here on this walk is the earlier one's
        if len(meetings) >= limit:
            break
    return meetings


# ---------------------------------------------------------------------------
# layout: a circle
# ---------------------------------------------------------------------------
#
# Every distinguished point sits on the rim and every walk comes in to it
# from inside: a node's distance from the rim is the number of iterations the
# walk still had to take from there, so a walk's seed sits as far in as the
# walk was long and every walker heads outward at the same rate.  Each walk
# owns an equal slice of the circle, in the order of its start in a
# depth-first pass over its tree, so a tree's walks sit side by side and
# their branches join without crossing; a node shared by several walks sits
# at the mean of their slices.  A small swirl, the same for every trail,
# turns the spokes into arcs that read as motion toward the rim.

RIM = 450.0          # radius of the rim, where the distinguished points sit
HUB = 36.0           # radius of `depth` iterations to go
SWIRL = 0.9          # radians of turn across the full depth
RINGS = 4            # guide rings at every quarter of the depth


def children_of(forest, root):
    """Depth-first order of a tree with the parent of each node."""
    order = [root]
    parent = {root: None}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in forest.pred.get(node, ()):
            parent[child] = node
            order.append(child)
            stack.append(child)
    return order, parent


def bbox(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def walk_hue(slot, walks):
    """A colour per walk, from its place on the circle: gold at the top
    through rose and violet to blue at the bottom and back, so neighbours
    differ and there is no seam.  Greens are left to the distinguished
    points."""
    t = (slot + 0.5) / walks
    tri = 1.0 - abs(2.0 * t - 1.0)
    return round((45.0 - 180.0 * tri) % 360.0, 1)


def layout(forest):
    """Positions centred on the origin, and the circle's geometry: `depth`
    is the cap the walks were sampled under, or the longest walk when there
    was none, and maps to the hub."""
    depth = int(forest.params.get("cap", 0)) or max(forest.steps)
    remaining = {}
    for (seed, orbits), steps in zip(forest.walks, forest.steps):
        last = len(orbits) - 1
        for i, orbit in enumerate(orbits):
            remaining.setdefault(orbit, steps - i * forest.every if i < last else 0)

    preorder = {}
    for r, root in enumerate(forest.roots):
        order, parent = children_of(forest, root)
        for i, node in enumerate(order):
            preorder.setdefault(node, (r, i))
    slots = sorted(range(len(forest.walks)), key=lambda w: (preorder[forest.walks[w][1][-1]][0],
                                                            preorder[forest.walks[w][1][0]][1], w))
    slot_of = {w: k for k, w in enumerate(slots)}
    n = len(forest.walks)
    through = {}
    for w, (seed, orbits) in enumerate(forest.walks):
        angle = -math.pi / 2 + 2 * math.pi * (slot_of[w] + 0.5) / n
        for orbit in orbits:
            through.setdefault(orbit, []).append(angle)

    positions = {}
    for node in forest.order:
        base = sum(through[node]) / len(through[node])
        frac = min(1.0, remaining[node] / depth)
        radius = RIM - (RIM - HUB) * frac
        angle = base + SWIRL * frac
        positions[node] = (radius * math.cos(angle), radius * math.sin(angle))
    geometry = {
        "rim": RIM,
        "hub": HUB,
        "depth": depth,
        "swirl": SWIRL,
        "hues": [walk_hue(slot_of[w], n) for w in range(n)],
    }
    return positions, geometry


def ring_radius(k):
    """Radius of the k-th guide ring, k quarters of the depth to go."""
    return RIM - (RIM - HUB) * k / RINGS


# ---------------------------------------------------------------------------
# svg
# ---------------------------------------------------------------------------

def fmt(v):
    text = "%.1f" % v
    return text[:-2] if text.endswith(".0") else text


def render(forest, positions, geometry, meetings, width, height, checked):
    pad = 2 * EDGE
    scale = min(1.0, (min(width, height) - 2 * pad) / (2 * RIM))
    cx, cy = width / 2, height / 2

    def sx(x):
        return cx + x * scale

    def sy(y):
        return cy + y * scale

    lines = []
    p = forest.params
    lines.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
                 'role="img" aria-labelledby="wf-title wf-desc">' % (width, height, width, height))
    lines.append('<title id="wf-title">The forest %s walks on GF(2^%s) build on their way to distinguished points</title>'
                 % (len(forest.walks), p.get("curve", "?")))
    lines.append('<desc id="wf-desc">%d orbits drawn, one every %d iterations, over %d iterations in all; '
                 '%d distinguished points, %d meetings between walks from different starts. Drawn as a circle: '
                 'each distinguished point on the rim, each walk coming in to it from as far inside as it had '
                 'iterations to go, up to %d at the hub. Generated by scripts/site/walk_forest.py from trails '
                 'written by ecc2k130/src/trailforest.cpp. %s</desc>' % (
                     len(positions), forest.every, total_steps(forest), len(forest.roots),
                     count_meetings(forest), geometry["depth"], forest.header))
    lines.append("<!-- %s; corpus records checked: %d -->" % (forest.header, checked))
    # An <img> evaluates the SVG's own media queries against the box it is
    # drawn in, not the page. At phone width the 1000-unit viewBox is ~360px
    # and a 0.8-unit stroke is a third of a pixel; thicken so the trails
    # stay visible without forcing the page to scroll sideways.
    lines.append("<style>"
                 ".e path{stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".ring{stroke:%s;stroke-width:%s;fill:none}"
                 ".rim{stroke:%s;stroke-width:%s;fill:none;opacity:.45}"
                 ".tick{fill:%s;font:11px ui-monospace,Menlo,monospace;text-anchor:middle;paint-order:stroke;stroke:%s;stroke-width:4px}"
                 ".n{fill:%s;stroke-width:%s}"
                 ".dp{fill:%s;stroke:%s;stroke-width:%s}"
                 ".a{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".b{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".s{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 "@media(max-width:720px){"
                 ".e path{stroke-width:2.8}"
                 ".n,.dp{stroke-width:2.2}"
                 ".a,.b,.s{stroke-width:4.2}"
                 ".tick{display:none}"
                 "circle{r:3.4px}"
                 "circle.dp{r:5px}"
                 "}"
                 "</style>" % (
                     fmt(STROKE * scale),
                     COLORS["ring"], fmt(STROKE * scale),
                     COLORS["dp"], fmt(STROKE * scale),
                     COLORS["tick"], COLORS["bg"],
                     COLORS["node_fill"], fmt(STROKE * scale),
                     COLORS["dp"], COLORS["dp"], fmt(STROKE * scale),
                     COLORS["walk_a"], fmt(HIGHLIGHT_STROKE * scale),
                     COLORS["walk_b"], fmt(HIGHLIGHT_STROKE * scale),
                     COLORS["shared"], fmt(HIGHLIGHT_STROKE * scale)))

    # Guide rings, as paths so every <circle> in the drawing is an orbit:
    # the rim where the walks end, and a ring at each quarter of the depth
    # labelled, on the right where the short walks leave room, with the
    # iterations still to go from it.
    def ring(radius):
        r = radius * scale
        return "M%s %sa%s %s 0 1 0 %s 0a%s %s 0 1 0 %s 0" % (
            fmt(cx - r), fmt(cy), fmt(r), fmt(r), fmt(2 * r), fmt(r), fmt(r), fmt(-2 * r))
    lines.append('<path class="rim" d="%s"/>' % ring(RIM))
    lines.append('<path class="ring" d="%s"/>' % "".join(ring(ring_radius(k)) for k in range(1, RINGS)))
    for k in range(1, RINGS):
        lines.append('<text class="tick" x="%s" y="%s">%s to go</text>' % (
            fmt(sx(ring_radius(k))), fmt(cy + 4), "{:,}".format(geometry["depth"] * k // RINGS)))

    # Edges, one path per walk in the walk's own colour; a stride two walks
    # share is drawn once, in the colour of the first to take it.
    by_walk = {}
    for w, (seed, orbits) in enumerate(forest.walks):
        for a, b in zip(orbits, orbits[1:]):
            if a not in by_walk:
                by_walk[a] = w
    lines.append('<g class="e">')
    for w in range(len(forest.walks)):
        d = []
        for a in forest.walks[w][1][:-1]:
            if by_walk.get(a) != w:
                continue
            b = forest.succ[a]
            ax, ay = positions[a]
            bx, by = positions[b]
            d.append("M%s %sL%s %s" % (fmt(sx(ax)), fmt(sy(ay)), fmt(sx(bx)), fmt(sy(by))))
        if d:
            lines.append('<path stroke="%s" d="%s"/>' % (hsl(geometry["hues"][w], 72, 66), "".join(d)))
    lines.append("</g>")

    # Where walks from different starts meet, on top of the plain edges.
    for m in meetings:
        for cls, trail in (("a", m["a"]), ("b", m["b"]), ("s", m["shared"])):
            pts = ["%s %s" % (fmt(sx(positions[o][0])), fmt(sy(positions[o][1]))) for o in trail]
            lines.append('<path class="%s" d="M%s"/>' % (cls, "L".join(pts)))

    # Nodes: hollow in the walk's colour, filled green at the rim.
    first = {}
    for w, (seed, orbits) in enumerate(forest.walks):
        for orbit in orbits:
            first.setdefault(orbit, w)
    lines.append("<g>")
    roots = set(forest.roots)
    for node in forest.order:
        x, y = positions[node]
        if node in roots:
            lines.append('<circle class="dp" cx="%s" cy="%s" r="%s"/>' % (fmt(sx(x)), fmt(sy(y)), fmt(DP_R * scale)))
        else:
            lines.append('<circle class="n" stroke="%s" cx="%s" cy="%s" r="%s"/>' % (
                hsl(geometry["hues"][first[node]], 60, 80), fmt(sx(x)), fmt(sy(y)), fmt(NODE_R * scale)))
    lines.append("</g>")
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def hsl(hue, sat, light):
    return "hsl(%s %d%% %d%%)" % (fmt(hue), sat, light)


def count_meetings(forest):
    return sum(1 for preds in forest.pred.values() if len(preds) > 1)


def total_steps(forest):
    return sum(forest.steps)


# ---------------------------------------------------------------------------

def build(trails_path, corpus_path=None, width=1000, height=1000):
    forest = read_trails(trails_path)
    checked = bind_to_corpus(forest, read_corpus(corpus_path)) if corpus_path else 0
    positions, geometry = layout(forest)
    meetings = find_meetings(forest, HIGHLIGHTS)
    return render(forest, positions, geometry, meetings, width, height, checked), forest


def export_graph(forest, positions, geometry, meetings, checked, title):
    """The same forest as a graph the page can explore: every drawn orbit a
    node with its layout position, every stride an edge, every walk the
    ordered list of its nodes and its colour.  Node names are whatever the
    trails carry (hash prefixes on the challenge curve), so nothing leaves
    the trails file that was not already in it.  Positions are the static
    figure's, measured from the circle's top-left corner so the circle's
    centre is (rim, rim), and the page needs no layout of its own."""
    import json
    index = {node: i for i, node in enumerate(forest.order)}
    roots = set(forest.roots)
    nodes = []
    for node in forest.order:
        x, y = positions[node]
        nodes.append([node, round(x + RIM, 1), round(y + RIM, 1), 1 if node in roots else 0])
    edges = [[index[a], index[b]] for a, b in forest.succ.items()]
    walks = []
    for w, ((seed, orbits), steps) in enumerate(zip(forest.walks, forest.steps)):
        walks.append({"id": seed, "steps": steps, "hue": geometry["hues"][w], "nodes": [index[o] for o in orbits]})
    highlights = [{"orbit": index[m["orbit"]], "a": [index[o] for o in m["a"]],
                   "b": [index[o] for o in m["b"]], "shared": [index[o] for o in m["shared"]]}
                  for m in meetings]
    graph = {
        "title": title,
        "header": forest.header,
        "params": forest.params,
        "every": forest.every,
        "checked": checked,
        "counts": {
            "walks": len(forest.walks),
            "nodes": len(forest.order),
            "edges": len(edges),
            "distinguished": len(forest.roots),
            "meetings": count_meetings(forest),
            "iterations": total_steps(forest),
            "longest": max(forest.steps),
        },
        "circle": {"rim": RIM, "hub": HUB, "depth": geometry["depth"], "rings": RINGS},
        "width": 2 * RIM,
        "height": 2 * RIM,
        "nodes": nodes,
        "edges": edges,
        "walks": walks,
        "highlights": highlights,
    }
    return json.dumps(graph, separators=(",", ":"), sort_keys=True) + "\n"


def build_graph(trails_path, corpus_path=None, title=""):
    forest = read_trails(trails_path)
    checked = bind_to_corpus(forest, read_corpus(corpus_path)) if corpus_path else 0
    positions, geometry = layout(forest)
    meetings = find_meetings(forest, HIGHLIGHTS)
    return export_graph(forest, positions, geometry, meetings, checked, title), forest


def hash_trails(forest, every):
    """Rewrite a raw (every step, orbits in the clear) trails file the way
    trailforest --sample writes one: sampled every `every` steps plus the
    endpoint, every orbit replaced by its public name.  This is how a scalar
    reference replay is compared against a sampled kernel replay without
    the orbits themselves leaving the machine that holds them."""
    if forest.every != 1:
        raise SystemExit("--hash-trails wants a raw trails file, not one already sampled")
    lines = ["# %s every %d hash sha256-16" % (forest.header[1:].strip(), every)]
    for index, (seed, orbits) in enumerate(forest.walks):
        steps = len(orbits) - 1
        picked = orbits[::every]
        if steps % every:
            picked.append(orbits[-1])
        lines.append("walk %d %d %s" % (index, steps, " ".join(orbit_name(o) for o in picked)))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description="Draw the walk forest for the status page.")
    parser.add_argument("--trails", required=True, help="trailforest output")
    parser.add_argument("--corpus", default=None,
                        help="what the trails must end on: 32-byte corpus records, or a .hashes file for sampled trails")
    parser.add_argument("--out", required=True, help="SVG to write")
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument("--check", action="store_true", help="fail if --out would change instead of writing it")
    parser.add_argument("--hash-trails", type=int, default=0, metavar="EVERY",
                        help="instead of drawing, write --trails re-sampled every EVERY steps with orbits hashed to --out")
    parser.add_argument("--json", action="store_true",
                        help="write the forest as a graph (nodes, edges, walks, layout) for the page's explorer instead of an SVG")
    parser.add_argument("--title", default="", help="the graph's title, with --json")
    args = parser.parse_args(argv)

    if args.hash_trails:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(hash_trails(read_trails(args.trails), args.hash_trails))
        return 0

    if args.json:
        text, forest = build_graph(args.trails, args.corpus, args.title)
        if args.check:
            with open(args.out, encoding="utf-8") as fh:
                if fh.read() != text:
                    raise SystemExit("%s is not what %s exports to; regenerate it" % (args.out, args.trails))
            print("%s is up to date" % args.out)
            return 0
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print("%s: %d walks, %d nodes -> %s" % (args.trails, len(forest.walks), len(forest.order), args.out), file=sys.stderr)
        return 0

    svg, forest = build(args.trails, args.corpus, args.width, args.height)
    if args.check:
        if not os.path.exists(args.out):
            raise SystemExit("%s does not exist" % args.out)
        with open(args.out, encoding="utf-8") as fh:
            if fh.read() != svg:
                raise SystemExit("%s is not what %s renders to; regenerate it" % (args.out, args.trails))
        print("%s is up to date" % args.out)
        return 0
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(svg)
    print("%s: %d walks, %d orbits, %d distinguished points, %d meetings -> %s" % (
        args.trails, len(forest.walks), len(forest.order), len(forest.roots), count_meetings(forest), args.out),
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
