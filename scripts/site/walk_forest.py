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

Layout: each tree is grown outward from its root with a persistent heading
per branch and a small hash-derived wobble, relaxed with springs and short-
range repulsion, and the trees are then packed onto the canvas largest first.
Pure Python, no dependencies, a few seconds for a couple of thousand nodes.
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
    "edge": "#8a97bb",
    "node_stroke": "#d3dbee",
    "node_fill": "#141a2f",
    "dp": "#4fe8ae",
    "walk_a": "#f8cd78",
    "walk_b": "#6f92ff",
    "shared": "#4fe8ae",
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
# layout
# ---------------------------------------------------------------------------

def wobble(key, salt):
    """A deterministic value in [-1, 1) from an orbit's name."""
    digest = hashlib.blake2b((key + ":" + salt).encode(), digest_size=4).digest()
    return int.from_bytes(digest, "big") / 2 ** 31 - 1.0


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


def subtree_sizes(order, parent):
    size = {node: 1 for node in order}
    for node in reversed(order):
        if parent[node] is not None:
            size[parent[node]] += size[node]
    return size


def grow(forest, root):
    """Initial positions: grow outward from the root with a persistent heading
    per branch that turns slowly, fanning at forks, so long chains meander
    instead of radiating or zigzagging."""
    order, parent = children_of(forest, root)
    size = subtree_sizes(order, parent)
    pos = {root: (0.0, 0.0)}
    heading = {root: None}
    turn = {root: 0.0}
    for node in order:
        kids = forest.pred.get(node, [])
        if not kids:
            continue
        n = len(kids)
        # The turn carries momentum, so the wobble is a slow drift rather
        # than a jitter at every step.  At a fork the branches take their
        # slots in drift order, so a branch that will curl left starts on
        # the left and the two do not cross just past the fork.
        drift = {}
        for kid in kids:
            d = 0.6 * turn[node] + 0.16 * wobble(kid, "turn")
            drift[kid] = max(-0.45, min(0.45, d))
        kids = sorted(kids, key=lambda k: (drift[k], -size[k]))
        for i, kid in enumerate(kids):
            if heading[node] is None:
                angle = wobble(node, "root") * math.pi + 2 * math.pi * i / n
            elif n == 1:
                angle = heading[node] + drift[kid]
            else:
                spread = min(2.2, 0.85 * n)
                angle = heading[node] + spread * ((i + 0.5) / n - 0.5) + drift[kid]
            heading[kid] = angle
            turn[kid] = drift[kid]
            x, y = pos[node]
            pos[kid] = (x + EDGE * math.cos(angle), y + EDGE * math.sin(angle))
    return order, parent, pos


def relax(nodes, parent, pos, iterations, repel_radius=2.0 * EDGE, step=0.35):
    """Springs along edges, short-range repulsion between nodes that are not
    joined, and a little stiffness that keeps chains smooth, on a grid so
    it stays linear in the node count."""
    if len(nodes) < 2:
        return pos
    cell = repel_radius
    pos = dict(pos)
    kids = {node: [] for node in nodes}
    for node in nodes:
        if parent[node] is not None:
            kids[parent[node]].append(node)
    for it in range(iterations):
        temp = step * (1.0 - it / iterations) + 0.05
        grid = {}
        for node in nodes:
            x, y = pos[node]
            grid.setdefault((int(x // cell), int(y // cell)), []).append(node)
        force = {node: [0.0, 0.0] for node in nodes}
        for node in nodes:
            x, y = pos[node]
            cx, cy = int(x // cell), int(y // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for other in grid.get((cx + dx, cy + dy), ()):
                        if other is node or other < node:
                            continue
                        if parent[other] is node or parent[node] is other:
                            continue
                        ox, oy = pos[other]
                        vx, vy = x - ox, y - oy
                        d2 = vx * vx + vy * vy
                        if d2 >= repel_radius * repel_radius:
                            continue
                        d = math.sqrt(d2) or 1e-6
                        push = (repel_radius - d) / repel_radius
                        push = push * push * EDGE * 0.9
                        fx, fy = vx / d * push, vy / d * push
                        force[node][0] += fx
                        force[node][1] += fy
                        force[other][0] -= fx
                        force[other][1] -= fy
        for node in nodes:
            p = parent[node]
            if p is None:
                continue
            x, y = pos[node]
            px, py = pos[p]
            vx, vy = px - x, py - y
            d = math.hypot(vx, vy) or 1e-6
            pull = (d - EDGE) * 0.5
            fx, fy = vx / d * pull, vy / d * pull
            force[node][0] += fx
            force[node][1] += fy
            force[p][0] -= fx
            force[p][1] -= fy
            # Stiffness: a node inside a chain is drawn toward the midpoint
            # of its two neighbours.
            if len(kids[node]) == 1:
                kx, ky = pos[kids[node][0]]
                force[node][0] += ((px + kx) / 2 - x) * 0.3
                force[node][1] += ((py + ky) / 2 - y) * 0.3
        for node in nodes:
            fx, fy = force[node]
            mag = math.hypot(fx, fy)
            if mag > EDGE:
                fx, fy = fx / mag * EDGE, fy / mag * EDGE
            x, y = pos[node]
            pos[node] = (x + fx * temp, y + fy * temp)
    return pos


def bbox(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def pack(trees, cell=1.7 * EDGE):
    """Place trees largest first, each at the first spot on a spiral from the
    centre where none of its nodes lands in a cell another tree occupies.
    Cell occupancy rather than bounding boxes lets a small tree settle into
    the bends of a large one, which is what keeps the drawing compact."""
    occupied = set()
    offsets = {}

    def cells(points, tx, ty):
        out = set()
        for x, y in points:
            cx, cy = int((x + tx) // cell), int((y + ty) // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    out.add((cx + dx, cy + dy))
        return out

    for root, points in trees:
        x0, y0, x1, y1 = bbox(points)
        cx, cy = -(x0 + x1) / 2, -(y0 + y1) / 2   # translation that centres the tree
        angle = wobble(root, "spiral") * math.pi
        k = 0
        while True:
            r = 1.3 * EDGE * math.sqrt(k)
            a = angle + k * 2.399963
            tx, ty = cx + r * math.cos(a), cy + r * math.sin(a)
            mine = cells(points, tx, ty)
            if not (mine & occupied):
                offsets[root] = (tx, ty)
                # A tree claims only its own cells; the dilation above is what
                # keeps the next one a cell away.
                occupied |= {(int((x + tx) // cell), int((y + ty) // cell)) for x, y in points}
                break
            k += 1
    return offsets


def layout(forest, iterations=160):
    positions = {}
    trees = []
    for root in forest.roots:
        order, parent, pos = grow(forest, root)
        pos = relax(order, parent, pos, iterations)
        trees.append((root, order, parent, pos))
    trees.sort(key=lambda t: -len(t[1]))
    offsets = pack([(root, list(pos.values())) for root, order, parent, pos in trees])
    for root, order, parent, pos in trees:
        ox, oy = offsets[root]
        for node in order:
            x, y = pos[node]
            positions[node] = (x + ox, y + oy)
    return positions


# ---------------------------------------------------------------------------
# svg
# ---------------------------------------------------------------------------

def fmt(v):
    text = "%.1f" % v
    return text[:-2] if text.endswith(".0") else text


def render(forest, positions, meetings, width, height, checked):
    pad = 2 * EDGE
    x0, y0, x1, y1 = bbox(list(positions.values()))
    scale = min((width - 2 * pad) / max(x1 - x0, 1.0), (height - 2 * pad) / max(y1 - y0, 1.0))
    scale = min(scale, 1.0)   # never blow a small forest up

    def sx(x):
        return (x - x0) * scale + (width - (x1 - x0) * scale) / 2

    def sy(y):
        return (y - y0) * scale + (height - (y1 - y0) * scale) / 2

    highlighted = set()
    lines = []
    p = forest.params
    lines.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
                 'role="img" aria-labelledby="wf-title wf-desc">' % (width, height, width, height))
    lines.append('<title id="wf-title">The forest %s walks on GF(2^%s) build on their way to distinguished points</title>'
                 % (len(forest.walks), p.get("curve", "?")))
    lines.append('<desc id="wf-desc">%d orbits drawn, one every %d iterations, over %d iterations in all; '
                 '%d distinguished points, %d meetings between walks from different starts. Generated by '
                 'scripts/site/walk_forest.py from trails written by ecc2k130/src/trailforest.cpp. %s</desc>' % (
                     len(positions), forest.every, total_steps(forest), len(forest.roots),
                     count_meetings(forest), forest.header))
    lines.append("<!-- %s; corpus records checked: %d -->" % (forest.header, checked))
    # An <img> evaluates the SVG's own media queries against the box it is
    # drawn in, not the page. At phone width the 1000-unit viewBox is ~360px
    # and a 0.8-unit stroke is a third of a pixel; thicken so the trails
    # stay visible without forcing the page to scroll sideways.
    lines.append("<style>"
                 ".e{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round}"
                 ".n{fill:%s;stroke:%s;stroke-width:%s}"
                 ".dp{fill:%s;stroke:%s;stroke-width:%s}"
                 ".a{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".b{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".s{stroke:%s;stroke-width:%s;fill:none;stroke-linecap:round;stroke-linejoin:round}"
                 ".na{fill:%s;stroke:%s;stroke-width:%s}"
                 ".nb{fill:%s;stroke:%s;stroke-width:%s}"
                 "@media(max-width:720px){"
                 ".e{stroke-width:2.8}"
                 ".n,.dp,.na,.nb{stroke-width:2.2}"
                 ".a,.b,.s{stroke-width:4.2}"
                 "circle{r:3.4px}"
                 "circle.dp{r:5px}"
                 "}"
                 "</style>" % (
                     COLORS["edge"], fmt(STROKE * scale),
                     COLORS["node_fill"], COLORS["node_stroke"], fmt(STROKE * scale),
                     COLORS["dp"], COLORS["dp"], fmt(STROKE * scale),
                     COLORS["walk_a"], fmt(HIGHLIGHT_STROKE * scale),
                     COLORS["walk_b"], fmt(HIGHLIGHT_STROKE * scale),
                     COLORS["shared"], fmt(HIGHLIGHT_STROKE * scale),
                     COLORS["walk_a"], COLORS["walk_a"], fmt(STROKE * scale),
                     COLORS["walk_b"], COLORS["walk_b"], fmt(STROKE * scale)))

    # Plain edges, one path per tree.
    lines.append('<g class="e">')
    for root in forest.roots:
        order, parent = children_of(forest, root)
        d = []
        for node in order:
            par = parent[node]
            if par is None:
                continue
            ax, ay = positions[par]
            bx, by = positions[node]
            d.append("M%s %sL%s %s" % (fmt(sx(ax)), fmt(sy(ay)), fmt(sx(bx)), fmt(sy(by))))
        if d:
            lines.append('<path d="%s"/>' % "".join(d))
    lines.append("</g>")

    # Highlighted trails on top of the plain edges.
    node_class = {}
    for m in meetings:
        for cls, trail in (("a", m["a"]), ("b", m["b"]), ("s", m["shared"])):
            pts = ["%s %s" % (fmt(sx(positions[o][0])), fmt(sy(positions[o][1]))) for o in trail]
            lines.append('<path class="%s" d="M%s"/>' % (cls, "L".join(pts)))
            if cls in ("a", "b"):
                for o in trail[:-1]:
                    node_class.setdefault(o, "n" + cls)
        highlighted.add(m["orbit"])

    lines.append("<g>")
    roots = set(forest.roots)
    for node in forest.order:
        x, y = positions[node]
        if node in roots:
            lines.append('<circle class="dp" cx="%s" cy="%s" r="%s"/>' % (fmt(sx(x)), fmt(sy(y)), fmt(DP_R * scale)))
        else:
            cls = node_class.get(node, "n")
            lines.append('<circle class="%s" cx="%s" cy="%s" r="%s"/>' % (cls, fmt(sx(x)), fmt(sy(y)), fmt(NODE_R * scale)))
    lines.append("</g>")
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def count_meetings(forest):
    return sum(1 for preds in forest.pred.values() if len(preds) > 1)


def total_steps(forest):
    return sum(forest.steps)


# ---------------------------------------------------------------------------

def build(trails_path, corpus_path=None, width=1000, height=900, iterations=160):
    forest = read_trails(trails_path)
    checked = bind_to_corpus(forest, read_corpus(corpus_path)) if corpus_path else 0
    positions = layout(forest, iterations)
    meetings = find_meetings(forest, HIGHLIGHTS)
    return render(forest, positions, meetings, width, height, checked), forest


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
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--iterations", type=int, default=160, help="relaxation steps per tree")
    parser.add_argument("--check", action="store_true", help="fail if --out would change instead of writing it")
    parser.add_argument("--hash-trails", type=int, default=0, metavar="EVERY",
                        help="instead of drawing, write --trails re-sampled every EVERY steps with orbits hashed to --out")
    args = parser.parse_args(argv)

    if args.hash_trails:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(hash_trails(read_trails(args.trails), args.hash_trails))
        return 0

    svg, forest = build(args.trails, args.corpus, args.width, args.height, args.iterations)
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
