// The walk forest as a graph you can explore.
//
// The static figure (walk-forest.svg) is drawn by scripts/site/walk_forest.py
// from the committed trails; this loads the same forest exported as a graph
// (walk-forest.json, and walk-forest-gf2-23.json for the test curve where
// trails merge) with the same layout, and puts it on a canvas. The layout is
// a circle: every distinguished point on the rim, every walk coming in to it
// from as far inside as it had iterations to go, each in its own colour.
// Drag to pan, wheel or pinch to zoom, hover a node to see which walk it is
// on and how far along, click a node to light up every path through it to
// its distinguished point, and press play to run every walk at once on one
// clock of iterations, each lighting its trail as it goes and arriving at
// the rim when its walk did.
//
// Progressive enhancement: the <img> stays until the graph has loaded and
// drawn, and stays if it cannot. Nothing here fetches anything but the two
// JSON files beside the page, and the node names in them are hash prefixes
// on the challenge curve, so the explorer shows exactly what the figure shows.
(function () {
  "use strict";

  var host = document.getElementById("forest");
  if (!host || !window.fetch || !document.createElement("canvas").getContext) return;

  var COLORS = {
    edgeDim: "#232b47",
    node: "#141a2f",
    ringDim: "#343d5e",
    guide: "#2a3354",
    tick: "#7d89aa",
    dp: "#4fe8ae",
    dpDim: "#2b6d55",
    walkA: "#f8cd78",
    walkB: "#6f92ff",
    walker: "#ffffff",
    label: "#e8edf7",
    bg: "#0b1020"
  };

  // One pass of every walk, rim to hub, takes this long however deep the
  // circle is; a walk's speed in iterations is the same for all of them.
  var PLAY_MS = 9000;
  var ARRIVE_MS = 700;

  var DATASETS = [
    { key: "real", file: "./walk-forest.json", label: "ECC2K-130, real walks" },
    { key: "test", file: "./walk-forest-gf2-23.json", label: "GF(2^23) test curve, where trails merge" }
  ];

  var canvas = document.createElement("canvas");
  canvas.className = "forest-canvas";
  canvas.setAttribute("aria-label", "The walk forest as an explorable graph");
  var ctx = canvas.getContext("2d");
  var controls = host.querySelector(".forest-controls");
  var select = host.querySelector("#forest-dataset");
  var playButton = host.querySelector("#forest-play");
  var resetButton = host.querySelector("#forest-reset");
  var readout = host.querySelector("#forest-readout");
  var image = host.querySelector("img");
  var tip = document.createElement("div");
  tip.className = "forest-tip";
  tip.hidden = true;

  var graph = null;          // the loaded dataset
  var view = { x: 0, y: 0, k: 1 };
  var dpr = Math.max(1, Math.min(2, window.devicePixelRatio || 1));
  var hover = -1;
  var selected = -1;
  var lit = null;            // Set of node indices on the selected paths
  var litEdges = null;       // Set of "a-b" keys
  var clock = null;          // iterations elapsed while playing, else null
  var started = 0;           // performance.now() at iteration 0
  var arrived = null;        // per walk: time it reached the rim, or 0
  var playing = false;
  var needsDraw = true;
  var succ = null, preds = null, walkOf = null, posOf = null;
  var owner = null;          // per node: the first walk through it, whose colour it takes
  var hue = null, hueDim = null, hueNode = null;

  function size() {
    var rect = host.getBoundingClientRect();
    var w = Math.max(300, Math.floor(rect.width));
    var h = Math.max(320, Math.min(860, Math.round(w * 0.92)));
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    canvas.style.width = w + "px";
    canvas.style.height = h + "px";
    needsDraw = true;
  }

  function fit() {
    if (!graph) return;
    var w = canvas.width / dpr, h = canvas.height / dpr;
    var pad = 16;
    view.k = Math.min((w - 2 * pad) / Math.max(graph.width, 1), (h - 2 * pad) / Math.max(graph.height, 1));
    view.x = (w - graph.width * view.k) / 2;
    view.y = (h - graph.height * view.k) / 2;
    needsDraw = true;
  }

  function toScreen(i) {
    var n = graph.nodes[i];
    return [n[1] * view.k + view.x, n[2] * view.k + view.y];
  }

  function index(g) {
    succ = new Array(g.nodes.length);
    preds = new Array(g.nodes.length);
    walkOf = new Array(g.nodes.length);
    posOf = new Array(g.nodes.length);
    owner = new Array(g.nodes.length);
    var i;
    for (i = 0; i < g.nodes.length; i++) { succ[i] = -1; preds[i] = []; walkOf[i] = []; posOf[i] = []; owner[i] = -1; }
    for (i = 0; i < g.edges.length; i++) {
      succ[g.edges[i][0]] = g.edges[i][1];
      preds[g.edges[i][1]].push(g.edges[i][0]);
    }
    hue = []; hueDim = []; hueNode = [];
    for (var w = 0; w < g.walks.length; w++) {
      var nodes = g.walks[w].nodes;
      for (var p = 0; p < nodes.length; p++) {
        walkOf[nodes[p]].push(w);
        posOf[nodes[p]].push(p);
        if (owner[nodes[p]] < 0) owner[nodes[p]] = w;
      }
      var h = g.walks[w].hue;
      hue.push("hsl(" + h + " 72% 66%)");
      hueDim.push("hsl(" + h + " 30% 30%)");
      hueNode.push("hsl(" + h + " 60% 80%)");
    }
  }

  // Every node on every path through `start`: back to every seed that feeds
  // it, forward to its distinguished point.
  function pathsThrough(start) {
    var nodes = new Set();
    var edges = new Set();
    var stack = [start];
    var i;
    while (stack.length) {
      var n = stack.pop();
      if (nodes.has(n)) continue;
      nodes.add(n);
      for (i = 0; i < preds[n].length; i++) { edges.add(preds[n][i] + "-" + n); stack.push(preds[n][i]); }
    }
    var n2 = start;
    while (succ[n2] >= 0) {
      edges.add(n2 + "-" + succ[n2]);
      n2 = succ[n2];
      nodes.add(n2);
    }
    return { nodes: nodes, edges: edges };
  }

  // Where walk w stands after t iterations: the index of the last node it
  // has passed and the fraction of the way to the next. Nodes sit every
  // `every` iterations except the last, which is wherever the walk ended.
  function standing(w, t) {
    var walk = graph.walks[w];
    var last = walk.nodes.length - 1;
    if (t >= walk.steps) return { i: last, f: 0, done: true };
    var every = graph.every || 1;
    var i = Math.min(last - 1, Math.floor(t / every));
    var span = i + 1 < last ? every : walk.steps - i * every;
    return { i: i, f: Math.max(0, Math.min(1, (t - i * every) / span)), done: false };
  }

  function lerp(a, b, f) {
    var pa = toScreen(a), pb = toScreen(b);
    return [pa[0] + (pb[0] - pa[0]) * f, pa[1] + (pb[1] - pa[1]) * f];
  }

  function drawGuides(w, h) {
    var c = graph.circle;
    if (!c) return;
    var cx = c.rim * view.k + view.x, cy = c.rim * view.k + view.y;
    ctx.lineWidth = 1;
    for (var q = 0; q < c.rings; q++) {
      var r = (c.rim - (c.rim - c.hub) * q / c.rings) * view.k;
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, 2 * Math.PI);
      ctx.strokeStyle = q === 0 ? COLORS.dpDim : COLORS.guide;
      ctx.stroke();
    }
    // Labels only when the rings are far enough apart to hold them.
    if ((c.rim - c.hub) / c.rings * view.k < 76) return;
    ctx.font = "11px ui-monospace, Menlo, monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.lineWidth = 4;
    ctx.strokeStyle = COLORS.bg;
    ctx.fillStyle = COLORS.tick;
    for (var q2 = 1; q2 < c.rings; q2++) {
      var r2 = (c.rim - (c.rim - c.hub) * q2 / c.rings) * view.k;
      var text = Math.floor(c.depth * q2 / c.rings).toLocaleString("en-US") + " to go";
      ctx.strokeText(text, cx + r2, cy);
      ctx.fillText(text, cx + r2, cy);
    }
  }

  function draw(now) {
    needsDraw = false;
    if (!graph) return;
    var w = canvas.width, h = canvas.height;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    var k = view.k;
    var r = Math.max(1.2, Math.min(4, 2.6 * k));
    var rDp = Math.max(1.8, Math.min(6, 3.6 * k));
    var dim = lit !== null || clock !== null;
    var i, a, b, pa, pb, wk;

    drawGuides(w, h);

    // Edges in their walk's colour, faded while a path is lit or the walks
    // are playing; a stride two walks share takes the first one's colour.
    ctx.lineWidth = Math.max(0.6, Math.min(2, 1.1 * k));
    ctx.lineCap = "round";
    for (wk = 0; wk < graph.walks.length; wk++) {
      var nodes = graph.walks[wk].nodes;
      ctx.strokeStyle = dim ? (clock !== null ? hueDim[wk] : COLORS.edgeDim) : hue[wk];
      ctx.beginPath();
      for (i = 0; i + 1 < nodes.length; i++) {
        a = nodes[i]; b = nodes[i + 1];
        if (owner[a] !== wk) continue;
        if (lit !== null && litEdges.has(a + "-" + b)) continue;
        pa = toScreen(a); pb = toScreen(b);
        ctx.moveTo(pa[0], pa[1]); ctx.lineTo(pb[0], pb[1]);
      }
      ctx.stroke();
    }
    if (lit !== null) {
      ctx.lineWidth = Math.max(1.4, Math.min(3.5, 2.4 * k));
      ctx.strokeStyle = COLORS.walkA;
      ctx.beginPath();
      for (i = 0; i < graph.edges.length; i++) {
        a = graph.edges[i][0]; b = graph.edges[i][1];
        if (!litEdges.has(a + "-" + b)) continue;
        pa = toScreen(a); pb = toScreen(b);
        ctx.moveTo(pa[0], pa[1]); ctx.lineTo(pb[0], pb[1]);
      }
      ctx.stroke();
    } else if (clock === null && graph.highlights.length) {
      // The static figure's meetings, in its colours.
      var cls = [["a", COLORS.walkA], ["b", COLORS.walkB], ["shared", COLORS.dp]];
      ctx.lineWidth = Math.max(1.2, Math.min(3, 2.4 * k));
      for (var m = 0; m < graph.highlights.length; m++) {
        for (var c = 0; c < cls.length; c++) {
          var trail = graph.highlights[m][cls[c][0]];
          ctx.strokeStyle = cls[c][1];
          ctx.beginPath();
          for (i = 0; i < trail.length; i++) {
            var p = toScreen(trail[i]);
            if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
          }
          ctx.stroke();
        }
      }
    }

    // Playing: each walk's trail so far, bright in its colour.
    var heads = null;
    if (clock !== null) {
      heads = [];
      ctx.lineWidth = Math.max(1.2, Math.min(3, 2 * k));
      for (wk = 0; wk < graph.walks.length; wk++) {
        var walkNodes = graph.walks[wk].nodes;
        var s = standing(wk, clock);
        ctx.strokeStyle = hue[wk];
        ctx.beginPath();
        var p0 = toScreen(walkNodes[0]);
        ctx.moveTo(p0[0], p0[1]);
        for (i = 1; i <= s.i; i++) {
          var pi = toScreen(walkNodes[i]);
          ctx.lineTo(pi[0], pi[1]);
        }
        var head = s.done ? toScreen(walkNodes[s.i]) : lerp(walkNodes[s.i], walkNodes[s.i + 1], s.f);
        ctx.lineTo(head[0], head[1]);
        ctx.stroke();
        heads.push(s.done ? null : head);
      }
    }

    // Nodes: hollow in their walk's colour, filled green on the rim.
    ctx.lineWidth = Math.max(0.6, Math.min(1.6, 1.1 * k));
    for (i = 0; i < graph.nodes.length; i++) {
      var n = graph.nodes[i];
      var p2 = toScreen(i);
      if (p2[0] < -8 || p2[1] < -8 || p2[0] > w / dpr + 8 || p2[1] > h / dpr + 8) continue;
      var isLit = lit === null ? clock === null : lit.has(i);
      ctx.beginPath();
      if (n[3]) {
        var home = clock !== null && walkOf[i].some(function (x) { return arrived[x] > 0; });
        ctx.arc(p2[0], p2[1], rDp, 0, 2 * Math.PI);
        ctx.fillStyle = isLit || home ? COLORS.dp : COLORS.dpDim;
        ctx.fill();
      } else {
        ctx.arc(p2[0], p2[1], r, 0, 2 * Math.PI);
        ctx.fillStyle = COLORS.node;
        ctx.fill();
        ctx.strokeStyle = isLit ? hueNode[owner[i]] : COLORS.ringDim;
        ctx.stroke();
      }
    }

    // Walkers, and a ring spreading from each distinguished point as its
    // walk arrives.
    if (heads) {
      for (wk = 0; wk < heads.length; wk++) {
        if (heads[wk]) {
          ctx.beginPath();
          ctx.arc(heads[wk][0], heads[wk][1], rDp + 1.5, 0, 2 * Math.PI);
          ctx.fillStyle = hue[wk];
          ctx.fill();
          ctx.beginPath();
          ctx.arc(heads[wk][0], heads[wk][1], Math.max(1.2, rDp - 1), 0, 2 * Math.PI);
          ctx.fillStyle = COLORS.walker;
          ctx.fill();
        } else if (arrived[wk] && now - arrived[wk] < ARRIVE_MS) {
          var age = (now - arrived[wk]) / ARRIVE_MS;
          var nodes2 = graph.walks[wk].nodes;
          var pd = toScreen(nodes2[nodes2.length - 1]);
          ctx.beginPath();
          ctx.arc(pd[0], pd[1], rDp + 2 + 14 * age, 0, 2 * Math.PI);
          ctx.strokeStyle = COLORS.dp;
          ctx.globalAlpha = 1 - age;
          ctx.lineWidth = 2;
          ctx.stroke();
          ctx.globalAlpha = 1;
        }
      }
    }

    // Hovered or selected node, ringed.
    var focus = hover >= 0 ? hover : selected;
    if (focus >= 0) {
      var p4 = toScreen(focus);
      ctx.beginPath();
      ctx.arc(p4[0], p4[1], rDp + 4, 0, 2 * Math.PI);
      ctx.strokeStyle = COLORS.label;
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }
  }

  function nearest(sx, sy) {
    if (!graph) return -1;
    var best = -1, bestD = 12 * 12;
    for (var i = 0; i < graph.nodes.length; i++) {
      var p = toScreen(i);
      var dx = p[0] - sx, dy = p[1] - sy;
      var d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }

  function describe(i) {
    var n = graph.nodes[i];
    var parts = [];
    var every = graph.every || 1;
    for (var j = 0; j < walkOf[i].length; j++) {
      var w = walkOf[i][j];
      var walk = graph.walks[w];
      var pos = posOf[i][j];
      var step = pos < walk.nodes.length - 1 ? pos * every : walk.steps;
      parts.push("walk " + walk.id + ", step " + step.toLocaleString("en-US") + " of " + walk.steps.toLocaleString("en-US"));
    }
    var what = n[3] ? "distinguished point" : (preds[i].length === 0 ? "start" : (preds[i].length > 1 ? "meeting of " + preds[i].length + " trails" : "orbit"));
    return what + " · " + n[0] + (parts.length ? " · " + parts.join(" · ") : "");
  }

  function setReadout(text) {
    if (readout) readout.textContent = text;
  }

  function summary() {
    var c = graph.counts;
    return graph.title + ": " + c.walks + " walks, " + c.nodes.toLocaleString("en-US") + " nodes, " +
      c.distinguished + " distinguished points on the rim, " + c.meetings + " meetings" +
      (graph.every > 1 ? ", one node every " + graph.every.toLocaleString("en-US") + " iterations" : "") +
      ". Distance in from the rim is iterations still to go. Drag to pan, wheel to zoom, hover a node, click one to light its paths.";
  }

  function select_(i) {
    selected = i;
    if (i < 0) { lit = null; litEdges = null; }
    else { var t = pathsThrough(i); lit = t.nodes; litEdges = t.edges; }
    if (clock === null) setReadout(i < 0 ? summary() : describe(i));
    needsDraw = true;
  }

  // --- animation ------------------------------------------------------
  // Every walk sets off from its seed at iteration 0 and moves at the same
  // number of iterations per second, so the short walks reach the rim first
  // and the order they arrive in is the order they finished.
  function depth() {
    return (graph.circle && graph.circle.depth) || graph.counts.longest;
  }
  function play() {
    if (!graph) return;
    select_(-1);
    arrived = [];
    for (var i = 0; i < graph.walks.length; i++) arrived.push(0);
    clock = 0;
    started = performance.now();
    playing = true;
    playButton.textContent = "Stop";
    needsDraw = true;
  }
  function stop() {
    playing = false;
    clock = null;
    arrived = null;
    if (playButton) playButton.textContent = "Play the walks";
    if (graph) setReadout(selected < 0 ? summary() : describe(selected));
    needsDraw = true;
  }
  function tick(now) {
    if (playing) {
      var total = depth();
      clock = Math.min(total, (now - started) / PLAY_MS * total);
      var home = 0;
      for (var i = 0; i < graph.walks.length; i++) {
        if (clock >= graph.walks[i].steps) {
          if (!arrived[i]) arrived[i] = now;
          home++;
        }
      }
      setReadout("iteration " + Math.floor(clock).toLocaleString("en-US") + " of " + total.toLocaleString("en-US") +
        " · " + home + " of " + graph.walks.length + " walks at their distinguished point");
      needsDraw = true;
      if (clock >= total && now - started > PLAY_MS + ARRIVE_MS) {
        playing = false;
        playButton.textContent = "Play again";
      }
    }
    if (needsDraw) draw(now);
    window.requestAnimationFrame(tick);
  }

  // --- events ---------------------------------------------------------
  var dragging = null;
  canvas.addEventListener("pointerdown", function (e) {
    dragging = { x: e.clientX, y: e.clientY, vx: view.x, vy: view.y, moved: false };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", function (e) {
    var rect = canvas.getBoundingClientRect();
    if (dragging) {
      var dx = e.clientX - dragging.x, dy = e.clientY - dragging.y;
      if (Math.abs(dx) + Math.abs(dy) > 3) dragging.moved = true;
      view.x = dragging.vx + dx;
      view.y = dragging.vy + dy;
      needsDraw = true;
      return;
    }
    var i = nearest(e.clientX - rect.left, e.clientY - rect.top);
    if (i !== hover) {
      hover = i;
      needsDraw = true;
      if (i >= 0) {
        tip.textContent = describe(i);
        tip.hidden = false;
      } else {
        tip.hidden = true;
      }
    }
    if (i >= 0) {
      tip.style.left = (e.clientX - rect.left + 12) + "px";
      tip.style.top = (e.clientY - rect.top + 12) + "px";
    }
  });
  canvas.addEventListener("pointerup", function (e) {
    if (dragging && !dragging.moved) {
      if (clock !== null && !playing) stop();
      var rect = canvas.getBoundingClientRect();
      var i = nearest(e.clientX - rect.left, e.clientY - rect.top);
      select_(i === selected ? -1 : i);
    }
    dragging = null;
  });
  canvas.addEventListener("pointerleave", function () { hover = -1; tip.hidden = true; needsDraw = true; });
  canvas.addEventListener("wheel", function (e) {
    e.preventDefault();
    var rect = canvas.getBoundingClientRect();
    var sx = e.clientX - rect.left, sy = e.clientY - rect.top;
    var factor = Math.exp(-e.deltaY * 0.0015);
    var k2 = Math.max(0.2, Math.min(12, view.k * factor));
    factor = k2 / view.k;
    view.x = sx - (sx - view.x) * factor;
    view.y = sy - (sy - view.y) * factor;
    view.k = k2;
    needsDraw = true;
  }, { passive: false });
  // Pinch: two pointers.
  var pinch = null;
  canvas.addEventListener("touchstart", function (e) {
    if (e.touches.length === 2) {
      pinch = { d: Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY), k: view.k };
    }
  }, { passive: true });
  canvas.addEventListener("touchmove", function (e) {
    if (pinch && e.touches.length === 2) {
      e.preventDefault();
      var d = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
      var rect = canvas.getBoundingClientRect();
      var cx = (e.touches[0].clientX + e.touches[1].clientX) / 2 - rect.left;
      var cy = (e.touches[0].clientY + e.touches[1].clientY) / 2 - rect.top;
      var k2 = Math.max(0.2, Math.min(12, pinch.k * d / pinch.d));
      var factor = k2 / view.k;
      view.x = cx - (cx - view.x) * factor;
      view.y = cy - (cy - view.y) * factor;
      view.k = k2;
      needsDraw = true;
    }
  }, { passive: false });
  canvas.addEventListener("touchend", function () { pinch = null; }, { passive: true });

  if (playButton) playButton.addEventListener("click", function () { if (playing) stop(); else play(); });
  if (resetButton) resetButton.addEventListener("click", function () { select_(-1); stop(); fit(); });
  window.addEventListener("resize", function () { size(); fit(); });

  function load(dataset) {
    setReadout("Loading " + dataset.label + "…");
    return fetch(dataset.file, { cache: "no-store" })
      .then(function (response) {
        if (!response.ok) throw new Error("HTTP " + response.status);
        return response.json();
      })
      .then(function (g) {
        graph = g;
        index(g);
        stop();
        select_(-1);
        fit();
        if (image && image.parentNode === host) {
          host.insertBefore(canvas, image);
          host.appendChild(tip);
          host.removeChild(image);
        }
        if (controls) controls.hidden = false;
        setReadout(summary());
      })
      .catch(function (err) {
        setReadout("The explorer could not load " + dataset.file + " (" + err.message + "); the figure above is the same forest.");
      });
  }

  if (select) {
    for (var d = 0; d < DATASETS.length; d++) {
      var option = document.createElement("option");
      option.value = DATASETS[d].key;
      option.textContent = DATASETS[d].label;
      select.appendChild(option);
    }
    select.addEventListener("change", function () {
      for (var d2 = 0; d2 < DATASETS.length; d2++) {
        if (DATASETS[d2].key === select.value) load(DATASETS[d2]);
      }
    });
  }

  size();
  host.style.position = "relative";
  load(DATASETS[0]);
  window.requestAnimationFrame(tick);
})();
