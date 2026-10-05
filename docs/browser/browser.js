/* Lab browser: a searchable index over docs/browser/data.json.
   Every node is built with createElement/textContent; no HTML strings. */
(function () {
  "use strict";
  var REPO = "https://github.com/aburan28/crypto/blob/main/";
  var DATA = null;
  var byId = { curve: {}, method: {}, fb: {}, candidate: {}, round: {}, session: {} };
  var view = document.getElementById("view");
  var loading = document.getElementById("loading");

  function el(tag, attrs, children) {
    var node = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        if (k === "text") node.textContent = attrs[k];
        else if (k === "class") node.className = attrs[k];
        else if (attrs[k] !== null && attrs[k] !== undefined) node.setAttribute(k, attrs[k]);
      });
    }
    (children || []).forEach(function (c) {
      if (c === null || c === undefined) return;
      node.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
    });
    return node;
  }
  function link(href, text, cls) { return el("a", { href: href, class: cls || null, text: text }); }
  function repoLink(path, text) { return link(REPO + path, text || path); }
  function chip(text, cls) { return el("span", { class: "chip" + (cls ? " " + cls : ""), text: text }); }
  function fmt(v, digits) {
    if (v === null || v === undefined || v === "") return "–";
    if (typeof v === "number") return Number.isInteger(v) ? String(v) : v.toFixed(digits === undefined ? 3 : digits);
    return String(v);
  }
  /* goodAbove: for rho / winner, a value above one means the candidate was
     cheaper than that round's rho; for IC / rho it is the reverse. */
  function ratioCell(v, goodAbove) {
    if (v === null || v === undefined) return el("td", { class: "n", text: "–" });
    var good = goodAbove ? v > 1 : v < 1;
    return el("td", { class: "n" }, [el("span", { class: "ratio" + (good ? " below" : ""), text: fmt(v, 2) + "×" })]);
  }
  function fact(label, value, mono) {
    var dd = el("dd", { class: mono ? "mono" : null });
    if (value instanceof Node) dd.appendChild(value); else dd.textContent = fmt(value);
    return el("div", { class: "fact" }, [el("dt", { text: label }), dd]);
  }
  function fieldLabel(c) {
    return c.field.characteristic === 2 ? "GF(2^" + c.field.degree + ")" : "GF(p), p of " + c.field.bits + " bits";
  }
  function bestRatio(c) {
    var vals = (c.leaderboard || []).map(function (b) { return b.ratio_rho; }).filter(function (v) { return typeof v === "number"; });
    return vals.length ? Math.min.apply(null, vals) : null;
  }
  function parseHash() {
    var h = location.hash.replace(/^#/, "") || "curves";
    var q = {};
    var qi = h.indexOf("?");
    if (qi >= 0) {
      h.slice(qi + 1).split("&").forEach(function (kv) {
        if (!kv) return;
        var p = kv.split("=");
        q[decodeURIComponent(p[0])] = decodeURIComponent(p.slice(1).join("=") || "");
      });
      h = h.slice(0, qi);
    }
    var parts = h.split("/");
    return { view: parts[0], id: parts.slice(1).join("/"), q: q };
  }
  function setQuery(viewName, q) {
    var pairs = Object.keys(q).filter(function (k) { return q[k] !== "" && q[k] !== null && q[k] !== undefined && q[k] !== false; })
      .map(function (k) { return encodeURIComponent(k) + "=" + encodeURIComponent(q[k]); });
    history.replaceState(null, "", "#" + viewName + (pairs.length ? "?" + pairs.join("&") : ""));
  }

  /* ---------- tables ---------- */
  function table(columns, rows, opts) {
    opts = opts || {};
    var state = { key: opts.sortKey || columns[0].key, dir: opts.sortDir || 1 };
    var wrap = el("div", { class: "scroll" });
    var tbl = el("table");
    var thead = el("thead");
    var tr = el("tr");
    columns.forEach(function (col) {
      var th = el("th", { class: col.numeric ? "n" : null, text: col.label, scope: "col" });
      th.addEventListener("click", function () {
        if (state.key === col.key) state.dir = -state.dir; else { state.key = col.key; state.dir = col.numeric ? -1 : 1; }
        render();
      });
      tr.appendChild(th);
    });
    thead.appendChild(tr);
    tbl.appendChild(thead);
    var tbody = el("tbody");
    tbl.appendChild(tbody);
    wrap.appendChild(tbl);
    function render() {
      var col = columns.filter(function (c) { return c.key === state.key; })[0];
      var sorted = rows.slice().sort(function (a, b) {
        var x = col.sort ? col.sort(a) : a[col.key], y = col.sort ? col.sort(b) : b[col.key];
        var xn = x === null || x === undefined, yn = y === null || y === undefined;
        if (xn && yn) return 0; if (xn) return 1; if (yn) return -1;
        if (typeof x === "number" && typeof y === "number") return (x - y) * state.dir;
        return String(x).localeCompare(String(y), undefined, { numeric: true }) * state.dir;
      });
      Array.prototype.forEach.call(thead.querySelectorAll("th"), function (th, i) {
        th.setAttribute("aria-sort", columns[i].key === state.key ? (state.dir > 0 ? "ascending" : "descending") : "none");
      });
      while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
      if (!sorted.length) {
        tbody.appendChild(el("tr", null, [el("td", { colspan: String(columns.length), class: "empty", text: opts.empty || "Nothing matches." })]));
      }
      sorted.forEach(function (row) {
        var r = el("tr");
        columns.forEach(function (c) {
          var cell = c.cell ? c.cell(row) : el("td", { class: c.numeric ? "n" : null, text: fmt(row[c.key], c.digits) });
          r.appendChild(cell);
        });
        tbody.appendChild(r);
      });
    }
    render();
    return wrap;
  }
  function tdLink(href, text) { return el("td", null, [link(href, text)]); }
  function tdText(text, cls) { return el("td", { class: cls || null, text: text }); }
  function tdList(items, cls) {
    var td = el("td", { class: cls || "wrap" });
    items.forEach(function (it, i) { if (i) td.appendChild(document.createTextNode(", ")); td.appendChild(it); });
    if (!items.length) td.textContent = "–";
    return td;
  }

  /* ---------- curves ---------- */
  function curvesView(q) {
    var head = el("div", null, [
      el("h1", { text: "Curves" }),
      el("p", { class: "lede", text: "Every curve in the ICV1 registry: " + DATA.counts.curves + " models. Search by name or identity, filter by family, field, trace, order, subgroup and what has been measured on it. Click a slug for its full record and everything that cites it. A conductor is not an invariant of a curve over a finite field; the endomorphism discriminant (−7 for every Koblitz curve) is the nearest thing recorded." })
    ]);
    var f = el("form", { class: "filters" });
    function input(name, label, type, extra) {
      var inp = el("input", { type: type || "text", name: name, value: q[name] || "", placeholder: extra && extra.placeholder || null, inputmode: type === "number" ? "numeric" : null });
      var lab = el("label", null, [label, inp]);
      if (extra && extra.wide) lab.className = "wide";
      f.appendChild(lab);
      return inp;
    }
    function range(name, label) {
      var lo = el("input", { type: "number", name: name + "_min", value: q[name + "_min"] || "", placeholder: "min" });
      var hi = el("input", { type: "number", name: name + "_max", value: q[name + "_max"] || "", placeholder: "max" });
      f.appendChild(el("label", null, [label, el("div", { class: "range" }, [lo, hi])]));
      return [lo, hi];
    }
    function select(name, label, options) {
      var sel = el("select", { name: name });
      options.forEach(function (o) {
        var opt = el("option", { value: o[0], text: o[1] });
        if ((q[name] || "") === o[0]) opt.selected = true;
        sel.appendChild(opt);
      });
      f.appendChild(el("label", null, [label, sel]));
      return sel;
    }
    function check(name, label) {
      var inp = el("input", { type: "checkbox", name: name });
      inp.checked = q[name] === "1";
      f.appendChild(el("label", { class: "check" }, [inp, label]));
      return inp;
    }
    var text = input("q", "Search slug, ICV1, EC1, standard or retired name", "search", { wide: true, placeholder: "e.g. ECC2K-130, f2m41, sect163k1, k0n31, EC1N7C…" });
    var family = select("family", "Family", [["", "any"], ["koblitz", "Koblitz"], ["binary", "random binary"], ["prime", "prime field"], ["subfield", "subfield"]]);
    var ch = select("char", "Characteristic", [["", "any"], ["2", "2"], ["p", "p (prime field)"]]);
    var deg = range("bits", "Field bits (degree m, or bits of p)");
    var tr = range("trace", "Trace of Frobenius");
    var ob = range("order_bits", "log₂ #E");
    var rb = range("r_bits", "log₂ r (prime subgroup)");
    var cof = input("cofactor", "Cofactor", "number");
    var a = select("a", "Curve a", [["", "any"], ["0", "0"], ["1", "1"]]);
    var endo = select("end", "Endomorphism discriminant", [["", "any"], ["-7", "−7 (Koblitz)"], ["unk", "not recorded"]]);
    var m1 = check("ecbench", "measured in an ecbench session");
    var m2 = check("board", "on the leaderboard");
    var m3 = check("tournament", "a tournament cell");
    var m4 = check("standard", "a published standard or challenge");
    var m5 = check("ec1", "has an EC1 representation");
    var reset = el("button", { type: "button", text: "Clear filters" });
    f.appendChild(el("label", null, ["", reset]));
    var count = el("p", { class: "count" });
    var holder = el("div");
    function current() {
      var out = {};
      Array.prototype.forEach.call(f.querySelectorAll("input,select"), function (i) {
        if (i.type === "checkbox") { if (i.checked) out[i.name] = "1"; }
        else if (i.value !== "") out[i.name] = i.value;
      });
      return out;
    }
    function num(v) { return v === "" || v === undefined ? null : Number(v); }
    function inRange(v, lo, hi) {
      if (v === null || v === undefined) return lo === null && hi === null;
      return (lo === null || v >= lo) && (hi === null || v <= hi);
    }
    function apply() {
      var s = current();
      setQuery("curves", s);
      var needle = (s.q || "").toLowerCase().trim();
      var rows = DATA.curves.filter(function (c) {
        if (s.family && c.family !== s.family) return false;
        if (s["char"] && String(c.field.characteristic) !== s["char"]) return false;
        if (!inRange(c.field.bits, num(s.bits_min), num(s.bits_max))) return false;
        if (!inRange(c.trace === null ? null : Number(c.trace), num(s.trace_min), num(s.trace_max))) return false;
        if (!inRange(c.order_bits, num(s.order_bits_min), num(s.order_bits_max))) return false;
        if ((s.r_bits_min || s.r_bits_max) && !inRange(c.r_bits, num(s.r_bits_min), num(s.r_bits_max))) return false;
        if (s.cofactor && String(c.cofactor) !== s.cofactor) return false;
        if (s.a && String(c.a) !== s.a) return false;
        if (s.end === "-7" && String(c.endomorphism_discriminant) !== "-7") return false;
        if (s.end === "unk" && c.endomorphism_discriminant !== null) return false;
        if (s.ecbench && !c.ecbench.length) return false;
        if (s.board && !c.on_leaderboard) return false;
        if (s.tournament && !c.tournament_cells.length) return false;
        if (s.standard && !c.standard_names.length) return false;
        if (s.ec1 && !c.ec1.length) return false;
        if (needle) {
          var hay = [c.slug, c.icv1, c.construction].concat(c.standard_names, c.retired_names, c.ec1, c.curve_uid).join(" ").toLowerCase();
          if (hay.indexOf(needle) < 0) return false;
        }
        return true;
      });
      count.textContent = rows.length + " of " + DATA.curves.length + " curves";
      while (holder.firstChild) holder.removeChild(holder.firstChild);
      holder.appendChild(table([
        { key: "slug", label: "Curve", cell: function (c) { return el("td", null, [link("#curve/" + c.slug, c.slug, "mono")].concat(c.standard_names.map(function (n) { return chip(n, "data"); }))); } },
        { key: "family", label: "Family" },
        { key: "field_bits", label: "Field", numeric: true, sort: function (c) { return c.field.bits; }, cell: function (c) { return tdText(fieldLabel(c), "n"); } },
        { key: "a", label: "a", cell: function (c) { return tdText(fmt(c.a)); } },
        { key: "trace", label: "Trace", numeric: true, sort: function (c) { return Number(c.trace); }, cell: function (c) { return tdText(fmt(c.trace), "n"); } },
        { key: "order_bits", label: "log₂ #E", numeric: true, digits: 2 },
        { key: "r_bits", label: "log₂ r", numeric: true, digits: 2 },
        { key: "cofactor", label: "Cofactor", numeric: true, sort: function (c) { return c.cofactor === null ? null : Number(c.cofactor); } },
        { key: "end", label: "End. disc.", cell: function (c) { return tdText(c.endomorphism_discriminant === null ? "–" : String(c.endomorphism_discriminant), "n"); } },
        { key: "cover", label: "Cover", cell: function (c) {
          var cert = c.hyperelliptic_cover && c.hyperelliptic_cover.certificate;
          return tdText(cert ? "genus " + cert.genus + ", degree " + cert.degree : "unresolved");
        } },
        { key: "best", label: "Best IC / rho", numeric: true, sort: bestRatio, cell: function (c) { return ratioCell(bestRatio(c)); } },
        { key: "measured", label: "Cited by", sort: function (c) { return c.ecbench.length + c.tournament_cells.length + (c.on_leaderboard ? 1 : 0); }, cell: function (c) {
          var items = [];
          if (c.on_leaderboard) items.push(chip("leaderboard", "ok"));
          if (c.ecbench.length) items.push(chip("ecbench ×" + new Set(c.ecbench.map(function (e) { return e.session_id; })).size, "ok"));
          if (c.tournament_cells.length) items.push(chip("tournament ×" + new Set(c.tournament_cells.map(function (t) { return t.round; })).size, "ok"));
          if (c.factor_bases.length) items.push(chip("factor bases " + c.factor_bases.length));
          if (!items.length) items.push(chip("registered only"));
          var td = el("td"); items.forEach(function (i) { td.appendChild(i); }); return td;
        } }
      ], rows, { sortKey: "field_bits", sortDir: 1 }));
    }
    f.addEventListener("input", apply);
    f.addEventListener("submit", function (e) { e.preventDefault(); apply(); });
    reset.addEventListener("click", function () {
      Array.prototype.forEach.call(f.querySelectorAll("input,select"), function (i) { if (i.type === "checkbox") i.checked = false; else i.value = ""; });
      apply();
    });
    view.appendChild(head); view.appendChild(f); view.appendChild(count); view.appendChild(holder);
    apply();
  }

  function curveDetail(slug) {
    var c = byId.curve[slug];
    if (!c) { view.appendChild(el("p", { class: "note", text: "No curve named " + slug + " in the registry." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#curves", "Curves"), " / " + c.slug]));
    var h = el("h1", { text: c.slug });
    view.appendChild(h);
    var chips = el("p");
    chips.appendChild(chip(c.family, "data"));
    c.standard_names.forEach(function (n) { chips.appendChild(chip(n, "data")); });
    if (c.on_leaderboard) chips.appendChild(chip("on the leaderboard", "ok"));
    if (c.ec1_unresolved) chips.appendChild(chip("no EC1 representation", "warn"));
    view.appendChild(chips);
    var facts = el("dl", { class: "facts" }, [
      fact("ICV1 identity", c.icv1, true),
      fact("Field", c.field.characteristic === 2 ? "GF(2^" + c.field.degree + "), modulus " + fmt(c.field.modulus) : "GF(p), p = " + c.field.p + " (" + c.field.bits + " bits)", false),
      fact("Coefficients", "a = " + fmt(c.a) + ", b = " + fmt(c.b), true),
      fact("Trace of Frobenius", c.trace),
      fact("Group order #E", c.order, true),
      fact("log₂ #E", fmt(c.order_bits, 3)),
      fact("Prime subgroup order r", c.r, true),
      fact("log₂ r", fmt(c.r_bits, 3)),
      fact("Cofactor", c.cofactor),
      fact("Generator", c.generator ? c.generator.join(", ") : null, true),
      fact("j-invariant", c.j, true),
      fact("Endomorphism discriminant", c.endomorphism_discriminant === null ? "not recorded" : c.endomorphism_discriminant),
      fact("Construction", c.construction, true)
    ]);
    view.appendChild(facts);

    view.appendChild(el("h2", { text: "Hyperelliptic cover" }));
    var cover = c.hyperelliptic_cover;
    if (cover) {
      var cert = cover.certificate;
      view.appendChild(el("dl", { class: "facts" }, [
        fact("Existence", cover.exists === true ? "Verified over the declared field" : (cover.status === "invalid_input" ? "Invalid input model" : "Unresolved")),
        fact("Cover direction", "H → E"),
        fact("Cover genus", cert ? cert.genus : null),
        fact("Map degree", cert ? cert.degree : null),
        fact("Field", cert ? "Same field as E" : "Not established"),
        fact("Field validation", cover.field_check || cover.reason),
        fact("Minimum genus / degree", "Not determined"),
        fact("Subfield descent / subgroup transfer", "Not tested"),
        fact("DLP advantage", "Not measured")
      ]));
      view.appendChild(el("p", { class: "note" }, [cert ? "An explicit geometric cover is recorded. Computational advantage requires separate evidence. " : "This result establishes neither existence nor nonexistence. ", repoLink("docs/curves/COVERS.md", "Proof and scope"), "; ", repoLink("docs/curves/covers.json", "replayable catalog findings"), "."]));
      if (cert) {
        var details = el("details", null, [el("summary", { text: "Cover equation and map certificate" })]);
        details.appendChild(el("p", { text: "Ascending powers of u; H: v² + h(u)v = f(u). Map: x = x(u), y = y_v(u)v + y_0(u). Coefficients use the curve's field representation." }));
        details.appendChild(el("pre", { text: JSON.stringify(cert, null, 2) }));
        view.appendChild(details);
      }
    } else {
      view.appendChild(el("p", { class: "note", text: "Cover existence has not been checked for this model." }));
    }

    view.appendChild(el("h2", { text: "Identities" }));
    var idl = el("dl", { class: "facts" });
    idl.appendChild(fact("EC1 alias", c.ec1.length ? c.ec1.join(", ") : (c.ec1_unresolved || "none"), true));
    idl.appendChild(fact("Curve UID", c.curve_uid.length ? c.curve_uid.join(", ") : "none", true));
    idl.appendChild(fact("Retired names (never write these)", c.retired_names.length ? c.retired_names.join(", ") : "none", true));
    idl.appendChild(fact("Registry sources", el("span", null, c.sources.map(function (s, i) {
      var isPath = /^docs\/|^research\/|^src\//.test(s);
      return el("span", null, [i ? "; " : "", isPath ? repoLink(s) : s]);
    }))));
    view.appendChild(idl);
    view.appendChild(el("h2", { text: "Model (the hashed record)" }));
    var model = c.model_json;
    try { model = JSON.stringify(JSON.parse(c.model_json), null, 1); } catch (e) { /* keep raw */ }
    view.appendChild(el("pre", { text: model }));

    view.appendChild(el("h2", { text: "On the leaderboard" }));
    if (c.leaderboard.length) {
      view.appendChild(table([
        { key: "regime", label: "Regime" },
        { key: "log2_r", label: "log₂ r", numeric: true, digits: 2 },
        { key: "best_recipe", label: "Best recipe", cell: function (b) { return tdText(fmt(b.best_recipe), "mono"); } },
        { key: "best_s", label: "S", numeric: true, digits: 3 },
        { key: "ratio_rho", label: "IC / rho", numeric: true, cell: function (b) { return ratioCell(b.ratio_rho); } },
        { key: "ratio_floor", label: "IC / floor", numeric: true, cell: function (b) { return ratioCell(b.ratio_floor); } },
        { key: "reference", label: "Reference", cell: function (b) { return tdText(fmt(b.reference), "wrap"); } }
      ], c.leaderboard));
      view.appendChild(el("p", { class: "note" }, ["Quoted from ", repoLink("docs/ic/leaderboard.json"), "; drawn on the ", link("../scoreboard/ic-leaderboard.html", "leaderboard page"), "."]));
    } else view.appendChild(el("p", { class: "note", text: "No whole-pipeline index-calculus row on the leaderboard for this curve." }));

    view.appendChild(el("h2", { text: "Measured by ecbench" }));
    if (c.ecbench.length) {
      view.appendChild(table([
        { key: "session_id", label: "Session", cell: function (e) { return tdLink("#session/" + e.session_id, e.session_id); } },
        { key: "arm", label: "Arm" },
        { key: "method", label: "Method", cell: function (e) { return tdLink("#method/" + e.method_id, e.method); } },
        { key: "verified", label: "Verified / runs", numeric: true, cell: function (e) { return tdText(e.verified + " / " + e.runs, "n"); } },
        { key: "mean_s", label: "Mean S", numeric: true, digits: 3 },
        { key: "mean_ratio_to_floor", label: "S / floor", numeric: true, cell: function (e) { return ratioCell(e.mean_ratio_to_floor); } }
      ], c.ecbench, { sortKey: "mean_s", sortDir: 1 }));
      view.appendChild(el("p", { class: "note", text: "Mean S over the session's own verified, non-warm-up runs on this curve; the session directory holds every record. S / floor compares to the generic floor √(π/2A) the harness records per run." }));
    } else view.appendChild(el("p", { class: "note", text: "No committed ecbench session measured this curve." }));

    view.appendChild(el("h2", { text: "Factor bases built on it" }));
    if (c.factor_bases.length) {
      view.appendChild(table(fbColumns(), c.factor_bases.map(function (id) { return byId.fb[id]; }).filter(Boolean)));
    } else view.appendChild(el("p", { class: "note", text: "No ecbench factor base recorded on this curve." }));

    view.appendChild(el("h2", { text: "Tournament cells" }));
    if (c.tournament_cells.length) {
      view.appendChild(table([
        { key: "round", label: "Round", cell: function (t) { return tdLink("#round/" + t.round, t.round); } },
        { key: "cell", label: "Cell", cell: function (t) { return tdText(t.cell, "mono"); } },
        { key: "stages", label: "Stages", cell: function (t) { return tdText(t.stages.join(", "), "wrap"); } },
        { key: "verdict", label: "Round verdict", sort: function (t) { return (byId.round[t.round] || {}).status; }, cell: function (t) { var r = byId.round[t.round] || {}; return tdText((r.status || "–") + (r.winner ? " · " + r.winner : "")); } },
        { key: "rho", label: "rho / winner", numeric: true, sort: function (t) { return (byId.round[t.round] || {}).rho_over_winner; }, cell: function (t) { return ratioCell((byId.round[t.round] || {}).rho_over_winner, true); } }
      ], c.tournament_cells));
    } else view.appendChild(el("p", { class: "note", text: "Not a cell in any tournament round." }));

    var cands = DATA.candidates.filter(function (k) { return k.slug === c.slug; });
    view.appendChild(el("h2", { text: "Candidate identities naming this curve (" + cands.length + ")" }));
    if (cands.length) view.appendChild(table(candidateColumns(), cands)); else view.appendChild(el("p", { class: "note", text: "No IC1 candidate identity in the repository names this curve." }));
  }

  /* ---------- methods ---------- */
  function paramsText(p) { return Object.keys(p || {}).map(function (k) { return k + "=" + p[k]; }).join(" ") || "(defaults)"; }
  function methodsView() {
    view.appendChild(el("h1", { text: "Methods" }));
    view.appendChild(el("p", { class: "lede", text: "Every ecbench method identity a committed session ran. The id ECM1h… is the first twelve hex digits of the SHA-256 of {schema, id, params} with defaults written out, so two arms with the same id ran the same algorithm with the same parameters, whatever they were called in the spec." }));
    view.appendChild(table([
      { key: "method_id", label: "Method id", cell: function (m) { return tdLink("#method/" + m.method_id, m.method_id); } },
      { key: "method", label: "Method" },
      { key: "family", label: "Family" },
      { key: "params", label: "Parameters", cell: function (m) { return tdText(paramsText(m.params), "mono wrap"); } },
      { key: "verified", label: "Verified / runs", numeric: true, cell: function (m) { return tdText(m.verified + " / " + m.runs, "n"); } },
      { key: "sessions", label: "Sessions", numeric: true, sort: function (m) { return m.sessions.length; }, cell: function (m) { return tdText(String(m.sessions.length), "n"); } },
      { key: "curves", label: "Curves", numeric: true, sort: function (m) { return m.curves.length; }, cell: function (m) { return tdText(String(m.curves.length), "n"); } }
    ], DATA.methods, { sortKey: "family" }));
  }
  function methodDetail(id) {
    var m = byId.method[id];
    if (!m) { view.appendChild(el("p", { class: "note", text: "No method " + id + " in any committed session." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#methods", "Methods"), " / " + m.method_id]));
    view.appendChild(el("h1", { text: m.method + " · " + m.method_id }));
    view.appendChild(el("dl", { class: "facts" }, [
      fact("Family", m.family), fact("Parameters", paramsText(m.params), true), fact("SHA-256 of the method record", m.method_sha256, true),
      fact("Runs (verified / all)", m.verified + " / " + m.runs), fact("Sessions", m.sessions.length), fact("Curves", m.curves.length)
    ]));
    view.appendChild(el("h2", { text: "By curve" }));
    var rows = [];
    DATA.curves.forEach(function (c) {
      c.ecbench.forEach(function (e) { if (e.method_id === id) rows.push({ slug: c.slug, r_bits: c.r_bits, session_id: e.session_id, arm: e.arm, runs: e.runs, verified: e.verified, mean_s: e.mean_s, floor: e.mean_ratio_to_floor }); });
    });
    view.appendChild(table([
      { key: "slug", label: "Curve", cell: function (r) { return tdLink("#curve/" + r.slug, r.slug); } },
      { key: "r_bits", label: "log₂ r", numeric: true, digits: 2 },
      { key: "session_id", label: "Session", cell: function (r) { return tdLink("#session/" + r.session_id, r.session_id); } },
      { key: "arm", label: "Arm" },
      { key: "verified", label: "Verified / runs", numeric: true, cell: function (r) { return tdText(r.verified + " / " + r.runs, "n"); } },
      { key: "mean_s", label: "Mean S", numeric: true, digits: 3 },
      { key: "floor", label: "S / floor", numeric: true, cell: function (r) { return ratioCell(r.floor); } }
    ], rows, { sortKey: "r_bits" }));
    view.appendChild(el("h2", { text: "Sessions" }));
    view.appendChild(el("ul", { class: "plain" }, m.sessions.map(function (s) { var ss = byId.session[s]; return el("li", null, [link("#session/" + s, s), ss ? " · " + ss.label : ""]); })));
  }

  /* ---------- factor bases ---------- */
  function fbColumns() {
    return [
      { key: "fb_id", label: "Factor base", cell: function (f) { return tdLink("#fb/" + f.fb_id, f.fb_id); } },
      { key: "family", label: "Family" },
      { key: "params", label: "Parameters", cell: function (f) { return tdText(paramsText(f.params), "mono"); } },
      { key: "curve", label: "Curve", cell: function (f) { return f.curve ? tdLink("#curve/" + f.curve, f.curve) : tdText("–"); } },
      { key: "signed_points", label: "Signed points", numeric: true },
      { key: "abscissae", label: "Abscissae", numeric: true },
      { key: "columns", label: "Columns", numeric: true },
      { key: "dimension", label: "Dim.", numeric: true },
      { key: "points_sha256", label: "Points SHA-256", cell: function (f) { return tdText((f.points_sha256 || "").slice(0, 16) + "…", "mono"); } }
    ];
  }
  function fbView() {
    view.appendChild(el("h1", { text: "Factor bases" }));
    view.appendChild(el("p", { class: "lede", text: "Every factor base a committed ecbench session built, identified by FB1h…, the hash of its curve, family, parameters, column count and the digest of its enumerated point set. The point digest is what makes two bases the same base; a count alone does not." }));
    view.appendChild(table(fbColumns(), DATA.factor_bases, { sortKey: "curve" }));
    view.appendChild(el("p", { class: "note" }, ["Tournament and ledger runs describe their bases as specs (kind, parent orbit, retained orbits) in their own frozen files, for example ", repoLink("docs/ic/params/"), "; those are not FB1 identities and are reached through each round."]));
  }
  function fbDetail(id) {
    var f = byId.fb[id];
    if (!f) { view.appendChild(el("p", { class: "note", text: "No factor base " + id + "." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#factor-bases", "Factor bases"), " / " + f.fb_id]));
    view.appendChild(el("h1", { text: f.fb_id }));
    view.appendChild(el("dl", { class: "facts" }, [
      fact("Family", f.family), fact("Parameters", paramsText(f.params), true), fact("Description", f.description),
      fact("Curve", f.curve ? link("#curve/" + f.curve, f.curve) : null, true),
      fact("Signed points", f.signed_points), fact("Abscissae", f.abscissae), fact("Columns", f.columns), fact("Dimension", f.dimension),
      fact("Points SHA-256", f.points_sha256, true), fact("Record SHA-256", f.fb_sha256, true)
    ]));
    view.appendChild(el("h2", { text: "Used by" }));
    view.appendChild(el("ul", { class: "plain" }, f.methods.map(function (m) { var mm = byId.method[m]; return el("li", null, [link("#method/" + m, m), mm ? " · " + mm.method + " " + paramsText(mm.params) : ""]); })
      .concat(f.sessions.map(function (s) { var ss = byId.session[s]; return el("li", null, [link("#session/" + s, s), ss ? " · " + ss.label : ""]); }))));
  }

  /* ---------- candidates ---------- */
  function candidateColumns() {
    return [
      { key: "candidate_id", label: "IC1 identity", cell: function (k) { return tdLink("#candidate/" + k.candidate_id, k.candidate_id); } },
      { key: "n", label: "n", numeric: true },
      { key: "slug", label: "Curve", cell: function (k) { return k.slug ? tdLink("#curve/" + k.slug, k.slug) : tdText(k.curve_tag ? "n" + k.n + " " + k.curve_tag + " (not in registry)" : "–"); } },
      { key: "factor_base_points", label: "Base points", numeric: true },
      { key: "summands", label: "Summands", numeric: true },
      { key: "solver", label: "Solver" },
      { key: "collector", label: "Collector" },
      { key: "linear_algebra", label: "Lin. alg." },
      { key: "descent", label: "Descent" },
      { key: "mentions", label: "Mentions", numeric: true },
      { key: "files", label: "Files", numeric: true }
    ];
  }
  function candidatesView() {
    view.appendChild(el("h1", { text: "Candidate identities" }));
    view.appendChild(el("p", { class: "lede", text: "Every tournament candidate identity (IC1N…h…) written anywhere in docs/ic, the tournament directory or the ECC2K-130 notes: " + DATA.counts.candidates + ". The code names the curve, factor-base size, decomposition summands and solver, relation collector, linear algebra, descent and isogeny step; the suffix is the first twelve hex digits of the candidate record's SHA-256. Mentions count occurrences; a candidate measured in a round appears in that round's files." }));
    view.appendChild(table(candidateColumns(), DATA.candidates, { sortKey: "n" }));
  }
  function candidateDetail(id) {
    var k = byId.candidate[id];
    if (!k) { view.appendChild(el("p", { class: "note", text: "No candidate " + id + "." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#candidates", "Candidates"), " / " + k.candidate_id]));
    view.appendChild(el("h1", { text: k.candidate_id }));
    view.appendChild(el("dl", { class: "facts" }, [
      fact("Curve", k.slug ? link("#curve/" + k.slug, k.slug) : "n" + k.n + " " + k.curve_tag + " (not in the registry)", true),
      fact("Factor-base points", k.factor_base_points), fact("Decomposition summands", k.summands), fact("Decomposition solver", k.solver),
      fact("Relation collector", k.collector), fact("Linear algebra", k.linear_algebra), fact("Target descent", k.descent), fact("Isogeny step", k.isogeny),
      fact("Record SHA-256 prefix", k.record_sha12, true), fact("Mentions / files", k.mentions + " / " + k.files)
    ]));
    view.appendChild(el("h2", { text: "Where it is written" }));
    view.appendChild(el("ul", { class: "plain" }, k.where.map(function (p) { return el("li", null, [repoLink(p)]); })));
  }

  /* ---------- rounds ---------- */
  function roundsView() {
    view.appendChild(el("h1", { text: "Tournament rounds" }));
    view.appendChild(el("p", { class: "lede", text: "Every directory under the candidate tournament's runs, with its decision as the round's own files state it. rho / winner is the round's observed ratio in its unit (an instruction count); above one read as the candidate cheaper than that round's rho. The reference changed across rounds, which the scoreboard's reference history documents; rounds before 23 September measured against a rho later shown unmatched." }));
    view.appendChild(table([
      { key: "round", label: "Round", cell: function (r) { return tdLink("#round/" + r.round, r.round); } },
      { key: "status", label: "Status", cell: function (r) { return el("td", null, [chip(r.status || "–", r.status === "promoted" ? "ok" : (r.status === "prepare failed" ? "warn" : ""))]); } },
      { key: "winner", label: "Winner" },
      { key: "classification", label: "Class", cell: function (r) { return el("td", null, [r.classification ? chip(r.classification) : document.createTextNode("–")]); } },
      { key: "rho_over_winner", label: "rho / winner", numeric: true, cell: function (r) { return ratioCell(r.rho_over_winner, true); } },
      { key: "beats_rho_strict", label: "Strict", cell: function (r) { return tdText(r.beats_rho_strict === true ? "yes" : r.beats_rho_strict === false ? "no" : "–"); } },
      { key: "cells", label: "Cells", numeric: true, sort: function (r) { return (r.cells || []).length; }, cell: function (r) { return tdText(String((r.cells || []).length), "n"); } },
      { key: "candidates", label: "Arms", sort: function (r) { return (r.candidates || []).length; }, cell: function (r) { return tdText((r.candidates || []).map(function (c) { return c.id; }).join(", "), "wrap"); } },
      { key: "unit", label: "Unit", cell: function (r) { return tdText(fmt(r.unit), "mono"); } }
    ], DATA.rounds, { sortKey: "round" }));
  }
  function roundDetail(name) {
    var r = byId.round[name];
    if (!r) { view.appendChild(el("p", { class: "note", text: "No round " + name + "." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#rounds", "Rounds"), " / " + r.round]));
    view.appendChild(el("h1", { text: r.round }));
    view.appendChild(el("dl", { class: "facts" }, [
      fact("Status", r.status), fact("Winner", r.winner), fact("Classification", r.classification),
      fact("rho / winner", r.rho_over_winner === null || r.rho_over_winner === undefined ? null : fmt(r.rho_over_winner, 4)),
      fact("Beats rho (point / strict)", (r.beats_rho === undefined ? "–" : String(r.beats_rho)) + " / " + (r.beats_rho_strict === undefined || r.beats_rho_strict === null ? "–" : String(r.beats_rho_strict))),
      fact("Rho parity", r.rho_parity === undefined ? null : String(r.rho_parity)), fact("Unit", r.unit, true), fact("Scope", r.scope),
      fact("Reason (if the round did not run)", r.reason),
      fact("Files", el("span", null, [repoLink(r.dir, "directory"), r.report ? " · " : "", r.report ? repoLink(r.report, "REPORT.md") : null]))
    ]));
    if (r.candidates && r.candidates.length) {
      view.appendChild(el("h2", { text: "Arms" }));
      view.appendChild(table([
        { key: "id", label: "Arm" },
        { key: "parent", label: "Parent" },
        { key: "config", label: "Configuration", cell: function (c) { return tdText(paramsText(c.config), "mono wrap"); } },
        { key: "hypothesis", label: "Hypothesis", cell: function (c) { return tdText(fmt(c.hypothesis), "wrap"); } },
        { key: "configuration_sha256", label: "Config SHA-256", cell: function (c) { return tdText((c.configuration_sha256 || "").slice(0, 16) + (c.configuration_sha256 ? "…" : ""), "mono"); } }
      ], r.candidates));
    }
    if (r.cells && r.cells.length) {
      view.appendChild(el("h2", { text: "Cells" }));
      view.appendChild(table([
        { key: "cell", label: "Cell", cell: function (c) { return tdText(c.cell, "mono"); } },
        { key: "slug", label: "Curve", cell: function (c) { return c.slug ? tdLink("#curve/" + c.slug, c.slug) : tdText("not in the registry"); } },
        { key: "degree", label: "Degree", numeric: true },
        { key: "subgroup_order", label: "r", cell: function (c) { return tdText(fmt(c.subgroup_order), "mono n"); } },
        { key: "cofactor", label: "Cofactor", numeric: true, sort: function (c) { return Number(c.cofactor); } },
        { key: "stages", label: "Stages", cell: function (c) { return tdText(c.stages.join(", "), "wrap"); } }
      ], r.cells, { sortKey: "degree" }));
    }
  }

  /* ---------- sessions ---------- */
  function sessionsView() {
    view.appendChild(el("h1", { text: "ecbench sessions" }));
    view.appendChild(el("p", { class: "lede", text: "Every committed ecbench session. CI re-audits each with full replays; an interrupted session is kept as evidence of the interruption. Counts are per record status as the session wrote them." }));
    view.appendChild(table([
      { key: "session_id", label: "Session", cell: function (s) { return tdLink("#session/" + s.session_id, s.session_id); } },
      { key: "label", label: "Label", cell: function (s) { return tdText(fmt(s.label), "wrap"); } },
      { key: "status", label: "Status", cell: function (s) { return el("td", null, [chip(s.status || "–", s.status === "complete" ? "ok" : "warn")]); } },
      { key: "records", label: "Records", numeric: true },
      { key: "verified", label: "Verified", numeric: true, sort: function (s) { return s.status_counts.verified || 0; }, cell: function (s) { return tdText(String(s.status_counts.verified || 0), "n"); } },
      { key: "curves", label: "Curves", numeric: true, sort: function (s) { return s.curves.length; }, cell: function (s) { return tdText(String(s.curves.length), "n"); } },
      { key: "arms", label: "Arms", numeric: true, sort: function (s) { return s.arms.length; }, cell: function (s) { return tdText(String(s.arms.length), "n"); } },
      { key: "dir", label: "Directory", cell: function (s) { return el("td", { class: "mono" }, [repoLink(s.dir, s.dir.replace("research/", ""))]); } }
    ], DATA.sessions, { sortKey: "dir" }));
  }
  function sessionDetail(id) {
    var s = byId.session[id];
    if (!s) { view.appendChild(el("p", { class: "note", text: "No session " + id + "." })); return; }
    view.appendChild(el("p", { class: "crumbs" }, [link("#sessions", "Sessions"), " / " + s.session_id]));
    view.appendChild(el("h1", { text: s.session_id }));
    view.appendChild(el("p", { class: "lede", text: s.label || "" }));
    view.appendChild(el("dl", { class: "facts" }, [
      fact("Status", s.status), fact("Records", s.records), fact("By status", paramsText(s.status_counts), true), fact("Spec", s.spec_id, true),
      fact("Host class", s.env_class_id, true), fact("Binary SHA-256", s.binary_sha256, true), fact("Commit", s.git_commit, true),
      fact("Directory", repoLink(s.dir), true)
    ]));
    view.appendChild(el("h2", { text: "Arms" }));
    view.appendChild(table([
      { key: "arm", label: "Arm" }, { key: "role", label: "Role" },
      { key: "method", label: "Method", cell: function (a) { return tdLink("#method/" + a.method_id, a.method); } },
      { key: "method_id", label: "Method id", cell: function (a) { return tdText(a.method_id, "mono"); } }
    ], s.arms));
    view.appendChild(el("h2", { text: "Curves" }));
    var rows = s.curves.map(function (slug) {
      var c = byId.curve[slug];
      var ref = c ? c.ecbench.filter(function (e) { return e.session_id === id; }) : [];
      var best = ref.length ? Math.min.apply(null, ref.map(function (e) { return e.mean_s === null ? Infinity : e.mean_s; })) : null;
      return { slug: slug, r_bits: c ? c.r_bits : null, arms: ref.length, best: best === Infinity ? null : best };
    });
    view.appendChild(table([
      { key: "slug", label: "Curve", cell: function (r) { return tdLink("#curve/" + r.slug, r.slug); } },
      { key: "r_bits", label: "log₂ r", numeric: true, digits: 2 },
      { key: "arms", label: "Arms measured", numeric: true },
      { key: "best", label: "Lowest mean S", numeric: true, digits: 3 }
    ], rows, { sortKey: "r_bits" }));
  }


  /* ---------- yield ledger ---------- */
  function yieldView(q) {
    var head = el("div", null, [
      el("h1", { text: "Yield ledger" }),
      el("p", { class: "lede", text: "What each index-calculus run yielded on one target: relations per trial, lookups per relation, the relation matrix's rank, and the decomposition solver's own statistics where the harness recorded them. One row per (curve, target, factor base, oracle, solver, round), quoted from the run's record. These are stage diagnostics: they rank bases and oracles on the same targets, and only a whole-pipeline S decides speed." })
    ]);
    var f = el("form", { class: "filters" });
    function select(name, label, options) {
      var sel = el("select", { name: name });
      options.forEach(function (o) { var opt = el("option", { value: o[0], text: o[1] }); if ((q[name] || "") === o[0]) opt.selected = true; sel.appendChild(opt); });
      f.appendChild(el("label", null, [label, sel]));
      return sel;
    }
    function uniq(key) { var s = {}; DATA.yields.forEach(function (y) { if (y[key] !== null && y[key] !== undefined) s[y[key]] = 1; }); return Object.keys(s).sort(); }
    var curve = select("curve", "Curve", [["", "any"]].concat(uniq("curve").map(function (c) { return [c, c]; })));
    var fbf = select("fb_family", "Factor-base family", [["", "any"]].concat(uniq("fb_family").map(function (c) { return [c, c]; })));
    var fb = select("fb_id", "Factor base", [["", "any"]].concat(uniq("fb_id").map(function (c) { return [c, c]; })));
    var oracle = select("oracle", "Oracle", [["", "any"]].concat(uniq("oracle").map(function (c) { return [c, c]; })));
    var solver = select("solver_name", "Solver", [["", "any"]].concat(uniq("solver_name").map(function (c) { return [c, c]; })));
    var status = select("status", "Outcome", [["", "any"], ["verified", "verified"], ["exhausted", "exhausted"], ["error", "error"], ["timeout", "timeout"]]);
    var session = select("session_id", "Session", [["", "any"]].concat(uniq("session_id").map(function (c) { return [c, c]; })));
    var withSolver = el("input", { type: "checkbox", name: "has_solver" }); withSolver.checked = q.has_solver === "1";
    f.appendChild(el("label", { class: "check" }, [withSolver, "has solver statistics"]));
    var reset = el("button", { type: "button", text: "Clear filters" });
    f.appendChild(el("label", null, ["", reset]));
    var count = el("p", { class: "count" });
    var holder = el("div");
    function current() {
      var out = {};
      Array.prototype.forEach.call(f.querySelectorAll("input,select"), function (i) { if (i.type === "checkbox") { if (i.checked) out[i.name] = "1"; } else if (i.value !== "") out[i.name] = i.value; });
      return out;
    }
    function apply() {
      var s = current();
      setQuery("yield", s);
      var rows = DATA.yields.filter(function (y) {
        for (var k in s) {
          if (k === "has_solver") { if (!y.solver) return false; continue; }
          if (String(y[k]) !== s[k]) return false;
        }
        return true;
      });
      var verified = rows.filter(function (y) { return y.status === "verified" && y["yield"] !== null; });
      var meanYield = verified.length ? verified.reduce(function (a, y) { return a + y["yield"]; }, 0) / verified.length : null;
      count.textContent = rows.length + " of " + DATA.yields.length + " runs" + (meanYield === null ? "" : " · mean yield over " + verified.length + " verified: " + (100 * meanYield).toFixed(3) + " %");
      while (holder.firstChild) holder.removeChild(holder.firstChild);
      holder.appendChild(table([
        { key: "curve", label: "Curve", cell: function (y) { return tdLink("#curve/" + y.curve, y.curve); } },
        { key: "target_index", label: "Target", numeric: true, cell: function (y) { return tdText(y.target_index === null ? "–" : "#" + y.target_index, "n"); } },
        { key: "fb_id", label: "Factor base", cell: function (y) { return y.fb_id ? el("td", null, [link("#fb/" + y.fb_id, y.fb_id, "mono"), " ", chip(y.fb_family || "")]) : tdText("–"); } },
        { key: "fb_columns", label: "Cols", numeric: true },
        { key: "oracle", label: "Oracle", cell: function (y) { return tdText(fmt(y.oracle), "mono"); } },
        { key: "solver_name", label: "Solver", cell: function (y) { return tdText(fmt(y.solver_name), "mono"); } },
        { key: "status", label: "Outcome", cell: function (y) { return el("td", null, [chip(y.status, y.status === "verified" ? "ok" : "warn")]); } },
        { key: "trials", label: "Trials", numeric: true },
        { key: "relations", label: "Relations", numeric: true },
        { key: "yield", label: "Yield", numeric: true, cell: function (y) { return tdText(y["yield"] === null ? "–" : (100 * y["yield"]).toFixed(3) + " %", "n"); } },
        { key: "lookups_per_relation", label: "Lookups / rel.", numeric: true, digits: 1 },
        { key: "matrix_rank", label: "Rank / rows", numeric: true, cell: function (y) { return tdText(y.matrix_rank === null ? "–" : y.matrix_rank + " / " + y.matrix_rows, "n"); } },
        { key: "solver", label: "Solver stats", sort: function (y) { return y.solver ? 1 : 0; }, cell: function (y) {
          if (!y.solver) return tdText("–");
          var sv = y.solver, parts = [];
          ["macaulay_rows", "macaulay_columns", "macaulay_degree", "macaulay_rank", "first_fall_degree", "sat_variables", "sat_clauses", "sat_conflicts", "sat_decisions", "sat_propagations"].forEach(function (k) { if (sv[k] !== undefined && sv[k] !== null) parts.push(k.replace(/_/g, " ") + " " + sv[k]); });
          if (sv.budget_exhausted) parts.push("budget exhausted");
          return tdText(parts.join(" · ") || "(empty)", "mono wrap");
        } },
        { key: "s", label: "S", numeric: true, digits: 3 },
        { key: "session_id", label: "Session", cell: function (y) { return tdLink("#session/" + y.session_id, y.session_id); } }
      ], rows, { sortKey: "curve", empty: "No runs match." }));
    }
    f.addEventListener("input", apply);
    f.addEventListener("submit", function (e) { e.preventDefault(); apply(); });
    reset.addEventListener("click", function () { Array.prototype.forEach.call(f.querySelectorAll("input,select"), function (i) { if (i.type === "checkbox") i.checked = false; else i.value = ""; }); apply(); });
    view.appendChild(head); view.appendChild(f); view.appendChild(count); view.appendChild(holder);
    apply();
  }
  /* ---------- vocabulary ---------- */
  function vocabularyView() {
    view.appendChild(el("h1", { text: "Vocabulary" }));
    view.appendChild(el("p", { class: "lede", text: "The named parts of the pipeline. Oracles, solvers and factor-base families carry no hashed identity of their own: a name enters an identity only through an ecbench method's parameters (ECM1) or a tournament candidate code (IC1)." }));
    [["Decomposition oracles", DATA.vocabulary.oracles], ["Solvers", DATA.vocabulary.solvers], ["Factor-base families", DATA.vocabulary.factor_base_families]].forEach(function (pair) {
      view.appendChild(el("h2", { text: pair[0] }));
      view.appendChild(el("dl", { class: "vocab" }, pair[1].reduce(function (acc, v) {
        acc.push(el("dt", { text: v.name + (v.params ? " · " + v.params : "") }));
        acc.push(el("dd", { text: v.what + (v.where ? " (" + v.where + ")" : "") }));
        return acc;
      }, [])));
    });
    view.appendChild(el("h2", { text: "Identity formats" }));
    view.appendChild(el("dl", { class: "vocab" }, [
      el("dt", { text: "icv1-<f2m<m>|fp<bits>>-t<trace>-<model8>" }), el("dd", { text: "Curve slug for display; the ICV1 string plus the model JSON is the record. A negative trace is written tm." }),
      el("dt", { text: "EC1<N<deg>|P<bits>|Q<p>D<n>>C<tag>h<sha12>" }), el("dd", { text: "A representation: field, model, subgroup and generator included. Join across repositories on this, never on the slug alone." }),
      el("dt", { text: "ECM1h<sha12>, FB1h<sha12>, ECS1h, ECBS1h, ECR1h, ECC1h" }), el("dd", { text: "ecbench method, factor base, spec, session, record, comparison." }),
      el("dt", { text: "IC1N<n>C<tag>fb<points>PDP<m><solver>RC<collector>LA<la>TD<descent>ISO<k>h<sha12>" }), el("dd", { text: "Tournament candidate identity; rho candidates use urn:ec-candidate and never an IC1 label." })
    ]));
  }

  /* ---------- router ---------- */
  function route() {
    var h = parseHash();
    while (view.firstChild) view.removeChild(view.firstChild);
    Array.prototype.forEach.call(document.querySelectorAll(".views a"), function (a) {
      var v = a.getAttribute("data-view");
      var active = h.view === v || (h.view === "curve" && v === "curves") || (h.view === "method" && v === "methods") || (h.view === "fb" && v === "factor-bases") || (h.view === "candidate" && v === "candidates") || (h.view === "round" && v === "rounds") || (h.view === "session" && v === "sessions") || (h.view === "yield" && v === "yield");
      if (active) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current");
    });
    switch (h.view) {
      case "curve": curveDetail(h.id); break;
      case "methods": methodsView(); break;
      case "method": methodDetail(h.id); break;
      case "factor-bases": fbView(); break;
      case "fb": fbDetail(h.id); break;
      case "candidates": candidatesView(); break;
      case "candidate": candidateDetail(h.id); break;
      case "rounds": roundsView(); break;
      case "round": roundDetail(h.id); break;
      case "sessions": sessionsView(); break;
      case "session": sessionDetail(h.id); break;
      case "yield": yieldView(h.q); break;
      case "vocabulary": vocabularyView(); break;
      default: curvesView(h.q);
    }
    window.scrollTo(0, 0);
  }

  fetch("./data.json", { cache: "no-store" }).then(function (r) {
    if (!r.ok) throw new Error("data.json: HTTP " + r.status);
    return r.json();
  }).then(function (d) {
    DATA = d;
    d.curves.forEach(function (c) { byId.curve[c.slug] = c; });
    d.methods.forEach(function (m) { byId.method[m.method_id] = m; });
    d.factor_bases.forEach(function (f) { byId.fb[f.fb_id] = f; });
    d.candidates.forEach(function (k) { byId.candidate[k.candidate_id] = k; });
    d.rounds.forEach(function (r) { byId.round[r.round] = r; });
    d.sessions.forEach(function (s) { byId.session[s.session_id] = s; });
    loading.hidden = true;
    view.hidden = false;
    window.addEventListener("hashchange", route);
    route();
  }).catch(function (e) {
    loading.textContent = "Could not load data.json (" + e.message + "). The index is a static file beside this page; open it from the published site or a local server rather than file://.";
  });
})();
