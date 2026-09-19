// Readable Pollard rho for GF(2^23). Same walk as ecc2k130/examples/rho_toy.py:
// R <- R + tau^j(R), j = 3 + ((HW(x)/2) mod 8). Recovers a planted discrete
// logarithm in the browser so the status site can show the method, not only
// the campaign counts.
(function (global) {
  "use strict";

  var M = 23;
  var MOD = (1 << 23) | (1 << 5) | 1;
  var DP_WEIGHT = 8;
  var J_MIN = 3;
  var N_BRANCHES = 8;

  function mul(a, b) {
    var r = 0;
    a |= 0;
    b |= 0;
    while (b) {
      if (b & 1) r ^= a;
      b >>>= 1;
      a <<= 1;
      if (a & (1 << M)) a ^= MOD;
    }
    return r >>> 0;
  }

  function sqr(a) { return mul(a, a); }

  function powField(a, e) {
    var r = 1;
    while (e) {
      if (e & 1) r = mul(r, a);
      a = sqr(a);
      e = Math.floor(e / 2);
    }
    return r;
  }

  function inv(a) { return powField(a, (1 << M) - 2); }

  function frob(a, n) {
    n = n % M;
    while (n--) a = sqr(a);
    return a;
  }

  function trace(a) {
    var t = 0, cur = a, i;
    for (i = 0; i < M; i++) {
      t ^= cur;
      cur = sqr(cur);
    }
    return t;
  }

  function halfTrace(c) {
    var z = 0, cur = c, i;
    for (i = 0; i <= (M - 1) / 2; i++) {
      z ^= cur;
      cur = sqr(sqr(cur));
    }
    return z;
  }

  function add(P, Q) {
    if (!P) return Q;
    if (!Q) return P;
    var x1 = P[0], y1 = P[1], x2 = Q[0], y2 = Q[1], lam, x3, y3;
    if (x1 === x2) {
      if ((y1 ^ y2) === x1) return null;
      lam = x1 ^ mul(y1, inv(x1));
      x3 = sqr(lam) ^ lam;
      y3 = sqr(x1) ^ mul(lam ^ 1, x3);
      return [x3, y3];
    }
    lam = mul(y1 ^ y2, inv(x1 ^ x2));
    x3 = sqr(lam) ^ lam ^ x1 ^ x2;
    y3 = mul(lam, x1 ^ x3) ^ x3 ^ y1;
    return [x3, y3];
  }

  function neg(P) { return P ? [P[0], P[0] ^ P[1]] : null; }

  function scalarMul(k, P) {
    if (k < 0) return scalarMul(-k, neg(P));
    var R = null;
    k = Math.floor(k);
    while (k) {
      if (k & 1) R = add(R, P);
      P = add(P, P);
      k = Math.floor(k / 2);
    }
    return R;
  }

  function frobPt(P, n) {
    return P ? [frob(P[0], n), frob(P[1], n)] : null;
  }

  function eqPt(A, B) {
    if (!A && !B) return true;
    if (!A || !B) return false;
    return A[0] === B[0] && A[1] === B[1];
  }

  function liftX(x) {
    if (!x) return null;
    var c = x ^ inv(sqr(x));
    if (trace(c)) return null;
    return [x, mul(x, halfTrace(c))];
  }

  function onCurve(P) {
    if (!P) return true;
    return (sqr(P[1]) ^ mul(P[0], P[1])) === (mul(sqr(P[0]), P[0]) ^ 1);
  }

  function groupOrder() {
    var t0 = 2, t1 = -1, i, n;
    for (i = 0; i < M - 1; i++) {
      n = -t1 - 2 * t0;
      t0 = t1;
      t1 = n;
    }
    return (1 << M) + 1 - t1;
  }

  function isPrime(n) {
    if (n < 2) return false;
    var d = 2;
    while (d * d <= n) {
      if (n % d === 0) return false;
      d++;
    }
    return true;
  }

  function modPow(a, e, m) {
    var r = 1;
    a %= m;
    while (e > 0) {
      if (e & 1) r = (r * a) % m;
      a = (a * a) % m;
      e = Math.floor(e / 2);
    }
    return r;
  }

  function modInv(a, m) {
    var t = 0, newt = 1, r = m, newr = ((a % m) + m) % m;
    while (newr !== 0) {
      var q = Math.floor(r / newr);
      var tmp = newt; newt = t - q * newt; t = tmp;
      tmp = newr; newr = r - q * newr; r = tmp;
    }
    if (r > 1) return null;
    return (t % m + m) % m;
  }

  function sqrtMod(a, p) {
    a = ((a % p) + p) % p;
    if (a === 0) return 0;
    if (modPow(a, (p - 1) / 2, p) !== 1) return null;
    var q = p - 1, s = 0;
    while (q % 2 === 0) { q /= 2; s++; }
    var z = 2;
    while (modPow(z, (p - 1) / 2, p) !== p - 1) z++;
    var m_ = s, c = modPow(z, q, p), t = modPow(a, q, p), r = modPow(a, (q + 1) / 2, p);
    while (t !== 1) {
      var i = 0, tt = t;
      while (tt !== 1) { tt = (tt * tt) % p; i++; }
      var b = modPow(c, 1 << (m_ - i - 1), p);
      m_ = i;
      c = (b * b) % p;
      t = (t * c) % p;
      r = (r * b) % p;
    }
    return r;
  }

  function gf2Rank(rows) {
    rows = rows.slice();
    var rank = 0, col, i, piv;
    for (col = 0; col < M; col++) {
      piv = -1;
      for (i = rank; i < rows.length; i++) {
        if ((rows[i] >>> col) & 1) { piv = i; break; }
      }
      if (piv < 0) continue;
      var tmp = rows[rank]; rows[rank] = rows[piv]; rows[piv] = tmp;
      for (i = 0; i < rows.length; i++) {
        if (i !== rank && ((rows[i] >>> col) & 1)) rows[i] ^= rows[rank];
      }
      rank++;
    }
    return rank;
  }

  function normalMasks() {
    var gamma = 0, cand, rows, g, i, k, masks, mask;
    for (cand = 2; cand < (1 << 12); cand++) {
      rows = [];
      g = cand;
      for (i = 0; i < M; i++) { rows.push(g); g = sqr(g); }
      if (gf2Rank(rows) === M) { gamma = cand; break; }
    }
    masks = [];
    g = gamma;
    for (i = 0; i < M; i++) {
      mask = 0;
      for (k = 0; k < M; k++) {
        if (trace(mul(1 << k, g))) mask |= 1 << k;
      }
      masks.push(mask);
      g = sqr(g);
    }
    return masks;
  }

  function hw(x, masks) {
    var n = 0, i;
    for (i = 0; i < masks.length; i++) {
      var bits = x & masks[i], c = 0;
      while (bits) { c ^= bits & 1; bits >>>= 1; }
      n += c;
    }
    return n;
  }

  function jOf(P, masks) {
    return ((Math.floor(hw(P[0], masks) / 2) % N_BRANCHES) + J_MIN);
  }

  function canonicalX(P) {
    var best = P[0], cur = P[0], i;
    for (i = 1; i < M; i++) {
      cur = sqr(cur);
      if (cur < best) best = cur;
    }
    return best;
  }

  function mulMod(a, b, m) { return ((a % m) * (b % m)) % m; }

  function setup(rng) {
    var n = groupOrder();
    var cofactor = 4;
    var r = n / cofactor;
    if (n % cofactor || !isPrime(r)) throw new Error("unexpected group order " + n);
    var G = null, P;
    while (!G) {
      P = liftX(rng() >>> 0 & ((1 << M) - 1));
      if (!P || !onCurve(P)) continue;
      G = scalarMul(cofactor, P);
      if (!G || !onCurve(G) || scalarMul(r, G)) G = null;
    }
    var disc = sqrtMod((1 - 8) % r, r);
    var inv2 = modInv(2, r);
    var s = null, cand, i;
    var opts = [mulMod(-1 + disc, inv2, r), mulMod(-1 - disc, inv2, r)];
    for (i = 0; i < opts.length; i++) {
      cand = (opts[i] % r + r) % r;
      if (eqPt(frobPt(G, 1), scalarMul(cand, G))) { s = cand; break; }
    }
    if (s === null) throw new Error("Frobenius eigenvalue mismatch");
    return { r: r, G: G, s: s };
  }

  function walkToDp(start, masks) {
    var P = start;
    var counts = [0, 0, 0, 0, 0, 0, 0, 0];
    var it, j;
    for (it = 0; it < 65536; it++) {
      if (hw(P[0], masks) <= DP_WEIGHT) return { P: P, steps: it, counts: counts };
      j = jOf(P, masks);
      counts[j - J_MIN]++;
      P = add(P, frobPt(P, j));
      if (!P) return null;
    }
    return null;
  }

  function multiplier(counts, s, r) {
    var mu = 1, j, factor;
    for (j = 0; j < counts.length; j++) {
      if (!counts[j]) continue;
      factor = (1 + modPow(s, j + J_MIN, r)) % r;
      mu = mulMod(mu, modPow(factor, counts[j], r), r);
    }
    return mu;
  }

  function recoverK(endA, aA, bA, endB, aB, bB, s, r, G, Q) {
    var c, rot, sc, num, den, k, invDen;
    for (c = 0; c < M; c++) {
      rot = frobPt(endB, c);
      if (eqPt(rot, endA)) sc = modPow(s, c, r);
      else if (eqPt(rot, neg(endA))) sc = (r - modPow(s, c, r)) % r;
      else continue;
      num = (aA - mulMod(sc, aB, r)) % r;
      if (num < 0) num += r;
      den = (mulMod(sc, bB, r) - bA) % r;
      if (den < 0) den += r;
      if (den === 0) continue;
      invDen = modInv(den, r);
      k = mulMod(num, invDen, r);
      if (eqPt(scalarMul(k, G), Q)) return k;
    }
    return null;
  }

  function mulberry32(seed) {
    return function () {
      seed |= 0;
      seed = seed + 0x6D2B79F5 | 0;
      var t = Math.imul(seed ^ seed >>> 15, 1 | seed);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0);
    };
  }

  function search(r, G, Q, s, masks, rng) {
    var store = {};
    var stats = { steps: 0, dps: 0, walks: 0 };
    var i, seed, start, hit, mu, a, b, key, prev, k;
    for (i = 0; i < 64; i++) {
      seed = 1 + rng() % (r - 1);
      start = add(scalarMul(seed, G), Q);
      stats.walks++;
      hit = walkToDp(start, masks);
      if (!hit) continue;
      stats.steps += hit.steps;
      stats.dps++;
      mu = multiplier(hit.counts, s, r);
      a = mulMod(mu, seed, r);
      b = mu % r;
      key = canonicalX(hit.P);
      prev = store[key];
      if (!prev) {
        store[key] = { P: hit.P, a: a, b: b, seed: seed, steps: hit.steps };
        continue;
      }
      if (prev.seed === seed) continue;
      k = recoverK(prev.P, prev.a, prev.b, hit.P, a, b, s, r, G, Q);
      if (k !== null) {
        return {
          k: k,
          stats: stats,
          seedA: prev.seed,
          seedB: seed,
          key: key,
          stepsA: prev.steps,
          stepsB: hit.steps
        };
      }
    }
    return { k: null, stats: stats };
  }

  function run(seed) {
    var rng = mulberry32(seed || 1);
    var curve = setup(rng);
    var planted = 1 + rng() % (curve.r - 1);
    var Q = scalarMul(planted, curve.G);
    var masks = normalMasks();
    var found = search(curve.r, curve.G, Q, curve.s, masks, rng);
    return {
      r: curve.r,
      planted: planted,
      recovered: found.k,
      verified: found.k === planted,
      stats: found.stats,
      seedA: found.seedA,
      seedB: found.seedB,
      key: found.key,
      stepsA: found.stepsA,
      stepsB: found.stepsB
    };
  }

  global.RhoToy = { run: run, M: M, DP_WEIGHT: DP_WEIGHT };
})(typeof window !== "undefined" ? window : this);
