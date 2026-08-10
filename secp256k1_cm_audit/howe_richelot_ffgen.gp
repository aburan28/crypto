\\ howe_richelot_ffgen.gp
\\
\\ Clean-room reimplementation of the Z/3Z-symmetric Richelot-dual
\\ construction from howe_5pairs_v2.gp, using PARI's NATIVE ffgen()
\\ finite-field type for F_{p^3} arithmetic instead of hand-rolled
\\ f3add/f3mul/f3inv routines.
\\
\\ Motivation: howe_5pairs_v2.gp itself was a bugfix for howe_5pairs.gp,
\\ which silently used the global 'p' (secp256k1 prime) inside its
\\ F_{p^3} helpers instead of a local toy prime (2026-07-27 log). Hand
\\ -rolled extension-field arithmetic is a recurring bug source; ffgen()
\\ lets PARI do it natively in ~10 lines instead of ~90.
\\
\\ Validated against the known-good p=43 case from howe_richelot_v5.gp /
\\ howe_5pairs_v2.gp Test 1: expects (aa,bb) = (41,5).
\\
\\ NOTE ON gp SCRIPT SYNTAX: this gp build (2.15.4) does not accept a
\\ for()/if() body split across multiple lines when read from a file --
\\ it parses line-by-line and throws "unexpected end of file". Keep
\\ control-flow bodies on ONE line (documented here since it cost time
\\ to discover; see RESEARCH_AUTOLAB_LOG.md 2026-08-10 entry).
\\
\\ Run: gp -q howe_richelot_ffgen.gp

\\ Classical Richelot dual: given a quadratic G1 = X^2 - sv*X + qv over
\\ F_{p^3}, and its Z/3-Galois-conjugate triple (G1,G2,G3) obtained by
\\ scaling by a primitive cube root of unity z3, compute the dual
\\ sextic y^2 = X^6 + aa*X^3 + bb via the standard bracket formula
\\ H_i = G_j * G_k' - G_j' * G_k (cyclic in i,j,k).
richelot_dual(sv, qv, z3) = {
  my(G1,G2,G3,H1,H2,H3,P,lc);
  G1 = 'X^2 - sv*'X + qv;
  G2 = 'X^2 - z3*sv*'X + z3^2*qv;
  G3 = 'X^2 - z3^2*sv*'X + z3*qv;
  H1 = G2*G3' - G2'*G3;
  H2 = G3*G1' - G3'*G1;
  H3 = G1*G2' - G1'*G2;
  P = H1*H2*H3;
  lc = polcoeff(P,6);
  [lc, polcoeff(P,3)/lc, polcoeff(P,0)/lc]
};

\\ Given rr = alpha^3 (the "cube" defining the source Z/3-symmetric
\\ curve y^2=(x^3-r1)(x^3-r2)), find a cube root alpha in F_{p^3}.
cube_root_in_fp3(rr, pp, t) = { my(c); c = rr + 0*t; polrootsmod('X^3 - c)[1] };

print("================================================================");
print("Test 1: p=43, expect (aa,bb) = (41,5) [howe_5pairs_v2.gp Test 1]");
print("================================================================");
{
  pp = 43;
  t = ffgen(pp^3, 't);
  alpha = cube_root_in_fp3(36, pp, t);            \\ alpha^3 = 36 = -7 mod 43
  z3p = 0; for(k=1,pp-1, if(Mod(k,pp)^3==1 && k!=1, z3p=k; break));
  z3 = z3p + 0*t;
  sv = 3*alpha; qv = 2*alpha^2;                    \\ matches sv=[0,3,0], qv=[0,0,2]
  res = richelot_dual(sv, qv, z3);
  print("  alpha^3 = ", alpha^3, "  (expect 36)");
  print("  (aa,bb) = (", res[2], ",", res[3], ")  (expect (41,5))");
  print("  MATCH: ", res[2]==41 && res[3]==5);
}
quit;
