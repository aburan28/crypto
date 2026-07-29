\\ glue_realization_verify.gp
\\ Thread 22: independent PARI cross-check of glue_realization_search.py.
\\
\\ The Python search computes the Frobenius charpoly of each candidate genus-2
\\ curve from #C(F_p) and #C(F_{p^2}) via Newton's identities.  This script
\\ re-derives the same charpoly with PARI's hyperellcharpoly() and checks it
\\ equals (T^2 - t_i T + p)(T^2 - t_j T + p) for the claimed pair (i,j).
\\
\\ Input: secp256k1_cm_audit/glue_realization_models.txt
\\        lines "p i j t_i t_j e1 e2 howe lead coeffs(hi->lo, comma separated)"
\\
\\ Run: gp -q secp256k1_cm_audit/glue_realization_verify.gp

default(parisize, 256000000);

fname = "secp256k1_cm_audit/glue_realization_models.txt";
lines = readstr(fname);

print("================================================================");
print("Thread 22: PARI cross-check of glue-realization models");
print("================================================================");

nok = 0; nbad = 0; nhowe = 0;
{
for(li = 1, #lines,
  s = lines[li];
  if(#s == 0, next);
  fields = strsplit(s, " ");
  if(fields[1] == "#", next);
  p    = eval(fields[1]);
  ii   = eval(fields[2]);
  jj   = eval(fields[3]);
  ti   = eval(fields[4]);
  tj   = eval(fields[5]);
  e1   = eval(fields[6]);
  e2   = eval(fields[7]);
  howe = eval(fields[8]);
  lead = eval(fields[9]);
  cf   = eval(concat(concat("[", fields[10]), "]"));   \\ hi -> lo

  \\ rebuild f(x) = lead * sum cf[k] x^(d-k+1)
  d = #cf - 1;
  f = 0;
  for(k = 1, #cf, f = f + cf[k] * x^(d - k + 1));
  f = lead * f;
  f = Mod(1,p) * f;

  cp  = hyperellcharpoly(f);                       \\ PARI's Frobenius charpoly
  tgt = (x^2 - ti*x + p) * (x^2 - tj*x + p);
  \\ also check the Newton-identity coefficients the Python code reported
  newton_ok = (polcoeff(tgt,3) == -e1) && (polcoeff(tgt,2) == e2);

  if(cp == tgt && newton_ok,
     nok++; if(howe, nhowe++),
     nbad++;
     print("  MISMATCH p=",p," pair(",ii,",",jj,")");
     print("    f       = ", lift(f));
     print("    charpoly= ", cp);
     print("    target  = ", tgt));
)
}

print("");
print("models checked : ", nok + nbad);
print("  agree with hyperellcharpoly : ", nok);
print("  mismatched                  : ", nbad);
print("  (of the agreeing, Howe-qualifying pairs: ", nhowe, ")");
print("");
{
if(nbad == 0,
   print("RESULT: all realizations confirmed independently by PARI. OK"),
   print("RESULT: FAILURES PRESENT -- see above."));
}
