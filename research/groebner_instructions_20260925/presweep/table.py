import json,glob,sys,pathlib
S=pathlib.Path(sys.argv[1])
SUITES=["frozen","chain","chain-holdout","chain-holdout-2","r2-holdout"]
P={"system":"crypto_lib::cryptanalysis::koblitz_groebner::build_decomposition_system",
   "roots":"crypto_lib::cryptanalysis::inherited_f4::ReducedBasis::from_system_with",
   "specialise":"crypto_lib::cryptanalysis::koblitz_groebner::InheritedBases::specialise_owned",
   "insert":"crypto_lib::cryptanalysis::inherited_f4::ReducedBasis::insert",
   "lin_elim":"crypto_lib::cryptanalysis::koblitz_groebner::eliminate_linear_generators",
   "echelon":["crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_suffix_counted","crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_m4ri_counted","crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_counted"]}
def g(d,k):
    v=P[k]
    if isinstance(v,list): return sum(d["inclusive_ir"].get(x,0) for x in v)
    return d["inclusive_ir"].get(v,0)
print("| suite | rung | word ops | Ir in groebner_decompose | Ir / word op | system % | roots % | specialise % (insert %) | lin. elim % | echelon kernels % | ref Ir → cand Ir | Ir ratio | word-op ratio |")
print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|")
T={"support":[0,0],"occurring":[0,0]}
for suite in SUITES:
    rungs=sorted((S/"out"/suite/"support").glob("rung*/summary.json"), key=lambda p:int(p.parent.name[4:]))
    st={"support":[0,0],"occurring":[0,0]}
    for p in rungs:
        d=json.load(open(p))
        if d.get("skipped"): continue
        r=json.load(open(S/"out"/suite/"occurring"/p.parent.name/"summary.json"))
        t=d["ir_total"]
        pc=lambda k: 100*g(d,k)/t
        print(f"| {suite} | {d['rung']} | {d['word_ops']:,} | {t:,} | {t/d['word_ops']:.0f} | {pc('system'):.0f} | {pc('roots'):.0f} | {pc('specialise'):.0f} ({pc('insert'):.0f}) | {pc('lin_elim'):.0f} | {pc('echelon'):.1f} | {r['ir_total']:,} → {t:,} | {r['ir_total']/t:.2f}× | {r['word_ops']/d['word_ops']:.2f}× |")
        for a,x in (("support",d),("occurring",r)):
            st[a][0]+=x["word_ops"]; st[a][1]+=x["ir_total"]; T[a][0]+=x["word_ops"]; T[a][1]+=x["ir_total"]
    print(f"| **{suite} total** | | {st['support'][0]:,} | {st['support'][1]:,} | {st['support'][1]/st['support'][0]:.0f} | | | | | | {st['occurring'][1]:,} → {st['support'][1]:,} | **{st['occurring'][1]/st['support'][1]:.2f}×** | {st['occurring'][0]/st['support'][0]:.2f}× |")
