#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= Run report: d3-ell1016839-rho

#text("# Run report: d3-ell1016839-rho\n\nSession ECBS1ha911e4f17fd2. Spec: ../../specs/d3-ell1016839-rho.json. Host: Apple M4 Pro, aarch64, environment ECBENV2h59034d4cd61b. Binary SHA-256: 1a87842e671358071b2950e6b780bb9d0c750c2a35684329fd5a1fc1967505fa; revision 466035c18f651bdc675806dfaa767c9a9c0073f5.\n\n| arm | curve | mean counted S | measured |\n|---|---|---:|---|\n| rho-neg-aa | icv1-fp40-tm1043-e27d97c8 | 0.897986 | 32 |\n| rho-neg | icv1-fp40-tm1043-e27d97c8 | 0.897986 | 32 |\n| rho-neg-aa | icv1-fp40-tm1043-a68449f5 | 0.850354 | 32 |\n| rho-neg | icv1-fp40-tm1043-a68449f5 | 0.850354 | 32 |\n\nMeasured failures: 0. Same-host exact replay only.\n\nUnpriced native work and wall figures stay descriptive; counted S is primary.\n")

#image("counts.svg", width: 100%)

