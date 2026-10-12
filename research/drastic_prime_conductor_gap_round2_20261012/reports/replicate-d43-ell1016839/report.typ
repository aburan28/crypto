#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= Run report: replicate-d43-ell1016839

#text("# Run report: replicate-d43-ell1016839\n\nSession ECBS1hc6ff884cac1d. Spec: ../../specs/replicate-d43-ell1016839.json. Host: Apple M4 Pro, aarch64, environment ECBENV2h59034d4cd61b. Binary SHA-256: 1a87842e671358071b2950e6b780bb9d0c750c2a35684329fd5a1fc1967505fa; revision 466035c18f651bdc675806dfaa767c9a9c0073f5.\n\n| arm | curve | mean counted S | measured |\n|---|---|---:|---|\n| rho-neg-aa | icv1-fp44-tm159-8117e199 | 1.06824 | 32 |\n| rho-neg | icv1-fp44-tm159-8117e199 | 1.06824 | 32 |\n| rho-neg-aa | icv1-fp44-tm159-0b4f7698 | 1.01143 | 32 |\n| rho-neg | icv1-fp44-tm159-0b4f7698 | 1.01143 | 32 |\n\nMeasured failures: 0. Same-host exact replay only.\n\nUnpriced native work and wall figures stay descriptive; counted S is primary.\n")

#image("counts.svg", width: 100%)

