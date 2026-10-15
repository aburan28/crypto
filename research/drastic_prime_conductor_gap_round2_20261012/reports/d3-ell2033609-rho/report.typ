#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= Run report: d3-ell2033609-rho

#text("# Run report: d3-ell2033609-rho\n\nSession ECBS1h4919bd7b691c. Spec: ../../specs/d3-ell2033609-rho.json. Host: Apple M4 Pro, aarch64, environment ECBENV2h59034d4cd61b. Binary SHA-256: 1a87842e671358071b2950e6b780bb9d0c750c2a35684329fd5a1fc1967505fa; revision 466035c18f651bdc675806dfaa767c9a9c0073f5.\n\n| arm | curve | mean counted S | measured |\n|---|---|---:|---|\n| rho-neg-aa | icv1-fp42-t763-f4cdb4ec | 0.911225 | 32 |\n| rho-neg | icv1-fp42-t763-f4cdb4ec | 0.911225 | 32 |\n| rho-neg-aa | icv1-fp42-t763-aebf527a | 0.966335 | 32 |\n| rho-neg | icv1-fp42-t763-aebf527a | 0.966335 | 32 |\n\nMeasured failures: 0. Same-host exact replay only.\n\nUnpriced native work and wall figures stay descriptive; counted S is primary.\n")

#image("counts.svg", width: 100%)

