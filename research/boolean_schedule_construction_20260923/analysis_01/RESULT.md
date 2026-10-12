# Discovery cost partition and unchanged-scan ceilings

24 cells, 10560 observations; complete and verified: True.

Ceilings are conditional timing diagnostics. No constructor optimization or speedup is measured. Null ceilings and failed comparability guards remain inconclusive.

| Variables | Family | Policy | Median optimistic ceiling | 95% interval | Comparability | Decision |
|---:|---|---|---:|---|---|---|
| 12 | cross_planted | 16 | None | None | False | INCONCLUSIVE |
| 12 | cross_planted | 64 | None | None | False | INCONCLUSIVE |
| 12 | cross_planted | dispatch | None | None | False | INCONCLUSIVE |
| 12 | planted | 16 | None | None | False | INCONCLUSIVE |
| 12 | planted | 64 | None | None | False | INCONCLUSIVE |
| 12 | planted | dispatch | None | None | False | INCONCLUSIVE |
| 12 | unplanted | 16 | None | None | False | INCONCLUSIVE |
| 12 | unplanted | 64 | None | None | False | INCONCLUSIVE |
| 12 | unplanted | dispatch | None | None | False | INCONCLUSIVE |
| 16 | cross_planted | 16 | None | None | False | INCONCLUSIVE |
| 16 | cross_planted | 64 | 1.3384903446064924 | [0.8066597831698503, 1.7845554834523036] | False | INCONCLUSIVE |
| 16 | cross_planted | dispatch | 1.8808541916382784 | [1.4786012526096033, 2.301841473178543] | False | INCONCLUSIVE |
| 16 | planted | 16 | None | None | False | INCONCLUSIVE |
| 16 | planted | 64 | None | None | False | INCONCLUSIVE |
| 16 | planted | dispatch | None | None | False | INCONCLUSIVE |
| 16 | unplanted | 16 | 1.0758095899551865 | [0.9660408163265306, 1.2712558139534884] | False | INCONCLUSIVE |
| 16 | unplanted | 64 | 1.0538051786140836 | [0.9451019066403682, 1.2210154079312958] | False | INCONCLUSIVE |
| 16 | unplanted | dispatch | 1.144628160747571 | [0.9162279343449985, 1.3145567169957413] | False | INCONCLUSIVE |
| 20 | cross_planted | 16 | 0.8831522523999955 | [0.7953755840421708, 0.918722585323865] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | cross_planted | 64 | 1.0787385646994538 | [1.0021636308745094, 1.138523681292571] | False | INCONCLUSIVE |
| 20 | cross_planted | dispatch | 0.9070953351940747 | [0.7760071690173771, 0.9433975782761338] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | planted | 16 | 0.8977611177175469 | [0.8307692307692308, 0.9662913773343221] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | planted | 64 | 1.1662202571116742 | [1.1093532952616145, 1.2648963345808828] | False | INCONCLUSIVE |
| 20 | planted | dispatch | 0.9077993350389514 | [0.8446524912726119, 0.9729135613471251] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | unplanted | 16 | 0.6092958118777874 | [0.6064014670917253, 0.6193154311259438] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | unplanted | 64 | 0.8698607594732648 | [0.8539057665674632, 0.8947033960925723] | False | INCONCLUSIVE |
| 20 | unplanted | dispatch | 0.6159416473924306 | [0.5942947368421052, 0.6366336345844057] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | cross_planted | 16 | 0.5350456585493635 | [0.5248769598724422, 0.5395872738100852] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | cross_planted | 64 | 0.795858639953201 | [0.7836518360976724, 0.8142567499704759] | False | INCONCLUSIVE |
| 24 | cross_planted | dispatch | 0.7916979584791796 | [0.7831587758359394, 0.7998720033789047] | False | INCONCLUSIVE |
| 24 | planted | 16 | 0.5474683004005213 | [0.539096919711738, 0.5515877814954664] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | planted | 64 | 0.8052138916417957 | [0.8016433849528045, 0.8158114894249948] | False | INCONCLUSIVE |
| 24 | planted | dispatch | 0.8138081624514646 | [0.8094357076780758, 0.8207587914809311] | False | INCONCLUSIVE |
| 24 | unplanted | 16 | 0.5390389794933172 | [0.5281522805091248, 0.5454329607145545] | True | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | unplanted | 64 | 0.7989559665231527 | [0.7859761318634272, 0.8079863035177093] | False | INCONCLUSIVE |
| 24 | unplanted | dispatch | 0.8099888911099468 | [0.8040893819161065, 0.8173977170061512] | False | INCONCLUSIVE |
