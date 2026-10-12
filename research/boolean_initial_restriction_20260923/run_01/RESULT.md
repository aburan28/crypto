# Complete standalone initial-restriction comparison

144 fixed cells, 28224 observations; all completed and verified: True.

The current reference includes all 22 predecessor methods, including its three SIMD treatments. Every cold solve and validation cost is charged. Ratios below use matched repetitions; censored costs remain null.

| Candidate | Split | Variables | Family | Reference | Median ratio | 95% interval | Pass |
|---|---|---:|---|---|---:|---|---|
| gray_delta_simd | regression | 16 | planted | retained_frontier | 0.8924717807473765 | [0.835838565482539, 0.9816527434600149] | False |
| gray_delta_simd | regression | 16 | planted | same_policy | 0.8924717807473765 | [0.835838565482539, 0.9816527434600149] | False |
| gray_delta_simd | regression | 16 | cross_planted | retained_frontier | 0.905234965034965 | [0.8510801760198744, 0.9652041161586647] | False |
| gray_delta_simd | regression | 16 | cross_planted | same_policy | 0.905234965034965 | [0.8510801760198744, 0.9652041161586647] | False |
| gray_delta_simd | regression | 16 | unplanted | retained_frontier | 1.009129129129129 | [0.9791299238623774, 1.0479681978798587] | False |
| gray_delta_simd | regression | 16 | unplanted | same_policy | 1.015502249913465 | [0.9864930404447437, 1.1094792717121913] | False |
| gray_delta_simd | regression | 20 | planted | retained_frontier | 1.040446542118978 | [0.9920344001389436, 1.071541096035111] | False |
| gray_delta_simd | regression | 20 | planted | same_policy | 1.0658487514206925 | [1.0286109793820521, 1.0975364716585787] | False |
| gray_delta_simd | regression | 20 | cross_planted | retained_frontier | 0.9837323028994819 | [0.9174087178684356, 1.0446941228020874] | False |
| gray_delta_simd | regression | 20 | cross_planted | same_policy | 1.0485950455210309 | [1.0406152846392993, 1.0614334669899292] | False |
| gray_delta_simd | regression | 20 | unplanted | retained_frontier | 1.090879237900519 | [1.0585396567447205, 1.237882170416441] | False |
| gray_delta_simd | regression | 20 | unplanted | same_policy | 1.096190175032726 | [1.0769840216288138, 1.2517375001147122] | True |
| gray_delta_simd | regression | 24 | planted | retained_frontier | 1.0708094850234913 | [1.0469407446349979, 1.0877219481802136] | False |
| gray_delta_simd | regression | 24 | planted | same_policy | 1.1076232530767358 | [1.0837480973915123, 1.1518925142676204] | True |
| gray_delta_simd | regression | 24 | cross_planted | retained_frontier | 0.8192334713209724 | [0.6390397264693413, 1.0607151747078114] | False |
| gray_delta_simd | regression | 24 | cross_planted | same_policy | 1.1044835290277357 | [1.0967304745911854, 1.1423854136744125] | True |
| gray_delta_simd | regression | 24 | unplanted | retained_frontier | 1.1012999577561255 | [1.0781456517957202, 1.118776255128259] | False |
| gray_delta_simd | regression | 24 | unplanted | same_policy | 1.1036406442503153 | [1.0927488063104385, 1.1580758360314416] | True |
| gray_delta_simd | holdout | 16 | planted | retained_frontier | 0.8256865815085519 | [0.7844446111152779, 0.9210686337354547] | False |
| gray_delta_simd | holdout | 16 | planted | same_policy | 0.9837500100642564 | [0.9295250971776238, 1.1437785944648615] | False |
| gray_delta_simd | holdout | 16 | cross_planted | retained_frontier | 0.9704137022397892 | [0.9051743119266055, 1.0032367569274492] | False |
| gray_delta_simd | holdout | 16 | cross_planted | same_policy | 0.987183595855029 | [0.9051743119266055, 1.0987153482082488] | False |
| gray_delta_simd | holdout | 16 | unplanted | retained_frontier | 1.0746492296732786 | [0.9801243441481957, 1.171302752293578] | False |
| gray_delta_simd | holdout | 16 | unplanted | same_policy | 1.0746492296732786 | [0.9801243441481957, 1.171302752293578] | False |
| gray_delta_simd | holdout | 20 | planted | retained_frontier | 0.9007621258392322 | [0.7661127547211547, 1.0933659798577104] | False |
| gray_delta_simd | holdout | 20 | planted | same_policy | 1.0987592226032474 | [1.043435479010129, 1.2586656727934051] | False |
| gray_delta_simd | holdout | 20 | cross_planted | retained_frontier | 0.9513260952650577 | [0.7387801418439717, 1.0535066330157172] | False |
| gray_delta_simd | holdout | 20 | cross_planted | same_policy | 1.043508817116631 | [1.0047375886524823, 1.081330921827003] | False |
| gray_delta_simd | holdout | 20 | unplanted | retained_frontier | 1.0738944898860847 | [1.0170198656587979, 1.0992369380315918] | False |
| gray_delta_simd | holdout | 20 | unplanted | same_policy | 1.1278906999864082 | [1.083336898395722, 1.2937066054452853] | True |
| gray_delta_simd | holdout | 24 | planted | retained_frontier | 0.7816943603114942 | [0.44109539703851086, 1.1035661386624918] | False |
| gray_delta_simd | holdout | 24 | planted | same_policy | 1.079668735234292 | [1.0560605416764852, 1.2206838480271025] | True |
| gray_delta_simd | holdout | 24 | cross_planted | retained_frontier | 1.0254697555106573 | [0.9301996127556094, 1.0968384279475982] | False |
| gray_delta_simd | holdout | 24 | cross_planted | same_policy | 1.096922199049251 | [1.0717815278221394, 1.3007302765473483] | True |
| gray_delta_simd | holdout | 24 | unplanted | retained_frontier | 1.107621238380999 | [1.0730458123330964, 1.310659335027026] | False |
| gray_delta_simd | holdout | 24 | unplanted | same_policy | 1.107621238380999 | [1.0730458123330964, 1.3120134793020166] | True |
| initial_simd | regression | 16 | planted | retained_frontier | 0.9893938485572751 | [0.9466688184589443, 1.0546509943924658] | False |
| initial_simd | regression | 16 | planted | same_policy | 4.389318391064825 | [3.988273793291519, 5.044168557328719] | True |
| initial_simd | regression | 16 | cross_planted | retained_frontier | 0.5318291700241741 | [0.481417147726648, 0.5862189140496142] | False |
| initial_simd | regression | 16 | cross_planted | same_policy | 4.439668039212481 | [4.2840395167682885, 4.5061507183053315] | True |
| initial_simd | regression | 16 | unplanted | retained_frontier | 1.0335761603422045 | [1.0186915887850467, 1.0754702514314165] | False |
| initial_simd | regression | 16 | unplanted | same_policy | 6.228122171945701 | [6.07839027457556, 6.300970873786408] | True |
| initial_simd | regression | 20 | planted | retained_frontier | 1.0589501097954408 | [1.0434907676551477, 1.0851532424340693] | False |
| initial_simd | regression | 20 | planted | same_policy | 7.050219992117116 | [6.909216334917913, 7.115204427854506] | True |
| initial_simd | regression | 20 | cross_planted | retained_frontier | 1.6486105539344647 | [1.2564288669497685, 1.8548745321985636] | False |
| initial_simd | regression | 20 | cross_planted | same_policy | 5.774326528351608 | [5.36575705326665, 6.199895424836601] | True |
| initial_simd | regression | 20 | unplanted | retained_frontier | 1.079255774896125 | [1.0592243304830449, 1.24583800400324] | False |
| initial_simd | regression | 20 | unplanted | same_policy | 7.17293472195297 | [7.124378811662311, 7.2457468074758635] | True |
| initial_simd | regression | 24 | planted | retained_frontier | 1.079834607662264 | [1.049931295008236, 1.0950833134107334] | False |
| initial_simd | regression | 24 | planted | same_policy | 7.250469468994344 | [7.214353553665494, 7.288997337375026] | True |
| initial_simd | regression | 24 | cross_planted | retained_frontier | 1.8569940208198603 | [1.3123840688041577, 2.190404870624049] | False |
| initial_simd | regression | 24 | cross_planted | same_policy | 7.124147641306402 | [7.082125455944211, 7.14384370617773] | True |
| initial_simd | regression | 24 | unplanted | retained_frontier | 1.0908927611879327 | [1.0742044942473865, 1.110761456193] | False |
| initial_simd | regression | 24 | unplanted | same_policy | 7.214868480255771 | [7.164867055265263, 7.260429077801611] | True |
| initial_simd | holdout | 16 | planted | retained_frontier | 0.9166915668311018 | [0.815694626474443, 1.0544479262512048] | False |
| initial_simd | holdout | 16 | planted | same_policy | 5.465949016100179 | [4.930232558139535, 6.152860958435985] | True |
| initial_simd | holdout | 16 | cross_planted | retained_frontier | 0.795049687796483 | [0.6009470512268618, 0.9768712871287129] | False |
| initial_simd | holdout | 16 | cross_planted | same_policy | 4.627129591992986 | [4.1074042186827375, 5.260752475247525] | True |
| initial_simd | holdout | 16 | unplanted | retained_frontier | 1.135611650485437 | [1.022247619047619, 1.222247619047619] | False |
| initial_simd | holdout | 16 | unplanted | same_policy | 6.273028571428571 | [6.198277664982836, 6.394448241816663] | True |
| initial_simd | holdout | 20 | planted | retained_frontier | 0.887836388902425 | [0.77842864048682, 1.094953488372093] | False |
| initial_simd | holdout | 20 | planted | same_policy | 7.02970535504176 | [6.831688799158302, 7.10384561664882] | True |
| initial_simd | holdout | 20 | cross_planted | retained_frontier | 1.9107446729501896 | [0.6107764705882353, 3.5458000760235056] | False |
| initial_simd | holdout | 20 | cross_planted | same_policy | 5.582181435718617 | [4.5475557524737855, 6.330141164786834] | True |
| initial_simd | holdout | 20 | unplanted | retained_frontier | 1.0534478924143045 | [0.9926304118956891, 1.0792507511667058] | False |
| initial_simd | holdout | 20 | unplanted | same_policy | 7.041061672283874 | [6.98036037216646, 7.171138829055976] | True |
| initial_simd | holdout | 24 | planted | retained_frontier | 0.7829786284996589 | [0.441055834767642, 1.1120998591477131] | False |
| initial_simd | holdout | 24 | planted | same_policy | 7.245116529137329 | [7.1831893995856015, 7.352432477923733] | True |
| initial_simd | holdout | 24 | cross_planted | retained_frontier | 1.9915815758584783 | [1.8354882590450325, 2.15225252157885] | False |
| initial_simd | holdout | 24 | cross_planted | same_policy | 7.122178718447813 | [6.99370841976655, 7.320346578572343] | True |
| initial_simd | holdout | 24 | unplanted | retained_frontier | 1.1050226087801596 | [1.073767249350535, 1.3097393786614835] | False |
| initial_simd | holdout | 24 | unplanted | same_policy | 7.177772880000634 | [7.162665614078461, 7.247316183845821] | True |
| packed_gray16_delta_simd | regression | 16 | planted | retained_frontier | 0.675045217759318 | [0.6431889027431421, 0.722831866949553] | False |
| packed_gray16_delta_simd | regression | 16 | planted | same_policy | 1.0242248677248678 | [0.9935318188314212, 1.055490909090909] | False |
| packed_gray16_delta_simd | regression | 16 | cross_planted | retained_frontier | 0.6797785869285373 | [0.6022275258552108, 0.7376073870790001] | False |
| packed_gray16_delta_simd | regression | 16 | cross_planted | same_policy | 0.9942788812968024 | [0.9629022452076627, 1.0379590587816914] | False |
| packed_gray16_delta_simd | regression | 16 | unplanted | retained_frontier | 0.8905887483853423 | [0.87929521745021, 0.9108838079957213] | False |
| packed_gray16_delta_simd | regression | 16 | unplanted | same_policy | 1.1474977697591657 | [1.0847157597209713, 1.1756548098871713] | True |
| packed_gray16_delta_simd | regression | 20 | planted | retained_frontier | 0.7053066471630658 | [0.35342143681100513, 0.8503600518689667] | False |
| packed_gray16_delta_simd | regression | 20 | planted | same_policy | 1.0412869405080065 | [1.0273835880301838, 1.0641541031176152] | False |
| packed_gray16_delta_simd | regression | 20 | cross_planted | retained_frontier | 0.5642469090737581 | [0.31948456417181964, 0.7922090458011455] | False |
| packed_gray16_delta_simd | regression | 20 | cross_planted | same_policy | 1.0328936062734728 | [1.0156636224419957, 1.047366741219883] | False |
| packed_gray16_delta_simd | regression | 20 | unplanted | retained_frontier | 0.9384969576874527 | [0.8434577006507592, 0.9956234626825375] | False |
| packed_gray16_delta_simd | regression | 20 | unplanted | same_policy | 1.024702199253286 | [1.0096206395995342, 1.0512729730244899] | False |
| packed_gray16_delta_simd | regression | 24 | planted | retained_frontier | 0.738834514836687 | [0.6443470738006282, 0.8839815970631743] | False |
| packed_gray16_delta_simd | regression | 24 | planted | same_policy | 1.0394290216164248 | [1.0160362294378165, 1.0617548819102347] | False |
| packed_gray16_delta_simd | regression | 24 | cross_planted | retained_frontier | 1.0018106261832544 | [0.932753284975967, 1.0299652175079639] | False |
| packed_gray16_delta_simd | regression | 24 | cross_planted | same_policy | 1.0313568344019393 | [1.0218470782547544, 1.0461721138940114] | False |
| packed_gray16_delta_simd | regression | 24 | unplanted | retained_frontier | 0.8575040688998645 | [0.8402004677942267, 0.9328249978898675] | False |
| packed_gray16_delta_simd | regression | 24 | unplanted | same_policy | 1.024203573079338 | [1.017701718204073, 1.0384090653751183] | False |
| packed_gray16_delta_simd | holdout | 16 | planted | retained_frontier | 0.731699126654829 | [0.6922749218971883, 0.7423070958489872] | False |
| packed_gray16_delta_simd | holdout | 16 | planted | same_policy | 1.0716040525263066 | [1.0420613868889732, 1.207425486237529] | False |
| packed_gray16_delta_simd | holdout | 16 | cross_planted | retained_frontier | 0.8039897700549035 | [0.7128813992084769, 0.8496913366164944] | False |
| packed_gray16_delta_simd | holdout | 16 | cross_planted | same_policy | 1.0363511514978807 | [0.9634330233727068, 1.1353842073857985] | False |
| packed_gray16_delta_simd | holdout | 16 | unplanted | retained_frontier | 1.0095769945356854 | [0.862634631983376, 1.046048780487805] | False |
| packed_gray16_delta_simd | holdout | 16 | unplanted | same_policy | 1.1723005086735359 | [1.0615548186011718, 1.2108619896114143] | True |
| packed_gray16_delta_simd | holdout | 20 | planted | retained_frontier | 0.578234023839398 | [0.26144572176620884, 1.005128205128205] | False |
| packed_gray16_delta_simd | holdout | 20 | planted | same_policy | 1.0396218532283037 | [1.0061100961007772, 1.1107173694705612] | False |
| packed_gray16_delta_simd | holdout | 20 | cross_planted | retained_frontier | 0.5799776027361773 | [0.240602844401557, 1.0775851367567344] | False |
| packed_gray16_delta_simd | holdout | 20 | cross_planted | same_policy | 1.0445714087848228 | [0.998196264701671, 1.180851815061883] | False |
| packed_gray16_delta_simd | holdout | 20 | unplanted | retained_frontier | 0.9942901858602803 | [0.8468664850136239, 1.0214869730622078] | False |
| packed_gray16_delta_simd | holdout | 20 | unplanted | same_policy | 1.0442663851829117 | [1.0041415116167647, 1.1143031051214787] | False |
| packed_gray16_delta_simd | holdout | 24 | planted | retained_frontier | 0.773030003372592 | [0.4820931443976617, 1.023084924224521] | False |
| packed_gray16_delta_simd | holdout | 24 | planted | same_policy | 1.0218499078905239 | [1.0124441595898093, 1.0383787475088402] | False |
| packed_gray16_delta_simd | holdout | 24 | cross_planted | retained_frontier | 0.8273343404654938 | [0.5900241203750576, 1.0181575330766588] | False |
| packed_gray16_delta_simd | holdout | 24 | cross_planted | same_policy | 1.0209722206586913 | [1.005771689497717, 1.0447157492541543] | False |
| packed_gray16_delta_simd | holdout | 24 | unplanted | retained_frontier | 0.8463944358086548 | [0.8173157263895335, 0.9890716468486818] | False |
| packed_gray16_delta_simd | holdout | 24 | unplanted | same_policy | 1.0207831894959936 | [1.0101934657270981, 1.0369362435037002] | False |

All gates: `{"gray_delta_simd": {"incremental": "REJECTED", "retained_frontier": "REJECTED", "same_policy": "REJECTED"}, "initial_simd": {"incremental": "REJECTED", "retained_frontier": "REJECTED", "same_policy": "PASS"}, "packed_gray16_delta_simd": {"incremental": "REJECTED", "retained_frontier": "REJECTED", "same_policy": "REJECTED"}}`

These are bounded generated Boolean solves. No calibrated-operation, production index-calculus, asymptotic or rho claim is established.
