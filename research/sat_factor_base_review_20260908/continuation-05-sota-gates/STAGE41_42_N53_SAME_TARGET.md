# Stages 41–42: n=53 boundary and same-target successor

Stage 41 executed the boundary autolab's first n=53 Koblitz-vs-rho probe. Run `20260912T151448Z-f8c50b0c35` constructed 19,928 factor-base points and 188 orbit columns, reached rank 189 plus 32 surplus relations, recovered the public synthetic scalar, and passed the schema-v2 claim check. The direct pointwise arm used 873.382 seconds wall versus 5.267 seconds for packed signed-Frobenius rho. Its status remains `PENDING_INDEPENDENT_VALIDATION`; the old autolab recorded the 16 GiB cap without a watchdog or peak-RSS measurement.

The measured tuning ladder kept the direct seed fixed. `fiber_batch_64` reduced the signed-quotient arm to 483.867 seconds and measured 669 MB peak RSS. Moving the signed Frobenius expansion into the support table removed 7.24 billion query canonicalization maps: signed-expanded `fiber_batch_64` took 63.872 seconds at 6.195 GB peak RSS. Four-summand `pair_pair_256` reduced the same eta-1/16 run to 35.606 seconds.

The m=4 base-size sweep measured eta denominators 96, 128, 192, and 1024. Eta 1/128 was the minimum: 94 orbit columns, 9,964 points, a 1.409 GB support table, 5.569 seconds setup, 15.204 seconds collection, and 21.520 seconds whole-process wall. Eta 1/96 took 25.176 seconds, 1/192 took 28.241 seconds, and 1/1024 took 66.078 seconds. Increasing the pair-pair batch from 256 to 1024 was also slower, 22.265 seconds, so neither change was retained.

The local same-target validation used a fresh public scalar and independent algorithm seeds. Direct and rho constructed the same point `(4831897703993312, 2257074655544889)` and recovered scalar `3499506542931`. Eta-1/128 signed-expanded `pair_pair_256` took 25.373 seconds wall, 25.358 core-seconds, and 1.561 GB peak RSS. Packed signed-Frobenius rho took 8.534 seconds wall, 8.517 core-seconds, and 174 MB peak RSS. Direct therefore remained 2.973 times slower.

Stage 42 adds an explicit public validation scalar to both producers. Each producer still consumes its ordinary seeded scalar draw before applying the explicit value, so the subsequent random stream is unchanged. The collector receives only the point `Q`; the scalar is retained by the validator and emitted with `fixture_scalar_source = explicit_public_validation_scalar`. The production runner uses untouched SHA-256-derived scalar and algorithm seeds, meters the clean build plus both arms with watchdogs and process-tree RSS, and refuses any target-coordinate or recovered-scalar mismatch.

These stages establish a finite n=53 construction, rank, tuning, and same-target comparison. They do not use imported points, establish an unknown-scalar n=53 attack, change the asymptotic exponent, satisfy licensed Magma or unaffiliated reproduction, or establish a Koblitz index-calculus SOTA result.

## Hosted same-target result

[Run `34705118094`](https://github.com/aburan28/crypto/actions/runs/34705118094) completed from merge commit `0654c7dd16b8dfd7a61b8770583ac39316b5ff6b` and passed an independent replay after download. Both arms constructed the public point `(2565091273463387, 5885236316843894)` and recovered validation scalar `476811900269`.

The direct arm used 9,964 factor-base points, 94 orbit columns, and a 95-column relation matrix. It reached full rank at relation 157 and retained 32 surplus relations, so 189 of 189 trials produced verified four-summand relations. Support-table setup took 7.399 seconds and collection took 31.132 seconds. Whole-process direct wall was 39.781 seconds, 39.776 core-seconds, and 1,526,292,480 bytes peak RSS.

Packed signed-Frobenius rho used the same target and took 3.279 seconds wall, 3.278 core-seconds, and 35,438,592 bytes peak RSS. Direct was therefore 12.131 times slower by whole-process wall and 11.993 times slower using the producer-reported algorithm charges. The clean build plus direct cost was 37.698 times rho wall.

The fresh build plus both scientific arms used 126.968 seconds sequential wall, 272.392 core-seconds, and at most 1,718,968,320 bytes sampled process-tree RSS. The tuning archive retains the original 873.382-second pointwise control, the signed-quotient batch candidates, the signed-expanded and eta ladders, the rejected 1024-entry query width, and the separate local same-target validation. The hosted artifact remains the admitted n=53 result.
