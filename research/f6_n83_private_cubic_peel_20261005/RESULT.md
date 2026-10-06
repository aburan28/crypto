# n83 F6 private cubic certificate after quartic peeling

The [preregistered sampled gate](PROTOCOL.md) completed on the exact
K0 five-summand system at ordinary public T001, torsion offset zero.
It first reproduced #1446's full `k=90` private-quartic certificate:
29,880 prolonged rows, 14,058 certified, 15,822 unresolved and the
frozen witness digest
`153342de74d4e40a2ea5f96ba35d24b7bdf2726700ed9b04714c9c89ca49dbdd`.

Among those 15,822 survivors, all contained degree-three terms. They
had 11,578,544 exact degree-three term occurrences. Sampling up to 32
evenly spaced candidates per row yielded 266,144 distinct candidate
columns. Counting those candidates against **every** surviving product
and every original equation certified only **719** more rows as having
a private cubic monomial. The 15,103 remaining rows exceed the frozen
2,000-row reduction threshold, so no core reduction was run. The
ordered cubic witness digest was
`264234a9932531de7a5cb02c302fe9db4150548e32cf8ca84401f4237616f4d5`.
This is an exact but sampled sufficient certificate; rows left
unresolved may still have private cubic terms outside the sample.

The planted `[0,2,4,6,8]` control checked all 29,880 generated
products as zero and replayed the full group sum. Both native release
processes exited zero and the RSS wrapper did not kill either. Peak
process-reported RSS for the ordinary run was 491,814,912 bytes.
The focused test for counting a candidate in original equations passed.
The release helper source SHA-256 was
`c986271f6717bf63f4fe0325a26afd0ac65ccc41d6938823d21d05ba4e548771`;
probe source SHA-256 was
`42db2ec7331b14e2f35acf267e076dafcc0408a50d7d867c11dfda10afde8fa0`;
frozen release binary SHA-256 was
`ab037db2ca74f8205346ab2b4a7631d83ae25b85ba7b048e9eaf0a776eb7211e`.
The [runner](run.sh), [status](status.tsv), [RSS samples](rss.tsv),
raw JSONL/stderr, build/test logs and [SHA-256 manifest](SHA256SUMS)
preserve the run.

This does not establish that the remaining core has no source-only
consequences. There is no ordinary decomposition, complete F6 solver,
one-target IC online interval, recovered target or paired rho run.
The complete candidate ID and speedup remain unknown. All wall-time
figures on this contended host are exploratory structural diagnostics.
