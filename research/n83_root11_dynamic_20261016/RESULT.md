# N83 dynamic one-root index: 4x smaller exact table with replayed one-point result

The dynamic five-byte one-root table retained exact four-summand recovery while reducing K1271 slot allocation from **5 GiB to 1.25 GiB**. On the frozen public point, all seven IC arms and the same-point rho run independently replayed. Both mirrored cold-process controls passed the frozen 10% A/A gate: fixed K1271 differed by 3.55%, dynamic K1271 by 0.009%. The complete screen status is `SCREEN_COMPLETE`; the selected factor base remains null until the four frozen reserve points confirm a candidate.

| Arm | Usable base points | Cold process, s | One-target online, ms | Exact root slots, GiB | Peak RSS, GB | Independent replay |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| K1271 fixed A | 210,986 | 158.352 | 262.842 | 5.00 | 6.26 | pass |
| K1271 dynamic A | 210,986 | 153.108 | 224.066 | 1.25 | 2.23 | pass |
| K1271 root22 | 210,986 | 152.779 | 223.362 | 1.25 | 2.97 | pass |
| K1024 dynamic | 169,984 | **151.016** | 74.419 | 0.625 | 1.27 | pass |
| K1536 dynamic | 254,976 | 170.652 | 436.327 | 1.25 | 2.61 | pass |
| K1271 dynamic B | 210,986 | 153.094 | 225.483 | 1.25 | 2.23 | pass |
| K1271 fixed B | 210,986 | 152.830 | 221.494 | 5.00 | 6.26 | pass |

K1024 dynamic has the lowest observed complete cold process interval in this allocation. K1271 dynamic averaged 153.101 s versus 155.591 s for fixed capacity, a 2.490 s difference. The fixed control itself varied by 5.522 s, so this panel establishes the memory reduction and correct dynamic behavior; it does not resolve a cold-time speed gain between the K1271 layouts. K1024 dynamic and K1271 root22 differ by 1.763 s, which also needs reserve-point confirmation. The same-point rho online interval was 633.415 ms. These one-point online intervals are descriptive until alternating pairs pass the Tier 1 gate.

The [IC1 row record](ic1_rows_38876385086.json) binds all seven verified runs to five immutable candidate manifests on workload `eebf98ee1031`, the same frozen public point as the smaller-base screen. K1024 dynamic is `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0ha513645c20c8`; K1271 dynamic is `IC1N83Ckb1a1fb210986PDP4rootRCguidedLAgaussTDdirectISO0hd2c722353afb`. The fixed and dynamic A/B controls each share their candidate ID and use distinct `R1`/`R2` run IDs. `materialize.py` checks the archived artifact, source, workload, executed table policy, exact base dump and canonical hash before writing those records. Running it again on the saved artifact is idempotent:

```sh
python3 research/n83_root11_dynamic_20261016/materialize.py \
  /path/to/extracted-artifact/results \
  --archive /path/to/n83-root11-dynamic-38876385086.tar.zst
```

The allocation was a GitHub Actions Ubuntu 24.04 runner exposing four AMD EPYC 9V74 CPUs and 15.61 GiB host memory. All completed arms met the 15.25 GiB RSS guard and recorded zero new cgroup throttle, OOM and OOM-kill events. This cloud instance has exploratory timing status; the run did not include the Tier 1 alternating five-point calibration panel. Cold process timing includes setup, rank acquisition, target handling and producer verification. Target-online timing begins after reusable base, index and factor-log preparation; independent scalar replay is retained separately.

The checked-in [analysis](analysis_38876385086.json), [executed-source receipt](source_38876385086.json) and [host receipt](host_38876385086.json) accompany [GitHub Actions run 38876385086](https://github.com/aburan28/crypto/actions/runs/38876385086), artifact ID `12036957402`. The full artifact includes stdout, stderr, exact base dumps, target relations, phase times, resource usage, source and binary hashes, and verifier output. A local copy is archived at `/Volumes/SSD990-2/n83-archive-20261016/raw/n83-root11-dynamic-38876385086.tar.zst` (31.1 MiB compressed), SHA-256 `ed24dcf04f751af6a7d832c21bb8e26295bd2cc4635eacfaec994751b3b82ff6`. The analysis SHA-256 is `91aa67e82993dc0ebf001768b0d91be77d79ba610c713d685dc1f43805c4e0a3`; the executed IC binary SHA-256 is `5287c6d1de3b71a8e2044c07db26ae31ac79cd2a42af3607469c2119a7c16d28`. The full serial release library suite passed 4,148 tests; its preserved log SHA-256 is `64a7dfe5feccc34eff85dac37e49e32eb4c13398dd299c1d5c600b55b242b266`.

Next, compare K1024 dynamic against K1271 root22 and K1271 dynamic on the four reserved public points, then freeze a separate five-point alternating IC/rho Tier 1 panel for the selected exact candidate. The first-witness equality bridge in Lean remains open as recorded in [PROOF_MAP.md](PROOF_MAP.md); the checked capacity arithmetic and independent full-base replay support this implementation result.
