# Native indexed-envelope discovery readback

This is a complete, resource-qualified **Boolean matrix-construction stage** comparison. It does not measure a complete index-calculus solve or a rho boundary. Full-IC and rho costs remain null.

All 128 cells were complete: 8960 A/B observations, 2560 A/A observations, and 244800 oracle-verified output matrices. Per-call cold totals include compilation, every application, validation, output allocation and destruction.

| Candidate | Dramatic groups | Incremental groups | Full primary gate |
|---|---:|---:|---|
| envelope_degree | 0/32 | 0/32 | discovery only |
| envelope_indexed | 0/32 | 16/32 | discovery only |

The primary ratio for each group uses the pointwise fastest old control, a paired 95% median-ratio interval, and the A/A noise floor. Every declared group is retained below.

| Cell | Candidate | Median ratio | 95% lower | 95% upper | A/A floor | >2x gate |
|---|---|---:|---:|---:|---:|---|
| n6-discovery-20261002-coefficients-b16 | envelope_degree | 0.7648 | 0.7605 | 0.7716 | 1.0453 | FAIL |
| n6-discovery-20261002-coefficients-b16 | envelope_indexed | 0.8504 | 0.8407 | 0.8838 | 1.0453 | FAIL |
| n6-discovery-20261002-coefficients-b64 | envelope_degree | 0.9731 | 0.9663 | 0.9857 | 1.0287 | FAIL |
| n6-discovery-20261002-coefficients-b64 | envelope_indexed | 1.3011 | 1.2838 | 1.3092 | 1.0287 | FAIL |
| n6-discovery-20261002-degree_cycle-b16 | envelope_degree | 0.8249 | 0.8207 | 0.8287 | 1.0559 | FAIL |
| n6-discovery-20261002-degree_cycle-b16 | envelope_indexed | 0.9200 | 0.9013 | 0.9266 | 1.0559 | FAIL |
| n6-discovery-20261002-degree_cycle-b64 | envelope_degree | 0.9911 | 0.9783 | 0.9996 | 1.0512 | FAIL |
| n6-discovery-20261002-degree_cycle-b64 | envelope_indexed | 1.2443 | 1.2312 | 1.2693 | 1.0512 | FAIL |
| n6-discovery-3407879-coefficients-b16 | envelope_degree | 0.7571 | 0.7526 | 0.7639 | 1.0299 | FAIL |
| n6-discovery-3407879-coefficients-b16 | envelope_indexed | 0.8433 | 0.8353 | 0.8759 | 1.0299 | FAIL |
| n6-discovery-3407879-coefficients-b64 | envelope_degree | 0.9777 | 0.9634 | 0.9937 | 1.0250 | FAIL |
| n6-discovery-3407879-coefficients-b64 | envelope_indexed | 1.2878 | 1.2713 | 1.2991 | 1.0250 | FAIL |
| n6-discovery-3407879-degree_cycle-b16 | envelope_degree | 0.8187 | 0.8143 | 0.8267 | 1.0480 | FAIL |
| n6-discovery-3407879-degree_cycle-b16 | envelope_indexed | 0.9068 | 0.8940 | 0.9240 | 1.0480 | FAIL |
| n6-discovery-3407879-degree_cycle-b64 | envelope_degree | 0.9759 | 0.9646 | 0.9834 | 1.0286 | FAIL |
| n6-discovery-3407879-degree_cycle-b64 | envelope_indexed | 1.2522 | 1.2500 | 1.2614 | 1.0286 | FAIL |
| n8-discovery-20261002-coefficients-b16 | envelope_degree | 0.7326 | 0.7295 | 0.7377 | 1.0582 | FAIL |
| n8-discovery-20261002-coefficients-b16 | envelope_indexed | 0.8193 | 0.8131 | 0.8259 | 1.0582 | FAIL |
| n8-discovery-20261002-coefficients-b64 | envelope_degree | 0.8950 | 0.8848 | 0.8988 | 1.0358 | FAIL |
| n8-discovery-20261002-coefficients-b64 | envelope_indexed | 1.1833 | 1.1743 | 1.1967 | 1.0358 | FAIL |
| n8-discovery-20261002-degree_cycle-b16 | envelope_degree | 0.7909 | 0.7880 | 0.7934 | 1.0453 | FAIL |
| n8-discovery-20261002-degree_cycle-b16 | envelope_indexed | 0.8686 | 0.8626 | 0.8795 | 1.0453 | FAIL |
| n8-discovery-20261002-degree_cycle-b64 | envelope_degree | 0.9320 | 0.9231 | 0.9435 | 1.0376 | FAIL |
| n8-discovery-20261002-degree_cycle-b64 | envelope_indexed | 1.1829 | 1.1677 | 1.1885 | 1.0376 | FAIL |
| n8-discovery-3407879-coefficients-b16 | envelope_degree | 0.7334 | 0.7300 | 0.7395 | 1.0452 | FAIL |
| n8-discovery-3407879-coefficients-b16 | envelope_indexed | 0.8137 | 0.8062 | 0.8226 | 1.0452 | FAIL |
| n8-discovery-3407879-coefficients-b64 | envelope_degree | 0.8962 | 0.8904 | 0.9008 | 1.0376 | FAIL |
| n8-discovery-3407879-coefficients-b64 | envelope_indexed | 1.1527 | 1.1476 | 1.1603 | 1.0376 | FAIL |
| n8-discovery-3407879-degree_cycle-b16 | envelope_degree | 0.7840 | 0.7797 | 0.7866 | 1.0498 | FAIL |
| n8-discovery-3407879-degree_cycle-b16 | envelope_indexed | 0.8589 | 0.8559 | 0.8682 | 1.0498 | FAIL |
| n8-discovery-3407879-degree_cycle-b64 | envelope_degree | 0.9343 | 0.9222 | 0.9390 | 1.0463 | FAIL |
| n8-discovery-3407879-degree_cycle-b64 | envelope_indexed | 1.1708 | 1.1615 | 1.1790 | 1.0463 | FAIL |
| n10-discovery-20261002-coefficients-b16 | envelope_degree | 0.7004 | 0.6934 | 0.7050 | 1.0568 | FAIL |
| n10-discovery-20261002-coefficients-b16 | envelope_indexed | 0.7611 | 0.7541 | 0.7689 | 1.0568 | FAIL |
| n10-discovery-20261002-coefficients-b64 | envelope_degree | 0.8629 | 0.8554 | 0.8680 | 1.0267 | FAIL |
| n10-discovery-20261002-coefficients-b64 | envelope_indexed | 1.1340 | 1.1291 | 1.1448 | 1.0267 | FAIL |
| n10-discovery-20261002-degree_cycle-b16 | envelope_degree | 0.7469 | 0.7453 | 0.7488 | 1.0486 | FAIL |
| n10-discovery-20261002-degree_cycle-b16 | envelope_indexed | 0.8106 | 0.8077 | 0.8162 | 1.0486 | FAIL |
| n10-discovery-20261002-degree_cycle-b64 | envelope_degree | 0.9156 | 0.9076 | 0.9176 | 1.0193 | FAIL |
| n10-discovery-20261002-degree_cycle-b64 | envelope_indexed | 1.1695 | 1.1584 | 1.1756 | 1.0193 | FAIL |
| n10-discovery-3407879-coefficients-b16 | envelope_degree | 0.6949 | 0.6931 | 0.6980 | 1.0553 | FAIL |
| n10-discovery-3407879-coefficients-b16 | envelope_indexed | 0.7510 | 0.7472 | 0.7561 | 1.0553 | FAIL |
| n10-discovery-3407879-coefficients-b64 | envelope_degree | 0.8620 | 0.8594 | 0.8653 | 1.0370 | FAIL |
| n10-discovery-3407879-coefficients-b64 | envelope_indexed | 1.1512 | 1.1472 | 1.1543 | 1.0370 | FAIL |
| n10-discovery-3407879-degree_cycle-b16 | envelope_degree | 0.7454 | 0.7430 | 0.7502 | 1.0440 | FAIL |
| n10-discovery-3407879-degree_cycle-b16 | envelope_indexed | 0.8147 | 0.8115 | 0.8164 | 1.0440 | FAIL |
| n10-discovery-3407879-degree_cycle-b64 | envelope_degree | 0.9071 | 0.8989 | 0.9215 | 1.0249 | FAIL |
| n10-discovery-3407879-degree_cycle-b64 | envelope_indexed | 1.1740 | 1.1603 | 1.1831 | 1.0249 | FAIL |
| n12-discovery-20261002-coefficients-b16 | envelope_degree | 0.6664 | 0.6648 | 0.6701 | 1.0532 | FAIL |
| n12-discovery-20261002-coefficients-b16 | envelope_indexed | 0.7124 | 0.7117 | 0.7174 | 1.0532 | FAIL |
| n12-discovery-20261002-coefficients-b64 | envelope_degree | 0.8514 | 0.8470 | 0.8538 | 1.0227 | FAIL |
| n12-discovery-20261002-coefficients-b64 | envelope_indexed | 1.1440 | 1.1410 | 1.1505 | 1.0227 | FAIL |
| n12-discovery-20261002-degree_cycle-b16 | envelope_degree | 0.7236 | 0.7184 | 0.7265 | 1.0507 | FAIL |
| n12-discovery-20261002-degree_cycle-b16 | envelope_indexed | 0.7760 | 0.7748 | 0.7787 | 1.0507 | FAIL |
| n12-discovery-20261002-degree_cycle-b64 | envelope_degree | 0.8960 | 0.8866 | 0.8983 | 1.0214 | FAIL |
| n12-discovery-20261002-degree_cycle-b64 | envelope_indexed | 1.1672 | 1.1573 | 1.1722 | 1.0214 | FAIL |
| n12-discovery-3407879-coefficients-b16 | envelope_degree | 0.6568 | 0.6551 | 0.6644 | 1.0663 | FAIL |
| n12-discovery-3407879-coefficients-b16 | envelope_indexed | 0.6976 | 0.6927 | 0.7019 | 1.0663 | FAIL |
| n12-discovery-3407879-coefficients-b64 | envelope_degree | 0.8536 | 0.8477 | 0.8561 | 1.0290 | FAIL |
| n12-discovery-3407879-coefficients-b64 | envelope_indexed | 1.1457 | 1.1377 | 1.1524 | 1.0290 | FAIL |
| n12-discovery-3407879-degree_cycle-b16 | envelope_degree | 0.7144 | 0.7062 | 0.7214 | 1.0520 | FAIL |
| n12-discovery-3407879-degree_cycle-b16 | envelope_indexed | 0.7659 | 0.7624 | 0.7716 | 1.0520 | FAIL |
| n12-discovery-3407879-degree_cycle-b64 | envelope_degree | 0.8963 | 0.8909 | 0.8973 | 1.0199 | FAIL |
| n12-discovery-3407879-degree_cycle-b64 | envelope_indexed | 1.1600 | 1.1551 | 1.1694 | 1.0199 | FAIL |

Host: aarch64 on 4 logical CPUs; reserved CPUs Array [Number(3)]. Whole-campaign wall time 12.645 s, other-process CPU 0.220 s, peak worker RSS 13848 KiB. The resource receipt covers the whole campaign; individual short calls have separate timers and A/A pairs, not individual resource certifications.

The frozen raw JSONL, input source, binaries, host records, resource receipt, native verifier result and this report are bound by the artifact manifest. A positive full result requires a further unchanged-source run on unused seeds.
