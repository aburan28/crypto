# Stage 30: bounded n=31 factor-base search

Stage 30 is an additive successor to a local n=31 search that reached its 900-second outer watchdog without selecting a factor base. The failed run used 2,044.192734 core-seconds, 900.010928 wall-seconds, and 294,076,416 bytes maximum process RSS. It produced no search result or recipe and supports no mathematical conclusion.

The successor preserves the exact search parameters: `K_0/F_(2^31)`, three summands, every degree-11 divisor-kernel candidate within the 2,048-abscissa guard, 256 seeded public census targets, one candidate validated by pair-table runs on two separate holdouts, 10,000 trials per holdout, and seed 29031. The only change is the outer watchdog, extended to 3,600 seconds.

The future scalar-blind target is unavailable to factor-base selection. The search uses no factor-base log labels and does not enumerate the target subgroup. Its wrapper records total core-seconds, wall time, average parallelism, process resources, and sampled aggregate process-tree RSS. Build time is retained separately.

Only a completed validated winner may produce `factor-base.json`. A timeout or failed validation is retained as incomplete and cannot launch the larger unknown-scalar run. Neither outcome establishes scaling or SOTA by itself.
