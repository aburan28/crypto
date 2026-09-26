# Pre-outcome release hold

The first main-based #804 release at `96fd75a59eb413f6dce995c4fd5cdc8c00b77855` pinned main `6c69273aa1b5011897ea3dd493a3ce00d55bd20a`. Before any measured child started, unrelated merged #829 advanced main to `66fde6a1e6990358e00b72a0cef54c9fe1c230d3`; the exact-main gate would reject it. `HOLD.json` records the changed paths and seals the prior freeze and protocol bytes. This is a superseded pre-outcome release, not a toy, solver, proof, or n131 result. The initial corrected-interface failure remains separately archived in `preoutcome_failure_0`.
