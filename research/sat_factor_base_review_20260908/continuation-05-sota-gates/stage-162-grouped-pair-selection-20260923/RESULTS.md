# Stage 162: grouped exact F4 critical-pair selection

Strong internal engineering and one finite public toy-PDP native-F4 pair-selection improvement. It is not a full solver panel, direct-MITM improvement, full IC/rho crossover, independent reproduction, or Koblitz index-calculus SOTA.

A charged post-hash profile put 695 of 3,327 top-of-stack samples in `State::insert`. Disassembly located the dominant loop in the quadratic scan that asks whether any other new-pair LCM divides the current LCM. Stage 162 groups equal LCMs, finds proper divisor LCMs by exact submask lookup, and then emits the same lowest-index representative in the same order as the prior UPDATE algorithm.

The clean selected process takes 57.862894 wall seconds, 57.736423 core-seconds, and 968146944 bytes peak RSS. Whole-target wall improves 1.048x over Stage 161. Equations, term count, F4 calls, matrices, pair counters, 87,513,949,370 elimination word XORs, roots, and witness are identical.

Direct MITM still wins by 19.29x wall, 19.30x CPU, and 21.12x RSS. Clean build plus run costs 225.354 wall seconds and 221.602 core-seconds.

Licensed Magma, a full native-F4 panel, end-to-end IC/rho cost, and independent reproduction remain open. This does not establish a SOTA.
