# Stage 163: shape-selected block-4 M4RI F4 elimination

Strong internal engineering and one finite public toy-PDP native-F4 M4RI improvement. It is not a full solver panel, direct-MITM improvement, full IC/rho crossover, independent reproduction, or Koblitz index-calculus SOTA.

Large dense matrices now use block-4 Method of Four Russians elimination when rows >= 128, columns >= 256, and columns <= 4*rows. Smaller or wider matrices retain streaming elimination. `PQ_F4_DISABLE_M4RI=1` is the same-binary control. Four large synthetic shapes pass exact pivot-column and canonical-row-space equality against streaming elimination.

The same-binary M4RI arm takes 54.538915 seconds versus 58.370052 for streaming, a 1.070x gain. The clean selected process takes 54.773779 wall seconds, 54.662604 core-seconds, and 915668992 bytes RSS.

Selected M4RI routes 196 matrices and reduces counted table-inclusive word XORs from 87,513,949,370 to 45,879,309,338. Direct MITM still wins by 18.26x wall and 18.27x CPU. Clean build plus run costs 223.248 wall seconds.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and independent external reproduction remain open. This does not establish a SOTA.
