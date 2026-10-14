# Worker B receipt

\`\`\`json
{
  "provider": "runpod-public-clone",
  "gpuType": "NVIDIA RTX PRO 6000 Blackwell Server Edition",
  "podId": "q43dhtqgu21c9w",
  "commit": "c10182d4a66550e04bb41168ba76e10446874e30",
  "exitCode": "0",
  "finishedAt": "2026-10-12T21:03:55Z",
  "terminationRequestedFor": "q43dhtqgu21c9w"
}
\`\`\`

Host: CUDA 13.3.73, driver 580.173.02, 188 SMs, 2430 MHz.

\`\`\`text
verify reports: 300; verify dp cap: 262144; verify steps: 16; verify launches: 6
benchmark steps: 64; benchmark launches: 4; repeats: 3
PASS: shipping sigma-table verification
sigma median 15760.478 M/s range 15739.614-15783.810 M/s
finished: 1.512 M it/s, 760232 distinguished points (300 verified against the reference, 0 dropped)
\`\`\`

Timing samples: 15739.614, 15783.810, 15760.478 M/s.
