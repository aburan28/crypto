# Worker A receipt

\`\`\`json
{
  "provider": "runpod-public-clone",
  "gpuType": "NVIDIA RTX PRO 6000 Blackwell Server Edition",
  "podId": "w7w5h1nz5nh5fo",
  "commit": "c10182d4a66550e04bb41168ba76e10446874e30",
  "exitCode": "0",
  "finishedAt": "2026-10-12T21:01:18Z",
  "terminationRequestedFor": "w7w5h1nz5nh5fo"
}
\`\`\`

Host: CUDA 13.3.73, driver 595.91.07, 188 SMs, 2430 MHz.

\`\`\`text
verify reports: 300; verify dp cap: 262144; verify steps: 16; verify launches: 6
benchmark steps: 64; benchmark launches: 4; repeats: 3
PASS: shipping sigma-table verification
sigma median 15994.424 M/s range 15980.796-16014.132 M/s
finished: 1.616 M it/s, 760232 distinguished points (300 verified against the reference, 0 dropped)
\`\`\`

Timing samples: 16014.132, 15994.424, 15980.796 M/s.
