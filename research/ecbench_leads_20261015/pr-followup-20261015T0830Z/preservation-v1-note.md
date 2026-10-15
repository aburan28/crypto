# Initial local process-permission failure

The first `node .../preserve.mjs` invocation exited1. Its stdout was redirected
to preservation-v1.json (empty); its stderr was returned by the execution tool.
This note records the diagnostic, rather than representing a transcription as a
raw captured log:

    Error: spawnSync git EPERM
    errno: -1
    code: EPERM
    syscall: spawnSync git
    spawnargs: [write-tree]

The subprocess diagnostic also reported status0 and tree output
`a81367b9a6846aeefd9320660a90383dae898b43`. No scientific file was edited.
The unchanged script succeeds in the network-enabled execution context; its
captured output and status are preservation-v2.log and preservation-v2.exit.
