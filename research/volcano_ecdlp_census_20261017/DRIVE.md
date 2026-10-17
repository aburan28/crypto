# Drive storage

Folder: [`volcano-ecdlp-census-20261017`](https://drive.google.com/drive/folders/11MPowjLRt25fg1RypJnbDDvWrNc_7PZ3) (id `11MPowjLRt25fg1RypJnbDDvWrNc_7PZ3`).

| phase | file | Drive id | bytes | SHA-256 | read back |
|:--|:--|:--|--:|:--|:--|
| A | FROZEN.json | 1_Y8auw7w4SsxeZ_Hpen7mcpyPw-RohDw | 3095 | 345028509a93be46f3b0867d33af90350a9dde001d90bbc51dd1771a5a98ed73 | size read back from the upload result = 3095 (matches); content not re-downloaded |

Staged but NOT uploaded (the only upload path in the executor is a base64 argument typed through the model context: about 120k tokens for the 94 KB specs archive and about 600k tokens for the 0.97 MB census archive, with no way to hash-check the copy beforehand). Their contents are committed in git on this branch. SHA-256 of the staged archives (built with `tar | xz -9e`, PNGs excluded from the census one; kept in the executor scratchpad only):

| archive | bytes | SHA-256 |
|:--|--:|:--|
| census-phaseA.tar.xz (census/, checks/, FROZEN.json; no PNG) | 971516 | 6d9be453cad9c63fc3cabee6e28257367486d7cbd569c26296b16001d98a9708 |
| specs-phaseA.tar.xz (specs/, runs/) | 93708 | 23ab2d128912d092063477f4a8398a11b7f7032f96e107b5199281f247c04a0d |
| A | phaseA-calibration-summary.md (= calibration/summary.md) | 1flai8c2P756HeyJFo33-dKShIHHlAfV2 | 2854 | 52ec924c415d454662a3713a4737ffffbc944a8616c8a402d89c05b3d913a6f3 | size matches |
