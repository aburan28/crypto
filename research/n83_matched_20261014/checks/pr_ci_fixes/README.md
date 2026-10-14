# PR 1828 CI repairs

The published head `12321a1855d0afada849d601559d1b5eaa6e0b39` had seven successful PR workflow runs and two failed runs in the retained GitHub snapshot. This directory preserves the failures and the local repair checks. [receipt.json](receipt.json) records exact commands, exit codes and source hashes.

The Lean workflow failed at its prerequisite check because the hosted runner lacked ripgrep. It now installs ripgrep before the unchanged proof-source rejection gate. The exact configured gate accepts the unchanged checked proof sources, rejects each forbidden construct, and fails closed when the executable is missing. This repair check exercises the source gate; the prior compiled proof records remain in the original study.

The curve-name diff check rejected 19 historical labels in the measured source snapshot even though the checker documents frozen archives as exempt. Diff validation now applies the existing normal-lint exemption contract to complete current contents and file paths. Seven tests cover the existing frozen directories, active paths, similar directory names, exact hash and manifest locality, modified contents, retired filenames, and the prior precise archived-label exception. The retained original checker reproduces all 19 failures; the repaired checker reports zero problems on the same published-head diff.

The old checker source is retained under [archives/](archives/), along with the GitHub snapshot, job payloads and failed-job logs. Both executions of that old source remain in the receipts: its initial local path and its final archive path. Existing measured source archives and Lean sources have zero working-tree diff.

A single remote-main fetch and a dry merge established that published head `12321a1` merges cleanly with `ed3780f4`. The dry merge created only a temporary Git tree object. No merge-conflict request was posted to Copilot, and no branch, index or working-tree merge occurred.

Cursor's review still requires human review because its analysis could not cover the large PR. The auxiliary autofix job also stopped before claiming the PR because a repository service secret was absent. These remain separate from the two repaired CI defects.

Next validation gate: commit/push these repairs with the new matched-study preparation and inspect applicable CI on that resulting head. The repository's Conductor command was unavailable (exit 127); its attempted command and output are retained in the receipt. No Cargo builds or performance measurements ran for this repair.
