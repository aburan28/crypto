# B7a's cases have no documents of their own

Every B7a case copies a document an earlier step froze (`../cases.json`,
`sources_sha256`). Some name it as `{here}/../../v2-b2/params/…`, a path
through this directory, so the directory must exist for the path to
resolve. This file is what keeps it in Git (B7a's protocol, amendment 1).
