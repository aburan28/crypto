# Anthropic workload federation for the autofix queue

Assessment date: 2026-10-16. Crypto snapshot: `ee9ae5a` (the starting
checkout). This records the migration plan; the workflow and Claude Console
configuration have not yet been changed or tested together.

## Scope and current state

The requested outcome is for `.github/workflows/autofix-queue.yml` to run its
Haiku and Sonnet Claude Code attempts using a GitHub-issued OIDC identity and
an Anthropic federation rule, without supplying `ANTHROPIC_API_KEY`. The
current `fix` job requires the repository secret and passes it to both
`anthropics/claude-code-action@v1` steps. The workflow runs on a schedule,
manual dispatch, and `pull_request: labeled`. It deliberately skips fork PRs
before queueing and checks out the selected same-repository PR branch in the
`fix` job. The GitHub token used to edit and push PRs is a separate credential;
Claude API federation does not replace it or the optional `AUTOFIX_GH_TOKEN`.

The `@v1` action already accepts `anthropic_federation_rule_id`,
`anthropic_organization_id`, `anthropic_service_account_id`, and optional
`anthropic_workspace_id`. It fetches and refreshes a GitHub OIDC token when
the job has `id-token: write`. No custom exchange client is needed for these
two action steps. The action uses a shared federation profile for its Claude
subprocesses; the live queue must still verify both attempts and a long run.

The screenshot supplied with the request showed an empty Issuers tab. It is a
point-in-time observation, not a live check of the current Console resources.
No real federation rule IDs, service account IDs, or successful exchange were
available for this assessment.

| Requirement | Status on 2026-10-16 | Evidence or remaining gate |
| --- | --- | --- |
| Identify the API-key dependency | Verified | `autofix-queue.yml` key guard and both action inputs. |
| Identify an OIDC-capable caller | Verified from upstream source | The action's `@v1` inputs and Anthropic's GitHub Actions guide. |
| Restrict the trust rule to the queue | Planned | Configure a rule and inspect actual OIDC claims from a trusted run. |
| Replace the key in both queue attempts | Not implemented | Workflow change after Console resources exist. |
| Prove live exchange and queue behavior | Not attempted | Run scheduled/manual and labeled-PR paths, inspect Console authentication history and job results. |
| Remove the repository API key | Not attempted | Do this after every consumer of that secret has been inventoried and migrated. |

## Recommended design

Use **GitHub Actions as the issuer for the autofix queue**. In the Claude
Console, choose **Connect workload → GitHub Actions**. Register issuer
`https://token.actions.githubusercontent.com` using OIDC discovery. Create a
dedicated Anthropic service account in the intended workspace and a rule with
audience `https://api.anthropic.com`, a workspace binding, and
`workspace:developer` scope for Claude Code. Record the resulting `fdrl_`,
organization UUID, `svac_`, and `wrkspc_` identifiers as repository variables
or directly in the workflow. They are identifiers rather than API secrets.
Keep the service account dedicated to this queue so its quota, access, and
authentication history are attributable to it.

**Run the credentialed `fix` job from the default branch for all three entry
points.** Schedule and `workflow_dispatch` already do that. For a labeled PR,
let the lightweight intake job use its GitHub token to call
`workflow_dispatch` on the default branch with the PR number; then triage and
fix run in that dispatched workflow. GitHub documents `workflow_dispatch` as
an exception that starts a new run even when called with `GITHUB_TOKEN`.
The intake job will need permission to dispatch Actions workflows; it must
not receive `id-token: write` or run the Claude action. Validate the PR number
again in the default-branch run before selecting a branch.
Preserve the current same-repository PR check, the quiet-period behavior, and
the existing concurrency guard. This route lets one rule match the exact
default-branch subject
`repo:aburan28/crypto:ref:refs/heads/main`, the exact audience, and the
`workflow_ref` for `.github/workflows/autofix-queue.yml@refs/heads/main`.
Confirm these claims from a real token payload before saving the rule. Decode
claims only; never print or store the signed JWT itself.

The existing labeled-PR run has a different OIDC subject,
`repo:aburan28/crypto:pull_request`, and runs from a PR merge ref. Granting a
rule to `repo:aburan28/crypto:*` just to cover it would also admit other
repository contexts unless further claim restrictions were correct. Keeping
the credentialed job on the default branch makes the rule simpler to audit.
This routing change is part of the migration, not a completed change.

After the rule exists, change only the `fix` job's permissions to include
`id-token: write`. Replace the `HAS_KEY`/`Require ANTHROPIC_API_KEY` guard with
a fail-closed check that the four federation identifiers are configured. In
both Claude action steps remove `anthropic_api_key` and pass the federation
IDs. Retain `github_token: ${{ github.token }}`: that input authorizes GitHub
operations and does not authenticate to the Claude API. Ensure neither
`ANTHROPIC_API_KEY` nor `CLAUDE_CODE_OAUTH_TOKEN` is exported in the job; the
action and SDK give static credentials precedence over federation, including
some empty-value configurations.
Anthropic's action documentation also says inline comment classification is
skipped under federation; verify that this queue does not depend on that
feature before cutover.

Keep the existing API key available outside the job only during cutover, then
delete the repository secret after verifying no other workflow uses it. Do
not silently fall back to the API key when OIDC exchange fails: the job should
fail with a visible authentication error and be diagnosed through the Console
authentication history. A failed exchange does not demonstrate a failure of
OIDC as a mechanism; inspect `iss`, `aud`, `sub`, claim match, service account
membership, scope, and token replay reason.

## Verification and rollout gates

1. Confirm the Claude Console organization and workspace where the queue is
   billed. Create and verify the issuer, service account, and rule with a
   Console administrator. Bind the rule to the queue's exact workload claims.
2. In a short, authorized main-branch `workflow_dispatch` run, verify one
   successful exchange in Console authentication history and one Claude API
   request attributed to the intended service account and workspace. Do not
   log the JWT, Anthropic access token, or old API key.
3. Run the queue's manual and scheduled paths, a same-repository labeled PR,
   both Haiku-first and Sonnet paths, and one run long enough to exercise
   refresh. Confirm the labeled event dispatches a default-branch run once,
   the existing PR filters still apply, and jobs fail closed without WIF IDs.
   The action's GitHub OIDC tokens carry `jti`; a reused assertion is rejected,
   so check for `jti_reused` and refresh errors in authentication history.
4. After the new workflow has passed applicable CI and review, remove the
   `ANTHROPIC_API_KEY` repository secret if an inventory confirms no other
   consumer. Keep the current secret until that check is complete, but never
   pass it to the WIF action during the pilot.

The documentation assessment can merge independently. The queue must not be
called migrated until the Console configuration, workflow change, and live
end-to-end runs have all passed.

## Other Anthropic callers

`tools/delphi` still requires `ANTHROPIC_API_KEY` for direct batch requests.
The separate `aburan28/conductor` checkout at `c6ecf61` has key-bound Delphi
and hosted inference paths. A local `codex/anthropic-wif` Conductor branch at
`88802e0` contains a
`conductor wif exec` launcher for Claude Code using GCP, AWS, or an issuer-run
token file, but that code is absent from the checked Conductor checkout and
its live Anthropic exchange was unverified. That launcher addresses Claude
processes on identity-bearing runners; it does not migrate the Conductor
server's batch-provider or vault-secret paths. Those consumers need their own
credential resolution, token refresh, and service-account tests before their
keys can be retired. A developer's local machine also needs a real issuer
identity; the Console does not turn a normal local shell into one.

## Sources

- Local workflow: `.github/workflows/autofix-queue.yml` and
  `.github/scripts/autofix_triage.py` at the assessment revision.
- [Anthropic WIF setup](https://platform.claude.com/docs/en/manage-claude/workload-identity-federation)
  and [GitHub Actions guide](https://platform.claude.com/docs/en/manage-claude/wif-providers/github-actions).
- [Anthropic WIF reference](https://platform.claude.com/docs/en/manage-claude/wif-reference)
  for credential precedence, scope, claim matching, and `jti` replay.
- [Claude Code action `@v1` inputs](https://github.com/anthropics/claude-code-action/blob/v1/action.yml)
  and [federation setup](https://github.com/anthropics/claude-code-action/blob/v1/docs/setup.md).
- [GitHub OIDC claims](https://docs.github.com/en/actions/reference/security/oidc)
  and [`GITHUB_TOKEN` dispatch behavior](https://docs.github.com/en/actions/concepts/security/github_token).
