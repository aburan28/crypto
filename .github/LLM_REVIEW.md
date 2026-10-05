# Repository-owned pull-request review

`.github/workflows/llm-review.yml` replaces reliance on a hosted review bot's
paid usage pool. It reviews pull requests with OpenRouter's free-model router
and can fall back to Fireworks. Provider outages and free-tier rate limits are
reported as warnings rather than blocking correctness CI.

Configure one or both providers in the repository settings:

1. Add the Actions secret `OPENROUTER_API_KEY`. The default model is
   `openrouter/free`; set the Actions variable `OPENROUTER_MODEL` to override
   it.
2. Optionally add the Actions secret `FIREWORKS_API_KEY` and the required
   Actions variable `FIREWORKS_MODEL`. Fireworks is attempted only when both
   are present, after OpenRouter fails or when OpenRouter is not configured.

The workflow runs under `pull_request_target` so that provider credentials are
available for fork pull requests. It checks out only the trusted default branch,
downloads the proposed patch through GitHub's API, sends that patch to the model
as untrusted text, and never executes pull-request code. The bot maintains one
comment per pull request instead of adding a new comment on every push.

Cursor Bugbot is a separate GitHub App and cannot be redirected to either API
from repository configuration. Disable its automatic reviews in the Cursor
dashboard after this workflow is enabled if duplicate reviews are undesirable.
