"""Claude Mythos 5.1 Message Batches runner.

``python3 -m tools.mythos_batch run --input requests.jsonl --out runs/x`` submits
the requests to the Anthropic Message Batches API, waits for the batch and writes
one result row per request.  Credentials come from the environment
(``ANTHROPIC_API_KEY``); see README.md.
"""

from .spec import (  # noqa: F401
    DEFAULT_MODEL,
    RequestError,
    build_params,
    load_requests,
    normalize_row,
)
from .runner import (  # noqa: F401
    collect,
    load_manifest,
    make_client,
    status,
    submit,
    wait,
)
