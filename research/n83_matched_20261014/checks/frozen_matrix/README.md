# Prepared execution matrix, before solver timing

Native controller prepare passed and froze32 ordered processes, their public
point inputs/seeds/config, the protocol and relevant source/binary/build hashes.
It launched no solver. The matrix SHA256 is
`896bd5d03e6f7d68998851aa7f9688f2fe0122546173cf74ab6b138eb19975d0`;
the machine protocol SHA256 is
`abc72a3e98367ea0ccbde40f44c68aa4e3355b43ace579f4f162355277d5b1bc`.

The independent native checker was then deliberately run on the prepared,
unexecuted matrix. It validated frozen identities/order/bindings, reports
zero invalid rows and218 expected missing answers, and returns exit2 with
`INCOMPLETE_OR_RETAINED_FAILURES`. This expected negative prerequisite check
is not a failed solver run or a successful completed campaign.

Separate old/new rho singleton equivalence preflight passed on a prior target
and fixed seed. Native independent scalar replay, nine common outcome fields,
preparation and total charges agree. In-memory altered-counter and altered-
scalar negative checks are rejected. This preflight is excluded from the new
panel's32-process matrix and its statistics. Its raw outputs/reservations/
resource commands are under `../equivalence_preflight/`.

The environment restarted before any main solver child. All five executable
hashes, source freeze and saved matrix/protocol were rechecked after startup;
files remained intact. Actual main-start and reservation records follow only
after the pre-timing metadata commit.
