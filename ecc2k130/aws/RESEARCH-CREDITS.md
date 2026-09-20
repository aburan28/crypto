# Cheap operation and research credits

This campaign exists to **measure index-calculus algorithms** against
a generic-group bound, so post-quantum parameter work has a number
instead of a phase quote. It is mathematical science. It is not a
key-recovery service.

Cost is iterations per dollar, not dollars per hour. The unit is the
same one the scoreboard uses: completed scalar updates. Wall-clock
belongs here as a practicality note, never as the metric.

Expected work remains \(2^{60.9} = 2.15 \times 10^{18}\) iterations
(Bailey et al.). That is an expectation, not a deadline.

## What the live fleet already is

Snapshot 2026-09-18 09:14 UTC, account `590183823895`:

- Campaign GPUs were **all spot**. The on-demand g7/g7e base from the
  first launch is gone.
- Costguard's $5k/month stop still prices a handful of tagged boxes
  against on-demand list. Most leftover spots are `CostGuardExempt` and
  do not hit that brake.
- A leftover `g7e.2xlarge` bench (`ecc2k130-iterfn-bench-g7e-r3`) and a
  `crypto-g7e-cursor` on-demand box are not walking the campaign and are
  the first things to stop.

## Iterations per dollar (matched spot quotes)

Rates are the tree's receipts: g7e collecting 14.1 B/s, g7 `--bench`
5.107 B/s, g6e CLMAD `--bench` 8.838 B/s, g6 CLMAD `--bench` 2.490 B/s,
g4dn software `--bench` 0.533 B/s. Collecting on Ada/Turing is still
unmeasured; if DP handling costs the same ~4% it does on the 6000, the
ranking below does not move.

Spot quotes taken in one pass at 2026-09-18T09:14Z
(`describe-spot-price-history`, cheapest AZ per type):

| type | GPU | cheapest $/h | it/$ (×10¹²) | region of that quote |
|---|---|---:|---:|---|
| g6.2xlarge | L4 | 0.1602 | **56.0** | eu-north-1 |
| g6e.2xlarge | L40S | 0.8812 | 36.1 | eu-north-1 |
| g7e.2xlarge | RTX PRO 6000 | 1.3271 | 38.2 | eu-north-1 |
| g7.2xlarge | RTX PRO 4500 | 0.6627 | 27.7 | us-east-1 |
| g6.2xlarge | L4 | 0.3194 | 28.1 | eu-west-3 (Paris, currently full) |
| g4dn.xlarge | T4 | 0.0894 | 21.5 | sa-east-1 |
| g4dn.2xlarge | T4 | 0.1098 | 17.5 | eu-north-1 |
| g4dn.xlarge | T4 | 0.2470 | 7.8 | ca-central-1 (leftover fill) |

So:

1. **Cheapest progress is leftover g6 in a cheap region** (Stockholm,
   São Paulo), then g7e / g6e where those spots exist.
2. **T4 leftover is the wrong spend.** Canada xlarge T4s buy about
   7.8×10¹² it/$ against 56×10¹² for a Stockholm L4. `CHEAP=1` on
   `launch_spot_all.sh` stops filling Turing.
3. **Do not buy a Compute Savings Plan or a capacity reservation** for
   this campaign. Promotional credit does not cover SP/CR upfront fees,
   spot dies on interrupt, and a 1-year commit is a bet the thread has
   not earned.
4. On-demand g7e at list (~$3.36/h in us-west-2) is ~15×10¹² it/$ —
   worse than every spot row above. Do not refill OD.

At the Stockholm g6 quote a full generic-group reference walk of the
Bailey expected work is about **$38k**. That number is a *ceiling on
the reference*, not a budget for a key-recovery job. The index-calculus
thread is priced in the same unit so a crossover, if one exists, can
be read off one table.

## Modal vs GCP (same unit, same L4 / 6000 receipts)

**GCP is cheaper than Modal. Leftover AWS g6 is cheaper than both.**
Neither is a reason to move the live campaign: the binary, S3 corpus,
and bootstrap already sit on account `590183823895`. Use Modal for
short benches and the academic credit; use GCP only if a research
award lands that cannot be spent on AWS.

Rates are the same receipts as the AWS table. Prices are public list
quotes, 2026-09-18. Modal GPU-only from
[modal.com/pricing](https://modal.com/pricing) (`$/s × 3600`). GCP
rows are `g2-standard-*` spot (VM + one L4 together) from published
aggregator list prices, not a receipt on our project.

| host | GPU | list $/h | it/$ (×10¹²) | note |
|---|---|---:|---:|---|
| AWS leftover g6.2xlarge | L4 | 0.1602 | **56.0** | Stockholm spot; our quote |
| GCP g2-standard-4 spot | L4 | 0.192 | 46.7 | cheapest published APAC |
| GCP g2-standard-8 spot | L4 | 0.229 | 39.2 | APAC; 8 vCPU, closer to g6.2xlarge |
| GCP g2-standard-8 spot | L4 | 0.340 | 26.3 | US-East matched SKU |
| AWS leftover g7e.2xlarge | RTX PRO 6000 | 1.327 | 38.2 | Stockholm spot; our quote |
| Modal RTX PRO 6000 | RTX PRO 6000 | 3.031 | 16.8 | $0.000842/s, GPU only |
| Modal L40S | L40S | 1.951 | 16.3 | $0.000542/s, GPU only |
| Modal L4 | L4 | 0.799 | 11.2 | $0.000222/s, GPU only |
| Modal L4 + 1 core + 8 GiB | L4 | 0.910 | 9.8 | same L4 plus Modal CPU/memory |

`it/$ = (B/s) × 3600 / ($/h) / 10¹²` with L4 at 2.490 B/s, L40S at
8.838 B/s, collecting 6000 at 14.1 B/s. Modal region multipliers
(1.15–1.75× on the Team plan) make those rows worse, not better.

So:

1. **On the same L4, GCP spot is about 2–4× Modal.** The gap is the
   serverless markup, not the silicon.
2. **On the 6000, leftover AWS g7e is about 2× Modal.** Do not rent
   Modal 6000s to walk.
3. **Do not migrate the campaign to GCP for a $1–5k research credit.**
   Standing up a second corpus eats the award. The cheap route is
   leftover AWS g6 plus the AWS Cloud Credit for Research form, or
   Modal's academic grant (up to $10k, faster review) spent on
   *benches*, not on the walk.
4. Modal stays the right host for a one-shot `--bench` / TOP_CLMAD
   receipt: per-second billing, no AMI, no leftover OD box.

## What this research is (and is not)

This tree studies **bleeding-edge index calculus** on algebraic groups
— Semaev-style summation polynomials, Weil-descent encodings, Gröbner /
WDSat solving, factor-base geometry — in order to say where those
methods actually sit relative to a generic-group algorithm. The
audience is post-quantum cryptography: if an algebraic method has a
real crossover, parameter writers need that number; if it does not,
they need that number too. The rule in `AGENTS.md` is the scientific
one: every claimed improvement is a **ratio to a boundary**, in one
unit, with a verified answer. A thread that cannot state its floor
has not started.

Pollard rho on a Koblitz curve is the **reference**, not the product.
It is run so the index-calculus variants can be placed on the same
axis (\(S = \mathrm{ops}/\sqrt{n}\)). It is not a service for
recovering keys, and a credit proposal must not read as one.

That framing is also the one the AWS program will accept. The FAQ
funds finite proofs of concept, benchmarks, and shareable research
tools. It refuses ongoing operations and “general funding of a lab.”

## Best credit route

Ranked for this scientific workload, not in general.

### 1. AWS Cloud Credit for Research — only if you have a `.edu` / institution address

<https://aws.amazon.com/government-education/research-and-technical-computing/cloud-credit-for-research/>

Apply at <https://pages.awscloud.com/aws-cloud-credit-for-research.html>.
The filled paste-ready form is
[`AWS-RESEARCH-CREDIT-APPLICATION.md`](AWS-RESEARCH-CREDIT-APPLICATION.md).

| | |
|---|---|
| Who | Full-time faculty or research staff (uncapped), or enrolled graduate / PhD student (**$5k cap**) |
| Email | Institution-issued only. Personal Gmail is rejected. |
| Account | Paid-tier AWS account number. As of 2026-02-16 Free Tier accounts cannot receive promotional credit. |
| What it pays | Promotional credit on **on-demand and spot** EC2. Not Savings Plan / reservation upfront. |
| Review | Rolling, **90–120 days** on the program page (some mirrors say 30–60). No expedite. |
| Contact | `aws-research-credit@amazon.com` |

**Apply as mathematical cryptanalysis in service of PQC.** Title it
as a finite, published measurement of algebraic index-calculus
against a generic-group bound, with open tooling other labs can
rerun. That is the program's first project type (“proof of concept
or benchmark for comparison”) plus the second (“repeatable, sharable
solutions”).

Proposal skeleton (paste-ready):

1. **Problem.** Post-quantum parameter selection needs to know whether
   algebraic index-calculus on small-characteristic or extension-field
   elliptic curves ever beats a generic-group algorithm, and at what
   size. The literature often prices one phase (relations, or Gröbner)
   and treats that as the method. This project prices the *whole*
   method — precomputation, encoding, solving, linear algebra, lifting,
   verification — against a counting floor and against Pollard rho on
   the same instances, in one unit.
2. **Scientific goal.** Locate the crossover (or show there is none
   in the sizes that fit) for the index-calculus variants already
   built in this repository, and publish every row so a later
   algebraic improvement can be classified as an advance, engineering,
   relabelling, or accounting. The work informs PQC: it says which
   algebraic attacks actually move the bound that a next-generation
   scheme has to clear.
3. **What it is not.** It is not a key-recovery service and not a
   campaign to finish a specific cryptanalytic challenge. Distinguished
   points and planted secrets exist so every row can be *checked*. A
   run without a verified answer is not a result.
4. **AWS work.** EC2 G/VT spot (g6 / g6e / g7 / g7e) for (a) the
   index-calculus solver stages and (b) a matched rho *reference* on
   the same curves, so \(S\) is comparable. S3 holds frozen experiment
   directories. Timeline: 90 days, four field sizes, then stop. The
   rho reference is sized to settle the ratio, not to exhaust
   \(2^{60.9}\) iterations.
5. **Share.** `docs/index-calculus-scoreboard.html`, the frozen
   experiment trees under `research/`, `AGENTS.md`'s boundary/table/
   ratio rule, and this repository. No private corpus. Anyone can
   rerun a row.
6. **After the credit.** The scoreboard and the method stay public.
   Future use is other researchers adding a variant as a new row
   against the same floor. The reference walk does not have to
   continue.
7. **Pricing Calculator URL.** Size a 90-day *measurement*, not a
   challenge attempt: e.g. 32× g6.2xlarge spot for IC solver +
   reference rho at the sizes that fit, plus a small g7e slice for
   the Blackwell receipt. Do not price 252× g7e.
8. **Student vs faculty.** A student award ($5k) covers the solver
   regression and a few reference sizes. A faculty/staff application
   can cover the full four-size exponent fit. If you are a student,
   put a faculty PI on the form and ask for the measurement slice.

Account to put on the form: `590183823895` (already paid-tier, already
has the experiment bucket). Do not open a second account just to apply.

### 2. AWS Activate — if this is a startup, not a thesis

<https://aws.amazon.com/startups/credits/>

- Founders / self-funded: **$1k immediately**, select participants up
  to $5k. Days, not months.
- Portfolio (accelerator / angel Org ID, pre-Series B): up to **$200k**.
- Paid-tier account, founded < 10 years, not a Free Tier account.

$1–5k is the fastest cash if the applicant is a company. $200k
covers a larger measurement campaign without a faculty appointment.
It is the wrong form if you cannot name an Activate Provider, and
the wrong form if the honest description is a thesis.

### 3. What not to chase for this thread

| program | why not |
|---|---|
| NVIDIA Academic Grant | Current CFPs are robotics / AV / 5G / federated learning, delivered as Saturn Cloud H100 hours, not EC2 T4/L4/Blackwell walking this client. |
| Google Cloud research credits | $1–5k. GCP L4 spot is cheaper than Modal L4 (~26–47 vs ~11 ×10¹² it/$) and still below leftover AWS g6 (56). Migration cost eats the award. |
| AWS Educate | Teaching credits, not a 90-day measurement. |
| Capacity reservations / Savings Plans | Upfront, 1-year, and promotional credit will not pay the reservation fee. |

## What to do this week (no credit in hand)

1. `CHEAP=1` leftover fills only. Leave Turing quota empty rather than
   spend it at 8×10¹² it/$.
2. Prefer regions whose g6/g6e/g7e spot is actually cheap (eu-north-1,
   sa-east-1, eu-west-3) over leftover T4 in ca-central-1 / Sydney.
3. Stop leftover bench/dev on-demand g7e when they are not measuring
   a named receipt.
4. Keep the $5k costguard as a hard stop on *tagged* boxes. Do not
   raise it while a credit application is in flight — a $5k student
   award and a $5k ceiling are the same number.
5. File the research-credit application as the 90-day index-calculus
   measurement (PQC-relevant algebraic vs generic-group bound), or
   Activate Founders if you have a company and need money this month.
   Size the calculator to the measurement, not to a challenge attempt.
