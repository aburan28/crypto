# Cheap operation and research credits

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

At the Stockholm g6 quote the expected work is about **$38k**. At the
Stockholm g7e quote it is about **$56k**. A week-solve (3.55 T it/s)
is the same dollars over 168 hours, not a different budget.

## Best credit route

Ranked for this workload, not in general.

### 1. AWS Cloud Credit for Research — only if you have a `.edu` / institution address

<https://aws.amazon.com/government-education/research-and-technical-computing/cloud-credit-for-research/>

Apply at <https://pages.awscloud.com/aws-cloud-credit-for-research.html>.

| | |
|---|---|
| Who | Full-time faculty or research staff (uncapped), or enrolled graduate / PhD student (**$5k cap**) |
| Email | Institution-issued only. Personal Gmail is rejected. |
| Account | Paid-tier AWS account number. As of 2026-02-16 Free Tier accounts cannot receive promotional credit. |
| What it pays | Promotional credit on **on-demand and spot** EC2. Not Savings Plan / reservation upfront. |
| Review | Rolling, **90–120 days** on the program page (some mirrors say 30–60). No expedite. |
| Contact | `aws-research-credit@amazon.com` |

**Do not apply as “run Pollard rho until ECC2K-130 falls.”** The FAQ
refuses ongoing research, lab operations, and general lab funding. A
year-long walk is exactly that.

**Do apply as a finite, published benchmark.** The repository already
is one: generic-group rho versus algebraic index-calculus on the same
instances, one unit \(S = \mathrm{ops}/\sqrt{n}\), a public scoreboard,
and open tooling. That is the program's first project type (“proof of
concept or benchmark for comparison”) plus the second (“repeatable,
sharable solutions”).

Proposal skeleton:

1. **Problem.** Where does index-calculus on \(E(\mathbb{F}_{p^3})\)
   actually cross Pollard rho? The literature quotes phase costs; this
   tree prices the whole method.
2. **AWS work.** EC2 G/VT spot (g6 / g6e / g7 / g7e) running the
   published `ecc2k130` client; S3 for the unversioned DP corpus;
   CloudWatch/costguard for the $ ceiling. Timeline: 90 days of
   measured comparison across four field sizes, then stop.
3. **Share.** `docs/index-calculus-scoreboard.html`, frozen experiment
   directories, and this repository. No private corpus.
4. **After the credit.** The scoreboard stays; the walk does not have
   to. Future use is other researchers repeating the comparison.
5. **Pricing Calculator URL.** Size it as a 90-day *benchmark*, not a
   week-solve: e.g. 32× g6.2xlarge spot for 90 days is the right
   order, not 252× g7e.
6. **Student vs faculty.** A student award ($5k) is a few days of the
   current leftover fleet. Only a faculty/staff application can cover
   the $38k expected-work figure. If you are a student, apply for the
   benchmark slice and get a faculty PI on the form.

Account to put on the form: `590183823895` (already paid-tier, already
has the campaign bucket). Do not open a second account just to apply.

### 2. AWS Activate — if this is a startup, not a thesis

<https://aws.amazon.com/startups/credits/>

- Founders / self-funded: **$1k immediately**, select participants up
  to $5k. Days, not months.
- Portfolio (accelerator / angel Org ID, pre-Series B): up to **$200k**.
- Paid-tier account, founded < 10 years, not a Free Tier account.

$1–5k is the fastest cash. $200k is the only AWS product that can
cover a week-solve without a faculty appointment. It is the wrong
form if you cannot name an Activate Provider.

### 3. What not to chase for this thread

| program | why not |
|---|---|
| NVIDIA Academic Grant | Current CFPs are robotics / AV / 5G / federated learning, delivered as Saturn Cloud H100 hours, not EC2 T4/L4/Blackwell walking this client. |
| Google Cloud research credits | $1–5k and a different GPU stack. Migration cost eats the award. |
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
5. File the research-credit application as the 90-day benchmark, or
   Activate Founders if you have a company and need money this month.
   Then wait. Do not scale to a week-solve on a card.
