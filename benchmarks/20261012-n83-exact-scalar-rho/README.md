# Exact scalar candidate and expanded rho policy

The selected reference-scalar Fold256/DP6 policy is compared against the
previously tuned reference-scalar Fold64/DP6 control on five fresh public
synthetic targets, with three paired seeds each. Native current/selected
process geometric ratio is **1.069651 [1.043178, 1.089238]**; mean process
cost is **2.987008 s versus 2.780545 s**. The lower interval edge overlaps
the observed A/A ratio range. This is a conditional one-host panel result.

Exact Barrett is verified but unselected: training means 2.923128 s versus
reference 2.903055 s. All 139 processes pass independent native answer
replay; internal capped-walk restarts remain in the charged work. Pooled
15-draw process medians are current 2.833963 s versus selected 2.850072 s,
an ordering different from the five target means. All draws and losing
policies are retained.

- [Run and results PDF](n83_faster_rho_evaluation.pdf).
- [Independent arithmetic/formalization audit PDF](independent_scalar_audit.pdf).
- [SVG and 150 dpi PNG copies](figures/).
- [Authoritative study](../../research/n83_speed_20261012/).
- [Report sources, copied native inputs and QA](../../research/n83_speed_20261012/reports/run-20261012T1543Z/).

These are byte-for-byte browsing copies. The report contains both the
descriptive run record and evidence-bounded results interpretation. Rebuild
the source report, then refresh every copy; never edit a rendered copy alone.
The study supplies raw records and native algorithms. Python renders recorded
values and documents only. Host: AMD EPYC 9V74 model, KVM guest, CPU3, one OS
thread, native core reservation, no contended campaign samples. No physical
host exclusivity, L2 admission, independent-host timing replay or calibrated
operation-equivalent S is asserted.
