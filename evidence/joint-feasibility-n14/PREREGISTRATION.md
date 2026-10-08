# N14: joint-protocol feasibility checks requested in the #59 review / #63 (preregistered before running)

Base: main `aab7f9b` (includes #58, released SIMSOPT). Environment: scratchpad env (Python 3.12, SIMSOPT 1.11.1, vmecpp 0.8.0), one thread. Four plasma modes: m = (rbc 1,1), (rbc 2,0), (zbs 1,1), (zbs 2,0).

**P1, selected401 parity:** regenerate selected401 at ns = 401 with vmecpp from the committed input. Qualified, with a distinct identity, if:
- the ideal narrow action score matches the #35 archived 1.100894343118703e-05 to ≤1e-6 relative;
- convergence is clean (`ier_flag` 0);
- the boundary matches the input to 1e-12.

The same policy already holds for reference401: 3.4e-8 relative, plus the starter samples.

**P2, comparability:** at ns = 25, evaluate reference, selected, and each single-mode extreme ±5e-4 (8 boundaries). Report volume, aspect ratio, major radius and B0 relative to reference.
- Preregistered tolerance for the J arm: each of these within ±1% of reference.
- If a single-mode extreme exceeds it, the bound for that mode shrinks to the largest value that complies, measured linearly. Otherwise the bounds stand.

**P3, derivative reliability:** at two base points (reference and selected), central-difference derivatives of the ideal narrow score A for each mode.
- ns = 25 with h ∈ {3e-6, 1e-5, 3e-5}; ns = 51 with h ∈ {1e-5, 3e-5}; ns = 101 with h = 1e-5.
- Pass if:
  - h-sensitivity at ns = 25: max/min over h within 10% (relative to the component's magnitude);
  - resolution agreement: ns = 25 (h = 1e-5) vs ns = 101 (h = 1e-5) within 20% for every component larger than 5% of the largest component, with the same sign for all of these.
- Stop/change rule: if ns = 25 fails but ns = 51 passes against ns = 101, the search uses ns = 51 (with the budget unchanged and evaluation cost recomputed). If both fail, the J arm is not ready and the pilot stops.

**P4, weight calibration** (no weight sweep): along the Step 3 direction d = selected − reference, at t ∈ {0, 0.5, 1}:
- measure the ideal narrow score A(t) at ns = 25;
- measure the fitter's flux objective J_flux(t), with the length-headroom coils held fixed and reference flux normalization;
- report the component changes.

Preregistered weight rule: w = ΔJ_flux(1) / (A(0) − A(1)), rounded to one significant figure. This is the break-even weight at which Step 3's own plasma change pays for its coil mismatch with these coils. It is frozen before any J run. If ΔA has the wrong sign (A does not decrease), stop: the objective is unusable.

Budget ≤40 min total; output ≤200 MiB. Failures are preserved. Exploratory feasibility only; not physics.
