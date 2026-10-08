# Joint-protocol feasibility probe N14 (issues #37, #63)

Exploratory feasibility evidence for [JOINT_PROTOCOL.md](../../docs/optimization/JOINT_PROTOCOL.md);
not a physics result. Preregistered in [PREREGISTRATION.md](PREREGISTRATION.md) before running.

- Producer: upstream `main` at `aab7f9b`, script [feasibility.py](feasibility.py), clean tree.
- Command: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
  PYTHONPATH=src python feasibility.py <length-headroom seed snapshot> <fresh output dir>`; the seed is
  `coil_check.candidate_snapshot` of `submissions/length-headroom-six-coil/candidate.json`.
- Environment: Python 3.12.15, SIMSOPT 1.11.1, vmecpp 0.8.0, NumPy 2.5.3, SciPy 1.18.1,
  netCDF4 1.7.4; macOS x86_64; one thread; 1270 s.
- [result.json](result.json): every solve (tag, ns, action A, solve time, volume, aspect,
  major radius, B0, ier), all derivative grids (both bases, ns 25/51/101, every h), per-extreme
  comparability, calibration rows and the P1 parity record.
- [generated-identities.json](generated-identities.json): SHA-256 of all 62 generated inputs and
  110 Wouts. Wouts are not distributed; rerun the command to regenerate them.

## Regenerated equilibria and the parity checks they passed

| Identity | Wout SHA-256 | Checks passed | Not established |
| --- | --- | --- | --- |
| reference401, vmecpp ns=401 from `evidence/plasma-design-v2/reference-input-401.json` | `0f2b626c4198f45ffa7a609f4f6f4808f4d8c64562bf4b1a6dd63c22667d5020` | boundary = input to 1e-12; edge flux; starter 64 interior points/target B 2.0e-10/3.4e-10; ideal narrow score vs #35 archive 3.4e-8; headroom fine/geometry checks reproduce to ≤1.5e-15 (#19) | dense identity with the original Wout |
| selected401, vmecpp ns=401 from `evidence/plasma-balanced-v1/selected-input-401.json` | `3d5df9a4df4d5225a29e3640064ef9e14974d141223b1c29bb3dd703c521b9c7` | ier 0; ideal narrow score 1.100894231e-5 vs #35 archive 1.100894343e-5 (1.0e-7) | dense identity; interior/starter samples (none published for this target) |

Scalar parity qualifies these two files only. It does not qualify a new J target, which needs
its own frozen input, Wout identity and checks before evaluation (see the protocol).
