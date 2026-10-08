"""N14: parity, comparability, derivative reliability and weight calibration (PREREGISTRATION.md).

Usage (repository root, one thread): python feasibility.py SEED_SNAPSHOT OUTDIR
"""

import copy
import json
import sys
import tempfile
import time
from pathlib import Path

import netCDF4
import numpy as np
import vmecpp

sys.path.insert(0, "src")
from fusion_baselines import coil_bounce as bounce  # noqa: E402
from fusion_baselines import coil_fit as fit  # noqa: E402
from fusion_baselines.vmec_trace import trace_geometry  # noqa: E402

MODES = (("rbc", 1, 1), ("rbc", 2, 0), ("zbs", 1, 1), ("zbs", 2, 0))
REF = json.loads(Path("evidence/plasma-design-v2/reference-input-401.json").read_text())
SEL = json.loads(Path("evidence/plasma-balanced-v1/selected-input-401.json").read_text())
ARCHIVED_SEL_NARROW = 1.100894343118703e-05
NS = {25: ([12, 25], [1e-8, 1e-11]), 51: ([12, 25, 51], [1e-8, 1e-10, 1e-11]),
      101: ([12, 25, 51, 101], [1e-8, 1e-10, 1e-11, 1e-11])}
seed_path, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=False)
work = Path(tempfile.mkdtemp(dir=out))
log, t0 = [], time.monotonic()


def mode(data, key, m, n):
    return next(r["value"] for r in data[key] if (r["m"], r["n"]) == (m, n))


def boundary(base, shift):
    data = copy.deepcopy(base)
    for (key, m, n), delta in zip(MODES, shift, strict=True):
        for row in data[key]:
            if (row["m"], row["n"]) == (m, n):
                row["value"] += float(delta)
    return data


def solve(data, ns, tag):
    path = work/f"{tag}-input.json"
    path.write_text(json.dumps(data))
    inp = vmecpp.VmecInput.from_file(path)
    if ns != 401:
        inp.ns_array = np.array(NS[ns][0])
        inp.ftol_array = np.array(NS[ns][1])
        inp.niter_array = np.array([10000]*len(NS[ns][0]))
    started = time.monotonic()
    wout = work/f"{tag}-ns{ns}.nc"
    result = vmecpp.run(inp, max_threads=1, verbose=False)
    result.wout.save(str(wout))
    seconds = time.monotonic() - started
    with netCDF4.Dataset(wout) as d:
        info = {k: float(d[k][...]) for k in ("volume_p", "aspect", "Rmajor_p", "b0")}
        info["ier"] = int(d["ier_flag"][...])
    return wout, info, seconds


def narrow(wout):
    cells = []
    for s in (0.25, 0.5, 0.75):
        m = bounce.measure(trace_geometry(str(wout), s, 801, 16, 2))
        cells += [c for c in m["cells"] if c["q"] in (0.1, 0.3, 0.5, 0.7, 0.9)]
    if len(cells) != 15:
        raise ValueError("incomplete narrow domain")
    return float(np.mean([c["score"] for c in cells]))


def A(data, ns, tag):
    wout, info, seconds = solve(data, ns, tag)
    value = narrow(wout)
    log.append(dict(tag=tag, ns=ns, A=value, solve_s=round(seconds, 2), **info))
    return value, info


report = {}
# P1: selected401 parity at ns=401.
wout, info, seconds = solve(SEL, 401, "selected401")
sel_narrow = narrow(wout)
report["P1"] = dict(narrow=sel_narrow, archived=ARCHIVED_SEL_NARROW,
                    rel=abs(sel_narrow/ARCHIVED_SEL_NARROW - 1), ier=info["ier"], solve_s=seconds,
                    qualified=abs(sel_narrow/ARCHIVED_SEL_NARROW - 1) <= 1e-6 and info["ier"] == 0)
report["P1"]["wout"] = str(wout)

# P2: comparability at ns=25.
_, ref_info = A(REF, 25, "ref")
rows = []
for tag, data in [("selected", SEL)] + [
        (f"{k}{m}{n}{sign:+d}", boundary(REF, [5e-4*sign*(i == j) for j in range(4)]))
        for i, (k, m, n) in enumerate(MODES) for sign in (1, -1)]:
    _, info = A(data, 25, tag)
    rows.append(dict(tag=tag, **{k: info[k]/ref_info[k] - 1 for k in ("volume_p", "aspect",
                                                                        "Rmajor_p", "b0")}))
report["P2"] = dict(reference=ref_info, relative_changes=rows,
                    within_1pct=all(abs(v) <= 0.01 for r in rows for k, v in r.items() if k != "tag"))

# P3: derivative reliability at reference and selected.
plan = [(25, h) for h in (3e-6, 1e-5, 3e-5)] + [(51, h) for h in (1e-5, 3e-5)] + [(101, 1e-5)]
derivs = {}
for base_tag, base in (("ref", REF), ("sel", SEL)):
    for ns, h in plan:
        g = []
        for i in range(4):
            step = [h*(i == j) for j in range(4)]
            plus = A(boundary(base, step), ns, f"{base_tag}-d{i}+h{h}")[0]
            minus = A(boundary(base, [-x for x in step]), ns, f"{base_tag}-d{i}-h{h}")[0]
            g.append((plus - minus)/(2*h))
        derivs[f"{base_tag}|ns{ns}|h{h}"] = g


def check(base_tag, ns_search):
    ref101 = np.array(derivs[f"{base_tag}|ns101|h1e-05"])
    big = abs(ref101) > 0.05*abs(ref101).max()
    search = np.array(derivs[f"{base_tag}|ns{ns_search}|h1e-05"])
    resolution = bool(np.all(np.sign(search[big]) == np.sign(ref101[big]))
                      and np.all(abs(search[big] - ref101[big]) <= 0.2*abs(ref101[big])))
    hs = (3e-6, 1e-5, 3e-5) if ns_search == 25 else (1e-5, 3e-5)
    stack = np.array([derivs[f"{base_tag}|ns{ns_search}|h{h}"] for h in hs])
    spread = (stack.max(axis=0) - stack.min(axis=0))/np.maximum(abs(stack).max(axis=0), 1e-300)
    return dict(resolution_ok=resolution, h_spread=spread.tolist(),
                h_ok=bool(np.all(spread[big] <= 0.10)))


report["P3"] = dict(derivatives=derivs, checks={f"{b}|ns{n}": check(b, n)
                                                 for b in ("ref", "sel") for n in (25, 51)})

# P4: weight calibration along the Step 3 direction with fixed headroom coils.
direction = [mode(SEL, k, m, n) - mode(REF, k, m, n) for k, m, n in MODES]
seed = json.loads(seed_path.read_text())
record = fit.Recorder(out/"recorder", time.monotonic() + 3600)
calibration = []
for t in (0.0, 0.5, 1.0):
    data = boundary(REF, [t*d for d in direction])
    a_value = A(data, 25, f"step3-t{t}")[0]
    model = fit.Model(seed, data, record)
    _, _, metrics = model.evaluate(model.x0)
    calibration.append(dict(t=t, A=a_value, J_flux=metrics["flux_objective"],
                            geometry_penalty=metrics["geometry_penalty"],
                            normal_rms=metrics["normal_rms"]))
dA = calibration[0]["A"] - calibration[2]["A"]
dJ = calibration[2]["J_flux"] - calibration[0]["J_flux"]
weight = float(f"{dJ/dA:.0e}") if dA > 0 else None
report["P4"] = dict(direction=direction, calibration=calibration, dA=dA, dJ_flux=dJ,
                    weight=weight, stop=dA <= 0)
report.update(log=log, seconds=round(time.monotonic() - t0, 1))
(out/"result.json").write_text(json.dumps(report, indent=1))
print(json.dumps({k: v for k, v in report.items() if k not in ("log",)}, indent=1)[:6000])
