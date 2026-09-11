"""Sequentially execute the experiment notebooks, logging status for each.

Runs each notebook in-process via nbclient with the results dir as cwd,
PYTHONPATH already pointing at the experiments/ root. Continues on failure
so we get a full pass/fail report in one go.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient.client import NotebookClient

RESULTS = Path(__file__).resolve().parent
OUTDIR = RESULTS / "_executed"
OUTDIR.mkdir(exist_ok=True)

# Run the kernel with cwd = the `01_phyllotaxis` junction (-> results), so the
# notebooks' `Path.cwd().name == '01_phyllotaxis'` path logic resolves exactly
# as in the original repo layout (REPO_ROOT, OUT, FIG all line up).
JUNCTION = RESULTS.parent / "01_phyllotaxis"
RUN_CWD = str(JUNCTION if JUNCTION.exists() else RESULTS)

ORDER = [
    "01_sec2.1_python_cpp_validation.ipynb",
    "02_sec2.1_descriptor_trait_validation.ipynb",
    "03_sec2.2_phenosuite_traits_vs_gold.ipynb",
    "04_sec2.3_geometric_fidelity.ipynb",
    "05_sec2.3_procedural_round_trip.ipynb",
    "06_sec2.4_heritability.ipynb",
    "07_sec2.5_distribution_match.ipynb",
    "08_sec2.5_data_driven_generation.ipynb",
    "09_sec2.5_holdout_validation.ipynb",
    "10_phase1_gwas_marker_replication.ipynb",
    "99_appendix_e1_trait_concordance.ipynb",
]

# Allow starting partway through: pass notebook filenames as args.
if len(sys.argv) > 1:
    ORDER = sys.argv[1:]

results = []
for name in ORDER:
    src = RESULTS / name
    if not src.exists():
        print(f"SKIP  {name}: file not found", flush=True)
        results.append((name, "MISSING", 0.0, ""))
        continue
    print(f"\n=== RUN  {name} ===", flush=True)
    t0 = time.time()
    nb = nbformat.read(src, as_version=4)
    client = NotebookClient(
        nb,
        timeout=7200,
        kernel_name="python3",
        resources={"metadata": {"path": RUN_CWD}},
        allow_errors=False,
    )
    status, err = "OK", ""
    try:
        client.execute()
    except Exception as e:  # noqa: BLE001
        status = "FAIL"
        err = f"{type(e).__name__}: {str(e)[:400]}"
    finally:
        nbformat.write(nb, OUTDIR / name)
    dt = time.time() - t0
    print(f"=== {status}  {name}  ({dt:.1f}s) ===", flush=True)
    if err:
        print(err, flush=True)
    results.append((name, status, dt, err))

print("\n\n========== SUMMARY ==========", flush=True)
for name, status, dt, err in results:
    print(f"{status:6}  {dt:7.1f}s  {name}", flush=True)
    if err:
        print(f"         {err}", flush=True)
n_ok = sum(1 for r in results if r[1] == "OK")
print(f"\n{n_ok}/{len(results)} notebooks OK", flush=True)
