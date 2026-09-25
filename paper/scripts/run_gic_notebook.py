"""Run experiment notebooks on the GIC 20-genotype maize voxel/skeleton dataset.

Unlike run_maize_notebook.py this does NOT patch clean_segment: the GIC
optim_skeleton.txt files are clean optimized curve skeletons (default
segment_skeleton already matches gold leaf counts), so the maize-reconstruction
medial-axis fix is unnecessary and would only hurt.

Isolation: the kernel runs with cwd = WORKDIR (a dedicated folder), so each
notebook's `outputs/` and `figures/` land there instead of clobbering the
sorghum results in experiments/results/outputs. PYTHONPATH makes the
`experiments` package importable regardless of cwd.
"""
import argparse
import os
import sys
import time
from pathlib import Path
import nbformat
from nbclient.client import NotebookClient

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "experiments" / "experiments" / "results"
REPO_EXP = REPO / "experiments"
WORKDIR = REPO / "paper" / "generated" / "maize_run"
WORKDIR.mkdir(exist_ok=True)
(WORKDIR / "outputs").mkdir(exist_ok=True)
(WORKDIR / "figures").mkdir(exist_ok=True)
EXEC = WORKDIR / "_executed"; EXEC.mkdir(exist_ok=True)

parser = argparse.ArgumentParser(description="Run tracked experiment notebooks on GIC maize skeletons")
parser.add_argument("--voxel-root", type=Path, required=True,
                    help="Directory containing reconstructed/ and skeletons/")
parser.add_argument("notebooks", nargs="+", help="Notebook filenames under experiments/experiments/results")
args = parser.parse_args()

os.environ["PHENOFRAME_VOXEL_PATH"] = str(args.voxel_root.resolve())
os.environ["PHENOFRAME_JENSINA_PATH"] = str(REPO / "paper" / "reference" / "phyllotaxy")
os.environ["PHENOFRAME_EXPERIMENT_ROOT"] = str(REPO_EXP)
os.environ["PHENOFRAME_OUTPUT_DIR"] = str(WORKDIR / "outputs")
os.environ["PHENOFRAME_FIGURE_DIR"] = str(WORKDIR / "figures")
os.environ["PYTHONPATH"] = os.pathsep.join((str(REPO), str(REPO_EXP), os.environ.get("PYTHONPATH", "")))

for name in args.notebooks:
    src = RESULTS / name
    nb = nbformat.read(src, as_version=4)
    print(f"=== RUN {name} on GIC maize (default segmentation) ===", flush=True)
    t0 = time.time()
    client = NotebookClient(
        nb, timeout=21600, kernel_name="python3",
        resources={"metadata": {"path": str(WORKDIR)}},
        allow_errors=False,
    )
    status = "OK"
    try:
        client.execute()
    except Exception as e:  # noqa: BLE001
        status = f"FAIL {type(e).__name__}: {str(e)[:300]}"
    finally:
        nbformat.write(nb, EXEC / name)
    print(f"=== {status}  {name}  ({time.time()-t0:.0f}s) ===", flush=True)
