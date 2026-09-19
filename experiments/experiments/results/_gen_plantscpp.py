"""Regenerate the PlantsC++ dataset (500 deterministic C++ plants, seeds 0-499)
used by notebooks 01 and 02. Same maize config the original 500 were sampled
from, so this is a faithful reproduction of the §2.1 validation population.
"""
import sys
from pathlib import Path

from phenoframe.generator import MaizeGenerator

REPO = Path(__file__).resolve().parents[3]
CFG = REPO / "paper" / "configs" / "maize_generator_config.xml"
OUT = REPO / "experiments" / "PlantsC++"
COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 500

gen = MaizeGenerator()  # uses PHENOFRAME_GENERATOR_EXE when set
print(f"exe: {gen.exe_path}", flush=True)
print(f"config: {CFG}", flush=True)
print(f"output: {OUT}", flush=True)
print(f"generating {COUNT} plants (seeds 0..{COUNT-1}) ...", flush=True)

paths = gen.generate(
    config_path=CFG,
    output_dir=OUT,
    count=COUNT,
    seed_start=0,
    output_prefix="plant",
    clean_first=True,
    capture_output=True,
)
print(f"DONE: wrote {len(paths)} descriptor XMLs to {OUT}", flush=True)
