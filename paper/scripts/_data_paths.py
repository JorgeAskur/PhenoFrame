"""Locations for the tracked paper inputs and locally generated outputs."""

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PAPER = REPO / "paper"
SORGHUM = PAPER / "source_data" / "sorghum"
MAIZE = PAPER / "source_data" / "maize"
GENERATION = PAPER / "source_data" / "generation"
SUMMARY = PAPER / "source_data" / "summary"
RESULTS = REPO / "experiments" / "experiments" / "results"
GENERATED = PAPER / "generated"
