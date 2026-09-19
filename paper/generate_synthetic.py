"""Regenerate the 324 descriptor instances used in the generation experiment."""

from __future__ import annotations

import argparse
from pathlib import Path

from phenoframe import MaizeGenerator


PAPER = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=324)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--exe", type=Path, help="Path to the ProceduralModelMaize executable")
    parser.add_argument("--config", type=Path, default=PAPER / "configs" / "sorghum_generator_config.xml")
    parser.add_argument("--output-dir", type=Path, default=PAPER / "generated" / "sorghum")
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")

    generator = MaizeGenerator(exe_path=args.exe)
    paths = generator.generate(
        config_path=args.config,
        output_dir=args.output_dir,
        count=args.count,
        seed_start=args.seed_start,
        clean_first=False,
    )
    if len(paths) != args.count:
        raise RuntimeError(f"Expected {args.count} descriptor XMLs, found {len(paths)}")
    print(f"Generated {len(paths)} descriptors in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
