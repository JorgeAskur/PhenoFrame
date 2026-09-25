"""Create the compact GIC maize Zenodo deposit archive and manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import re
import zipfile
from pathlib import Path


PAPER = Path(__file__).resolve().parent
INCLUDE = ("voxels.txt", "optim_skeleton.txt", "angles.txt", "error.txt")
SCAN_RE = re.compile(
    r"_\d+-(?P<plant_number>\d+)-(?P<genotype>[A-Za-z0-9]+)-"
    r"(?P<replicate>\d+)_(?P<timestamp>\d{4}-\d{2}-\d{2}_"
    r"\d{2}-\d{2}-\d{2}\.\d+)_(?P<reconstruction_id>\d+)$"
)


def declared_count(path: Path) -> int | None:
    if not path.is_file() or path.stat().st_size <= 2:
        return None
    try:
        with path.open("r", encoding="utf-8") as stream:
            return int(stream.readline().strip())
    except (OSError, UnicodeError, ValueError):
        return None


def angle_count(path: Path) -> int | None:
    if not path.is_file() or path.stat().st_size <= 2:
        return None
    try:
        values = path.read_text(encoding="utf-8").split()
        return len(values) // 3 if len(values) % 3 == 0 else None
    except (OSError, UnicodeError):
        return None


def build_manifest(source: Path) -> tuple[list[dict[str, object]], list[Path]]:
    rows: list[dict[str, object]] = []
    files: list[Path] = []
    for scan in sorted(path for path in source.iterdir() if path.is_dir()):
        match = SCAN_RE.search(scan.name)
        ids = match.groupdict() if match else {}
        present = {name: (scan / name).is_file() for name in INCLUDE}
        files.extend(scan / name for name in INCLUDE if present[name])
        error = ""
        if present["error.txt"]:
            error = (scan / "error.txt").read_text(encoding="utf-8", errors="replace").strip()
        rows.append({
            "scan_id": scan.name,
            "plant_number": ids.get("plant_number", ""),
            "genotype": ids.get("genotype", ""),
            "replicate": ids.get("replicate", ""),
            "timestamp": ids.get("timestamp", "").replace("_", "T", 1),
            "reconstruction_id": ids.get("reconstruction_id", ""),
            "has_voxels": present["voxels.txt"],
            "has_optimized_skeleton": present["optim_skeleton.txt"],
            "has_angles": present["angles.txt"],
            "n_voxels": declared_count(scan / "voxels.txt") or "",
            "n_skeleton_points": declared_count(scan / "optim_skeleton.txt") or "",
            "n_leaves": angle_count(scan / "angles.txt") or "",
            "voxels_bytes": (scan / "voxels.txt").stat().st_size if present["voxels.txt"] else 0,
            "skeleton_bytes": (scan / "optim_skeleton.txt").stat().st_size if present["optim_skeleton.txt"] else 0,
            "angles_bytes": (scan / "angles.txt").stat().st_size if present["angles.txt"] else 0,
            "status_or_error": error,
        })
    return rows, files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=PAPER / "generated" / "deposit" / "gic_maize_voxels_skeletons_2024.zip",
    )
    args = parser.parse_args()
    source = args.source.resolve()
    if not source.is_dir():
        parser.error(f"Source directory does not exist: {source}")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    rows, files = build_manifest(source)
    if len(rows) != 4171:
        raise RuntimeError(f"Expected 4,171 scan directories, found {len(rows):,}")
    fields = list(rows[0])
    manifest = io.StringIO(newline="")
    writer = csv.DictWriter(manifest, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

    prefix = "gic_maize_voxels_skeletons_2024"
    readme = (PAPER / "GIC_DATASET_README.md").read_text(encoding="utf-8")
    with zipfile.ZipFile(
        output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True
    ) as archive:
        archive.writestr(f"{prefix}/README.md", readme)
        archive.writestr(f"{prefix}/manifest.csv", manifest.getvalue())
        for index, path in enumerate(files, 1):
            relative = path.relative_to(source)
            archive.write(path, f"{prefix}/scans/{relative.as_posix()}")
            if index % 1000 == 0:
                print(f"Archived {index:,}/{len(files):,} files", flush=True)

    digest = hashlib.sha256()
    with output.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    checksum = output.with_suffix(output.suffix + ".sha256")
    checksum.write_text(f"{digest.hexdigest()}  {output.name}\n", encoding="ascii")
    success = sum(bool(row["has_optimized_skeleton"] and row["has_angles"]) for row in rows)
    print(f"Created {output} ({output.stat().st_size / 1024**3:.2f} GiB)")
    print(f"Scans: {len(rows):,}; optimized skeleton + angles: {success:,}")
    print(f"SHA-256: {digest.hexdigest()}")


if __name__ == "__main__":
    main()
