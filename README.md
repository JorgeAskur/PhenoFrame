# PhenoFrame

Toolkit for descriptor-based procedural maize plant modeling and phenotyping trait computation.

## Overview

PhenoFrame provides three main capabilities:

1. **Procedural Plant Generation** — Build 3D maize plant meshes (OBJ) from parametric descriptions via a Python wrapper around a C++ geometry engine. The Python wrapper exposes all geometry-related parameters from the C++ engine (26 per-leaf attributes, including surface noise, midrib geometry, ligule shaping, and wave parameters) and defaults to maximum mesh resolution (200x200 tessellation density).
2. **Trait Computation** — Compute leaf-level phenotyping traits (leaf length, angle, connection point, tip position) from plant descriptor XML files using either the C++ engine or a pure-Python forward model. Both pipelines produce functionally identical results (length MAE = 0.0001 mm, inclination MAE = 0.0015°).
3. **Stochastic Plant Generation** — Sample synthetic descriptor populations from a hierarchical model (4 latent factors + per-node canopy curves + whorl compression) implemented in the sibling C++ `MaizeGenerator` project. A thin Python wrapper (`phenoframe.MaizeGenerator`) handles binary discovery, config-XML derivation from population statistics, and batch invocation of `Maize.exe --headless` — no manual subprocess plumbing required.

The descriptor format is documented in [`DESCRIPTOR_FORMAT.md`](DESCRIPTOR_FORMAT.md), with a formal XSD at [`phenoframe/schemas/descriptor.xsd`](phenoframe/schemas/descriptor.xsd).

An interactive application for the Procedural Model and Stochastic Plant Generator is available at [Maize Procedural Model](https://github.com/JorgeAskur/ProceduralModelMaize).
## Requirements

- Python ≥ 3.10
- A platform with a compiled `maize_c_api` shared library:
  - **Windows**: pre-built DLL is bundled with the package — no build step required.
  - **Linux / macOS**: build from source via `MaizeModel/build_api.sh` (needs `g++` ≥ 7 with C++17 support).

The package targets Python 3.10–3.13 on Linux, Windows, and macOS; the test suite (`pytest`) runs on any platform with a built `maize_c_api` library.

## Installation

The core package — mesh generation (`phenoframe.Maize`), trait extraction (`compute_traits_from_descriptor`), and the stochastic-generator wrapper (`MaizeGenerator`) — has no required runtime Python dependencies. The skeleton-to-descriptor conversion module (`phenoframe.skeleton_to_descriptor`, `phenoframe.skeleton_traits`) additionally requires **NumPy**, **SciPy**, and **pandas**; install them with the `[conversion]` extra (`pip install phenoframe[conversion]`), also included in `[experiments]`. `pytest` is only needed to run the test suite, and `matplotlib` / `jupyter` are only needed for the example notebook (both pulled in via the `[examples]` extra).

### From source

```bash
pip install .
```

### Development install with tests

```bash
pip install -e ".[test]"
python -m pytest          # 64 tests, runs in ~1 s on a developer machine
```

### Notebook examples

```bash
pip install -e ".[examples]"
jupyter notebook examples/quickstart.ipynb
```

### Building the C++ library (Linux / macOS)

```bash
chmod +x MaizeModel/build_api.sh
./MaizeModel/build_api.sh
```

The build script auto-detects the platform and writes `libmaize_c_api.{so,dylib}` directly into `phenoframe/libs/`, where the wrapper looks for it first.

For an explicit manual build:

```bash
g++ -shared -fPIC -O2 -std=c++17 -DMAIZE_C_API_BUILD \
    MaizeModel/maize_c_api.cpp MaizeModel/Maize.cpp MaizeModel/Descriptor.cpp \
    MaizeModel/vect3d.cpp MaizeModel/tinyxml2.cpp \
    -IMaizeModel -o phenoframe/libs/libmaize_c_api.so
```

### Building the C++ library (Windows / MSYS2)

```cmd
MaizeModel\build_api.bat
```

The script expects `g++` from MSYS2 / MinGW on the `PATH`. The pre-built DLL shipped in `phenoframe/libs/maize_c_api.dll` is functionally equivalent.

## Coordinate frame

All world-space coordinates returned by both pipelines (C++ and pure-Python) are **Y-up**:

- **+Y** is the stem growth direction.
- **X / Z** form the horizontal ground plane.
- Lengths are in **meters**; angles are in **degrees**.

For a leaf trait dict, `connection_point` is the leaf-stem junction, `tip_position` is the blade tip, and `leaf_angle` decomposes into:

- `azimuth_deg`: orientation of the base tangent in the XZ plane, in `[0, 360)`.
- `inclination_deg`: angle of the base tangent above the horizontal, in `[-90, 90]`.

## Quick start

### Generate a plant from code

```python
from phenoframe import Maize, Tiller, Leaf

with Maize() as gen:
    gen.reset()
    gen.set_species("maize")
    gen.set_global_settings(
        density_v=120.0, density_u=120.0,
        stem_density_scale=0.15, stem_row_override=0,
        stem_ribbon_spacing=0.15,
    )
    tiller_idx = gen.add_tiller(Tiller().to_ctype())
    gen.add_leaf(tiller_idx, Leaf(
        id=0, distance=0.24, leaf_length=0.45,
        leaf_width=0.15, leaf_angle=35, droopiness=-30.0,
    ).to_ctype())
    gen.rebuild()
    gen.save_obj("plants/my_plant")
```

### Load a plant from XML and export

```python
from phenoframe import from_xml_to_obj

from_xml_to_obj("plants/plant_0.xml", "plants/output")
```

### Compute phenotyping traits (no DLL required)

```python
from pathlib import Path
from phenoframe import compute_traits_from_descriptor

traits = compute_traits_from_descriptor(Path("plants/plant_0.xml"))
for t in traits:
    print(f"Leaf {t['leaf_index']}: length={t['leaf_length']:.4f}")
```

## API reference

### Procedural model — `phenoframe.Maize`

| Method | Description |
|---|---|
| `Maize(library_path=None)` | Load the shared library and create a fresh plant handle. Raises `CApiVersionMismatch` on ABI mismatch. |
| `reset()` | Reset the plant to an empty state. |
| `set_species(name)` | Set the free-form species identifier. |
| `set_global_settings(density_v, density_u, stem_density_scale, stem_row_override, stem_ribbon_spacing)` | Configure mesh tessellation and stem ribbon parameters. |
| `default_tiller() -> TillerDesc` | Get the engine's default tiller parameters as a ctypes struct. |
| `default_leaf() -> LeafDesc` | Get the engine's default leaf parameters as a ctypes struct. |
| `add_tiller(desc=None) -> int` | Append a tiller; returns its index. |
| `set_tiller(index, desc) -> bool` | Replace the tiller at `index`. |
| `tiller_count() -> int` | Number of tillers attached. |
| `add_leaf(tiller_index, desc=None) -> int` | Append a leaf to a tiller; returns its index. |
| `set_leaf(tiller_index, leaf_index, desc) -> bool` | Replace a specific leaf. |
| `leaf_count(tiller_index) -> int` | Number of leaves on a tiller. |
| `rebuild() -> bool` | Regenerate the mesh from the current descriptor (call before reading geometry/traits). |
| `triangle_count() -> int` | Triangle count of the most recently rebuilt mesh. |
| `leaf_spline_trait_count() -> int` | Number of leaves with available spline traits. |
| `get_leaf_spline_traits(leaf_index) -> dict` | World-space trait dict (see *Coordinate frame* above). |
| `save_obj(path, leaf_texture_path, stem_texture_path, include_stem, include_leaves, separate_leaves) -> bool` | Export Wavefront OBJ + MTL. |
| `load_xml(path, rebuild=True) -> bool` | Load a descriptor XML. |
| `save_xml(path) -> bool` | Serialize the current descriptor to XML. |
| `close()` | Release the C++ handle (also runs on `__exit__`). |

### Builders

| Class | Purpose |
|---|---|
| `Leaf(...)` | Pythonic dataclass for a leaf descriptor. Call `.to_ctype()` to get a `LeafDesc` payload. |
| `Tiller(...)` | Pythonic dataclass for a tiller descriptor. Call `.to_ctype()` to get a `TillerDesc` payload. |

### Trait pipeline (pure Python)

| Function | Description |
|---|---|
| `compute_traits_from_descriptor(xml_path) -> list[dict]` | Parse a descriptor XML and return one trait dict per leaf. |
| `write_traits_xml(out_path, source_descriptor, traits)` | Serialize traits to a flat XML for downstream pipelines. |

### Stochastic generator wrapper — `phenoframe.MaizeGenerator`

Wraps the C++ `MaizeGenerator` from the sibling [`MaizeProceduralModel`](https://github.com/) project via its headless CLI. PhenoFrame does not ship the binary; the wrapper discovers an existing build at construction time.

| Method | Description |
|---|---|
| `MaizeGenerator(exe_path=None)` | Locate `Maize.exe`. Resolution order: explicit arg → `PYMAIZE_GENERATOR_EXE` env var → conventional project paths (`./Maize.exe`, `../MaizeProceduralModel/x64/Release/Maize.exe`, etc.) → system `PATH`. Raises `MaizeGeneratorError` with the full search list if nothing is found. |
| `derive_config(canopy_params, base_config_path, output_config_path, ...)` | Build a fitted generator config XML by overriding values in a base template. Handles mean/stdDev pairs for plant-level parameters and bottom/middle/top values for per-rank canopy curves. Curve element names are inferred from the `<param>Base` convention (`leafLengthBase` → `leafLengthScaleCurve`, `leafAngleBase` → `leafAngleOffsetCurve`, etc.) and can be overridden per parameter. |
| `generate(config_path, output_dir, count, seed_start=0, ...)` | Invoke `Maize.exe --headless`; returns the list of generated descriptor XML paths so they can be fed directly to `compute_traits_from_descriptor`. |
| `find_maize_exe(explicit_path=None)` | Module-level path-discovery helper, also exported. |
| `MaizeGeneratorError` | Raised when the binary cannot be located or exits non-zero. |

Minimal example:

```python
from phenoframe import MaizeGenerator, compute_traits_from_descriptor

gen = MaizeGenerator()  # locates Maize.exe automatically
gen.derive_config(
    canopy_params={
        "leafAngleBase":  dict(mean=59.0, std=15.0,
                                bottom_offset=-6.0, top_offset=+3.0),
        "droopinessBase": dict(mean=-76.6, std=63.0,
                                bottom_offset=-32.0, top_offset=+30.0),
        "leafLengthBase": dict(mean=0.317, std=0.114,
                                bottom_scale=0.76, middle_scale=1.0, top_scale=0.94),
    },
    base_config_path="MaizeProceduralModel/maize_generator_config.xml",
    output_config_path="MaizeProceduralModel/sorghum_config.xml",
    tiller_azimuth_noise=53.0,
    leaf_jitter_scale=0.60,
)
xml_paths = gen.generate(
    config_path="MaizeProceduralModel/sorghum_config.xml",
    output_dir="MaizeProceduralModel/plant_export",
    count=324,
    seed_start=0,
)
traits = [compute_traits_from_descriptor(p) for p in xml_paths]
```

### Convenience

| Function | Description |
|---|---|
| `build_example(output_obj="plants/example") -> bool` | Build a small five-leaf demo plant and export OBJ. |
| `from_xml_to_obj(xml_path, output_obj, ...) -> bool` | One-shot: load descriptor, rebuild, export OBJ. |

### Versioning

| Symbol | Description |
|---|---|
| `phenoframe.__version__` | Package version (currently `"0.1.0"`). |
| `phenoframe.wrapper.EXPECTED_C_API_VERSION` | The C-ABI version this Python wrapper expects (currently `(1, 1, 0)`). |
| `phenoframe.wrapper.CApiVersionMismatch` | Raised on major-version mismatch or struct-layout mismatch between the loaded shared library and the Python wrapper. |

## Project structure

```
phenoframe/                   # Importable Python package
  __init__.py              # Public API re-exports
  _version.py              # __version__
  wrapper.py               # ctypes wrapper around the C++ mesh engine
  generator.py             # Wrapper around the C++ stochastic generator (Maize.exe --headless)
  traits.py                # Pure-Python trait computation from descriptor XML
  skeleton_to_descriptor.py # Inverse fitting: skeleton → procedural / spline-override descriptor
  skeleton_traits.py       # Trait extraction directly from voxel skeletons
  libs/                    # Bundled platform binaries
    maize_c_api.dll        # Windows
    libmaize_c_api.so      # Linux (built from source)
    libmaize_c_api.dylib   # macOS (built from source)
  schemas/
    descriptor.xsd         # Formal XML schema for plant descriptors
MaizeModel/                # C++ source for the geometry engine
  maize_c_api.h / .cpp     # C ABI (v1.1.0)
  Maize.cpp / Maize.h      # Core geometry engine
  Descriptor.cpp / .h      # XML descriptor parser
  vect3d.cpp / .h          # Vector math
  tinyxml2.cpp / .h        # Bundled XML parser
  build_api.bat / .sh      # Build scripts (output goes to phenoframe/libs/)
plants/                    # Sample data
  plant_0.xml              # Sample plant descriptor
  maize_leaf.png           # Leaf texture
  maize_stem_texture.png   # Stem texture
tests/                     # Pytest suite (64 tests)
experiments/               # Validation experiments (reproduce the paper)
  experiments/
    README.md              # Experiment overview and status
    EXPERIMENT_RESULTS.md  # Detailed results report
    results/               # Notebooks + scripts (generated outputs are gitignored)
DESCRIPTOR_FORMAT.md       # Human-readable descriptor spec
LICENSE                    # MIT
pyproject.toml             # Packaging metadata
```

## License

[MIT](LICENSE).
