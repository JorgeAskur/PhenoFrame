# PyMaize

Toolkit for procedural maize plant modeling and phenotyping trait computation.

## Overview

PyMaize provides two main capabilities:

1. **Procedural Plant Generation** — Build 3D maize plant meshes (OBJ) from parametric descriptions via a Python wrapper around a C++ geometry engine.
2. **Trait Computation** — Compute leaf-level phenotyping traits (leaf length, angle, connection point, tip position) from plant descriptor XML files using spline reconstruction.

The descriptor format is documented in [`DESCRIPTOR_FORMAT.md`](DESCRIPTOR_FORMAT.md), with a formal XSD at [`pymaize/schemas/descriptor.xsd`](pymaize/schemas/descriptor.xsd).

## Requirements

- Python ≥ 3.10
- A platform with a compiled `maize_c_api` shared library:
  - **Windows**: pre-built DLL is bundled with the package — no build step required.
  - **Linux / macOS**: build from source via `MaizeModel/build_api.sh` (needs `g++` ≥ 7 with C++17 support).

The CI workflow ([`.github/workflows/test.yml`](.github/workflows/test.yml)) exercises Ubuntu, Windows, and macOS on Python 3.10 and 3.13.

## Installation

PyMaize has zero runtime Python dependencies. `pytest` is only needed to run the test suite, and `matplotlib` / `jupyter` are only needed for the example notebook (both pulled in via the `[examples]` extra).

### From source

```bash
pip install .
```

### Development install with tests

```bash
pip install -e ".[test]"
python -m pytest          # 40 tests, runs in < 1 s on a developer machine
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

The build script auto-detects the platform and writes `libmaize_c_api.{so,dylib}` directly into `pymaize/libs/`, where the wrapper looks for it first.

For an explicit manual build:

```bash
g++ -shared -fPIC -O2 -std=c++17 -DMAIZE_C_API_BUILD \
    MaizeModel/maize_c_api.cpp MaizeModel/Maize.cpp MaizeModel/Descriptor.cpp \
    MaizeModel/vect3d.cpp MaizeModel/tinyxml2.cpp \
    -IMaizeModel -o pymaize/libs/libmaize_c_api.so
```

### Building the C++ library (Windows / MSYS2)

```cmd
MaizeModel\build_api.bat
```

The script expects `g++` from MSYS2 / MinGW on the `PATH`. The pre-built DLL shipped in `pymaize/libs/maize_c_api.dll` is functionally equivalent.

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
from pymaize import Maize, Tiller, Leaf

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
from pymaize import from_xml_to_obj

from_xml_to_obj("plants/plant_0.xml", "plants/output")
```

### Compute phenotyping traits (no DLL required)

```python
from pathlib import Path
from pymaize import compute_traits_from_descriptor

traits = compute_traits_from_descriptor(Path("plants/plant_0.xml"))
for t in traits:
    print(f"Leaf {t['leaf_index']}: length={t['leaf_length']:.4f}")
```

Or from the command line:

```bash
python test_TraitComputation.py --xml plants/plant_0.xml --out-xml traits.xml
```

## API reference

### Procedural model — `pymaize.Maize`

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

### Convenience

| Function | Description |
|---|---|
| `build_example(output_obj="plants/example") -> bool` | Build a small five-leaf demo plant and export OBJ. |
| `from_xml_to_obj(xml_path, output_obj, ...) -> bool` | One-shot: load descriptor, rebuild, export OBJ. |

### Versioning

| Symbol | Description |
|---|---|
| `pymaize.__version__` | Package version (currently `"0.1.0"`). |
| `pymaize.wrapper.EXPECTED_C_API_VERSION` | The C-ABI version this Python wrapper expects (currently `(1, 0, 0)`). |
| `pymaize.wrapper.CApiVersionMismatch` | Raised on major-version mismatch or struct-layout mismatch between the loaded shared library and the Python wrapper. |

## Project structure

```
pymaize/                   # Importable Python package
  __init__.py              # Public API re-exports
  _version.py              # __version__
  wrapper.py               # ctypes wrapper around the C++ engine
  traits.py                # Pure-Python trait computation from descriptor XML
  libs/                    # Bundled platform binaries
    maize_c_api.dll        # Windows
    libmaize_c_api.so      # Linux (built from source)
    libmaize_c_api.dylib   # macOS (built from source)
  schemas/
    descriptor.xsd         # Formal XML schema for plant descriptors
MaizeModel/                # C++ source for the geometry engine
  maize_c_api.h / .cpp     # C ABI
  Maize.cpp / Maize.h      # Core geometry engine
  Descriptor.cpp / .h      # XML descriptor parser
  vect3d.cpp / .h          # Vector math
  tinyxml2.cpp / .h        # Bundled XML parser
  build_api.bat / .sh      # Build scripts (output goes to pymaize/libs/)
plants/                    # Sample data
  plant_0.xml              # Sample plant descriptor
  maize_leaf.png           # Leaf texture
  maize_stem_texture.png   # Stem texture
tests/                     # Pytest suite (40 tests)
examples/                  # Notebooks
DESCRIPTOR_FORMAT.md       # Human-readable descriptor spec
LICENSE                    # MIT
pyproject.toml             # Packaging metadata
test_MaizePM.py            # Demo: build plants from code and from XML
test_TraitComputation.py   # Demo: compute traits from descriptor XML
```

## License

[MIT](LICENSE).
