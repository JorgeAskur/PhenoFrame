# PyMaize

Toolkit for procedural maize plant modeling and phenotyping trait computation.

## Overview

PyMaize provides two main capabilities:

1. **Procedural Plant Generation** — Build 3D maize plant meshes (OBJ) from parametric descriptions via a Python wrapper around a C++ geometry engine.
2. **Trait Computation** — Compute leaf-level phenotyping traits (leaf length, angle, connection point, tip position) from plant descriptor XML files using spline reconstruction.

## Project Structure

```
MaizeWrapper.py          # Python ctypes wrapper for the C++ procedural model
Traits/
  trait_lib.py           # Pure-Python trait computation from descriptor XML
  __init__.py
MaizeModel/              # C++ procedural model source and pre-built DLL
  maize_c_api.h          # C API header
  maize_c_api.cpp        # C API implementation
  Maize.cpp / Maize.h    # Core geometry engine
  Descriptor.cpp / .h    # XML descriptor parser
  vect3d.cpp / .h        # Vector math utilities
  tinyxml2.cpp / .h      # XML parsing (bundled dependency)
  maize_c_api.dll        # Pre-built Windows DLL
  build_api.bat / .sh    # Build scripts
plants/
  plant_0.xml            # Sample plant descriptor
  maize_leaf.png         # Leaf texture
  maize_stem_texture.png # Stem texture
test_MaizePM.py          # Demo: build plants from code and XML
test_TraitComputation.py # Demo: compute traits from descriptor XML
```

## Requirements

- Python 3.10+
- Windows (pre-built DLL included) or Linux/macOS (compile from source)

## Quick Start

### Generate a plant from code

```python
from MaizeWrapper import Maize, Tiller, Leaf

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
from MaizeWrapper import from_xml_to_obj

from_xml_to_obj("plants/plant_0.xml", "plants/output")
```

### Compute phenotyping traits

```python
from Traits import compute_traits_from_descriptor
from pathlib import Path

traits = compute_traits_from_descriptor(Path("plants/plant_0.xml"))
for t in traits:
    print(f"Leaf {t['leaf_index']}: length={t['leaf_length']:.4f}")
```

Or from the command line:

```bash
python test_TraitComputation.py --xml plants/plant_0.xml --out-xml traits.xml
```

## Building the C++ Library from Source

### Windows (MSYS2 / MinGW)

```bash
cd MaizeModel
./build_api.bat
```

### Linux / macOS

```bash
cd MaizeModel
chmod +x build_api.sh
./build_api.sh
```

Or manually:

```bash
g++ -shared -O2 -std=c++17 -DMAIZE_C_API_BUILD \
    maize_c_api.cpp Maize.cpp Descriptor.cpp vect3d.cpp tinyxml2.cpp \
    -I. -o libmaize_c_api.so
```
