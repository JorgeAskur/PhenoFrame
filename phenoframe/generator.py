"""Python wrapper for the C++ stochastic plant descriptor generator.

This module exposes a :class:`MaizeGenerator` class that drives the C++
``MaizeGenerator`` implementation (from the sibling ``MaizeProceduralModel``
project) via its headless command-line mode. It is the stochastic-sampling
counterpart to :class:`phenoframe.Maize`, which wraps the procedural mesh engine.

Two pieces of functionality are exposed:

1. **Configuration synthesis.** ``derive_config`` reads a base config XML
   (typically the maize template shipped with MaizeProceduralModel), overrides
   the per-parameter means / standard deviations / canopy-offset curves with
   values supplied by the caller, and writes the resulting fitted config.

2. **Batch generation.** ``generate`` invokes the C++ executable in headless
   mode to sample a batch of plant descriptors (and optionally OBJ meshes) into
   a chosen output directory, then returns the generated descriptor paths so
   they can be read back by :func:`phenoframe.compute_traits_from_descriptor`.

Locating the binary
-------------------

The wrapper does not ship its own ``Maize.exe``. It locates an existing build
via the following resolution order:

1. Explicit path passed to :class:`MaizeGenerator`.
2. Environment variable ``PHENOFRAME_GENERATOR_EXE`` (or the legacy
   ``PYMAIZE_GENERATOR_EXE``).
3. A short list of conventional locations relative to the current working
   directory and the user's home, including the sibling
   ``MaizeProceduralModel/x64/Release/Maize.exe`` layout used by the project.
4. ``shutil.which`` on the system ``PATH``.

Set the env var (or pass the explicit path) when none of the conventional
locations apply on your machine.

Typical usage
-------------

::

    from phenoframe import MaizeGenerator, compute_traits_from_descriptor

    gen = MaizeGenerator()  # locates Maize.exe automatically
    fitted = gen.derive_config(
        canopy_params={
            "leafAngleBase":  dict(mean=59.0, std=15.0,
                                   bottom_offset=-6.0, top_offset=+3.0),
            "droopinessBase": dict(mean=-76.6, std=63.0,
                                   bottom_offset=-32.0, top_offset=+30.0),
            "leafLengthBase": dict(mean=0.317, std=0.114),
        },
        base_config_path="MaizeProceduralModel/maize_generator_config.xml",
        output_config_path="MaizeProceduralModel/sorghum_config.xml",
        tiller_azimuth_noise=53.0,
        leaf_jitter_scale=0.60,
    )
    xml_paths = gen.generate(
        config_path=fitted,
        output_dir="MaizeProceduralModel/plant_export",
        count=324,
        seed_start=0,
    )
    traits = [compute_traits_from_descriptor(p) for p in xml_paths]
"""
from __future__ import annotations

import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence


class MaizeGeneratorError(RuntimeError):
    """Raised when the C++ generator binary cannot be located or returns
    non-zero exit status."""


_CONVENTIONAL_LOCATIONS: tuple[str, ...] = (
    # Project-relative (when running from a project repo)
    "Maize.exe",
    "x64/Release/Maize.exe",
    "../MaizeProceduralModel/x64/Release/Maize.exe",
    "../ProceduralModelMaize/x64/Release/Maize.exe",
    "external/ProceduralModelMaize/x64/Release/Maize.exe",
    # User-home layout used by the PhenoFrame project setup
    "Documents/Research/MaizeProceduralModel/x64/Release/Maize.exe",
)


def find_maize_exe(explicit_path: Optional[str | os.PathLike] = None) -> Path:
    """Locate the C++ generator binary.

    Resolution order: explicit argument -> ``PHENOFRAME_GENERATOR_EXE`` (or
    legacy ``PYMAIZE_GENERATOR_EXE``) env var
    -> conventional project layouts -> system ``PATH``. Raises
    :class:`MaizeGeneratorError` with the full search list when nothing is
    found, so the error message is actionable.
    """
    tried: list[str] = []

    if explicit_path is not None:
        p = Path(explicit_path).expanduser().resolve()
        tried.append(str(p))
        if p.is_file():
            return p

    env = os.environ.get("PHENOFRAME_GENERATOR_EXE") or os.environ.get("PYMAIZE_GENERATOR_EXE")
    if env:
        p = Path(env).expanduser().resolve()
        tried.append(f"generator env var={p}")
        if p.is_file():
            return p

    for rel in _CONVENTIONAL_LOCATIONS:
        for root in (Path.cwd(), Path.home()):
            p = (root / rel).resolve()
            tried.append(str(p))
            if p.is_file():
                return p

    on_path = shutil.which("Maize.exe") or shutil.which("Maize")
    if on_path:
        tried.append(f"PATH:{on_path}")
        return Path(on_path)

    raise MaizeGeneratorError(
        "Could not locate the C++ generator binary (Maize.exe).\n"
        "Searched:\n  - " + "\n  - ".join(tried) +
        "\nSet PHENOFRAME_GENERATOR_EXE (or PYMAIZE_GENERATOR_EXE) or pass "
        "exe_path=... when constructing MaizeGenerator."
    )


class MaizeGenerator:
    """Wrapper around the C++ stochastic descriptor generator.

    Parameters
    ----------
    exe_path:
        Optional explicit path to ``Maize.exe``. If omitted, the binary is
        located via the resolution order documented in
        :func:`find_maize_exe`.

    Attributes
    ----------
    exe_path: Path
        The resolved path to the generator binary, captured at construction.
    """

    def __init__(self, exe_path: Optional[str | os.PathLike] = None):
        self.exe_path: Path = find_maize_exe(exe_path)

    # ----------------------------------------------------------------- config
    def derive_config(
        self,
        canopy_params: Mapping[str, Mapping[str, float]],
        base_config_path: str | os.PathLike,
        output_config_path: str | os.PathLike,
        tiller_azimuth_noise: Optional[float] = None,
        leaf_jitter_scale: Optional[float] = None,
        extra_scalars: Optional[Mapping[str, float]] = None,
        root_attrs: Optional[Mapping[str, str]] = None,
        indent: bool = True,
    ) -> Path:
        """Build a fitted generator config by overriding values in a base XML.

        For each entry in ``canopy_params`` the function sets the
        ``mean``/``stdDev`` attributes of the named element and (optionally)
        the bottom/middle/top values of the corresponding ``<X>OffsetCurve``
        element (or ``<X>Curve`` for multiplicative parameters).

        Parameters
        ----------
        canopy_params:
            Mapping ``parameter_name -> spec``. Each spec is itself a mapping
            with the keys:

            - ``mean`` (required), ``std`` (required) — per-plant Gaussian
              parameters.
            - ``bottom_offset``, ``middle_offset``, ``top_offset`` (optional)
              — additive canopy curve offsets written to the named
              ``OffsetCurve`` element.
            - ``bottom_scale``, ``middle_scale``, ``top_scale`` (optional)
              — multiplicative canopy curve factors written to the named
              ``Curve`` element.
            - ``curve_element`` (optional) — explicit override of the curve
              element name; defaults to ``"<param>OffsetCurve"`` or
              ``"<param>Curve"`` based on which keys are present.

        base_config_path:
            Path to the template XML config (e.g.
            ``MaizeProceduralModel/maize_generator_config.xml``).

        output_config_path:
            Where to write the fitted XML. Parent directories are created if
            needed.

        tiller_azimuth_noise:
            Optional override of the ``<tillerAzimuthNoise value="..."/>``
            scalar (uniform amplitude in degrees for per-leaf azimuth jitter).

        leaf_jitter_scale:
            Optional override of the ``<leafJitterScale value="..."/>``
            multiplier on per-leaf trait jitter.

        extra_scalars:
            Mapping ``element_name -> value`` for any additional scalar
            overrides (``value`` attribute). Useful for tweaking visual
            defaults such as ``widthTaper``.

        root_attrs:
            Mapping ``attr_name -> str`` of attributes to set on the root
            element (e.g. ``{"speciesName": "Sorghum"}``). Useful for tagging
            the population the config represents.

        indent:
            If ``True`` (default), pretty-print the output XML using 4-space
            indentation. Set ``False`` to preserve a single-line format that
            byte-matches a hand-edited template.

        Returns
        -------
        Path
            Resolved path to the written config file.
        """
        base = Path(base_config_path).expanduser().resolve()
        out = Path(output_config_path).expanduser().resolve()
        if not base.is_file():
            raise FileNotFoundError(f"base config not found: {base}")

        tree = ET.parse(str(base))
        root = tree.getroot()

        if root_attrs:
            for k, v in root_attrs.items():
                root.set(k, str(v))

        for name, spec in canopy_params.items():
            self._set_mean_std(root, name, float(spec["mean"]), float(spec["std"]))
            curve = self._infer_curve_element(name, spec)
            if curve is not None:
                bot, mid, top = self._curve_triple(spec)
                self._set_curve(root, curve, bot, mid, top)

        if tiller_azimuth_noise is not None:
            self._set_value(root, "tillerAzimuthNoise", float(tiller_azimuth_noise))
        if leaf_jitter_scale is not None:
            self._set_value(root, "leafJitterScale", float(leaf_jitter_scale))

        if extra_scalars:
            for name, value in extra_scalars.items():
                self._set_value(root, name, float(value))

        if indent:
            try:
                ET.indent(tree, space="    ")
            except AttributeError:
                pass

        out.parent.mkdir(parents=True, exist_ok=True)
        tree.write(str(out), encoding="utf-8", xml_declaration=True)
        return out

    @staticmethod
    def _set_mean_std(root: ET.Element, name: str, mean: float, std: float) -> None:
        for el in root.iter(name):
            el.set("mean", f"{mean:.6f}")
            el.set("stdDev", f"{std:.6f}")
            return
        raise KeyError(f"<{name}> element not found in base config")

    @staticmethod
    def _set_curve(root: ET.Element, name: str, bottom: float, middle: float, top: float) -> None:
        for el in root.iter(name):
            el.set("bottom", f"{bottom:.6f}")
            el.set("middle", f"{middle:.6f}")
            el.set("top", f"{top:.6f}")
            return
        raise KeyError(f"<{name}> curve element not found in base config")

    @staticmethod
    def _set_value(root: ET.Element, name: str, value: float) -> None:
        for el in root.iter(name):
            el.set("value", f"{value:.6f}")
            return
        raise KeyError(f"<{name}> scalar element not found in base config")

    @staticmethod
    def _infer_curve_element(name: str, spec: Mapping[str, float]) -> Optional[str]:
        # Naming convention used by the maize_generator_config.xml template:
        #   leafAngleBase    -> leafAngleOffsetCurve   (additive)
        #   leafLengthBase   -> leafLengthScaleCurve   (multiplicative)
        #   droopinessBase   -> droopinessOffsetCurve  (additive)
        # Callers can override via spec["curve_element"].
        if "curve_element" in spec:
            return str(spec["curve_element"])
        has_offset = any(k in spec for k in ("bottom_offset", "middle_offset", "top_offset"))
        has_scale = any(k in spec for k in ("bottom_scale", "middle_scale", "top_scale"))
        stem = name.removesuffix("Base") if name.endswith("Base") else name
        if has_offset:
            return f"{stem}OffsetCurve"
        if has_scale:
            return f"{stem}ScaleCurve"
        return None

    @staticmethod
    def _curve_triple(spec: Mapping[str, float]) -> tuple[float, float, float]:
        if "bottom_offset" in spec or "top_offset" in spec or "middle_offset" in spec:
            return (
                float(spec.get("bottom_offset", 0.0)),
                float(spec.get("middle_offset", 0.0)),
                float(spec.get("top_offset", 0.0)),
            )
        return (
            float(spec.get("bottom_scale", 1.0)),
            float(spec.get("middle_scale", 1.0)),
            float(spec.get("top_scale", 1.0)),
        )

    # -------------------------------------------------------------- generate
    def generate(
        self,
        config_path: str | os.PathLike,
        output_dir: str | os.PathLike,
        count: int,
        seed_start: int = 0,
        output_prefix: str = "plant",
        clean_first: bool = False,
        extra_args: Optional[Sequence[str]] = None,
        capture_output: bool = True,
    ) -> list[Path]:
        """Generate ``count`` synthetic plants via the C++ generator.

        Parameters
        ----------
        config_path:
            Path to the generator config XML.
        output_dir:
            Directory to write descriptor XML / OBJ files into. Created if
            needed; pass ``clean_first=True`` to remove any existing contents
            first (useful when re-running experiments).
        count:
            Number of plants to sample.
        seed_start:
            Starting RNG seed. Plant ``i`` uses ``seed_start + i``.
        output_prefix:
            File-stem prefix for generated files (e.g. ``plant`` yields
            ``plant_0000.xml`` ...).
        clean_first:
            If ``True``, ``rmtree`` the output directory before generating.
        extra_args:
            Additional command-line arguments appended to the headless
            invocation (e.g. flags this wrapper doesn't yet expose).
        capture_output:
            If ``False``, the generator's stdout/stderr are streamed live to
            the terminal instead of being captured. Useful when generating
            very large batches and wanting progress visible.

        Returns
        -------
        list[Path]
            Paths of generated descriptor XMLs, sorted by name.
        """
        config = Path(config_path).expanduser().resolve()
        out = Path(output_dir).expanduser().resolve()
        if not config.is_file():
            raise FileNotFoundError(f"config not found: {config}")
        if clean_first and out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True, exist_ok=True)

        # The C++ generator writes outputs relative to its working directory.
        # We run it from the config's parent so paths in the config XML resolve
        # consistently with the existing project layout.
        cwd = config.parent
        rel_out = os.path.relpath(out, cwd)
        cmd = [
            str(self.exe_path),
            "--headless",
            "--config", config.name,
            "--output", f"{rel_out}/{output_prefix}",
            "--count", str(int(count)),
            "--seed-start", str(int(seed_start)),
        ]
        if extra_args:
            cmd.extend(str(a) for a in extra_args)

        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=capture_output,
            text=True,
        )
        if result.returncode != 0:
            stderr = result.stderr if capture_output else "(stderr was streamed)"
            raise MaizeGeneratorError(
                f"Maize.exe headless failed (exit {result.returncode}).\n"
                f"Command: {' '.join(cmd)}\n"
                f"Working dir: {cwd}\n"
                f"--- stderr ---\n{stderr}"
            )

        return sorted(out.glob(f"{output_prefix}_*.xml"))

    # -------------------------------------------------------------- helpers
    def version(self) -> str:
        """Return whatever the generator binary prints under ``--version`` or
        a string indicating the version flag is unsupported."""
        try:
            result = subprocess.run(
                [str(self.exe_path), "--version"],
                capture_output=True, text=True, timeout=5,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return "unknown (version probe failed)"
        if result.returncode != 0:
            return "unknown (binary does not support --version)"
        return result.stdout.strip() or "unknown (empty version output)"

    def __repr__(self) -> str:
        return f"MaizeGenerator(exe_path={self.exe_path!s})"
