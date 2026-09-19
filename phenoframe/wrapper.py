"""Python wrapper for the C++ Maize procedural plant generator.

This module exposes a `Maize` class that mirrors the C ABI declared in
``MaizeModel/maize_c_api.h``. It loads the platform's pre-built shared library
from ``phenoframe/libs/``, validates the ABI version on construction, and offers
a Pythonic interface for assembling tillers and leaves, rebuilding the mesh,
and exporting OBJ / descriptor XML.

Coordinate frame for all geometry returned by this module:

- World-space, **Y-up**.
- ``+Y`` is the stem growth direction.
- ``X`` and ``Z`` form the horizontal ground plane.
- Lengths are in meters; angles are in degrees.
"""

from __future__ import annotations

import ctypes
import os
import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path
import platform
from typing import Optional, Sequence


# C-API version this Python wrapper was built against. Must stay in sync with
# MAIZE_C_API_VERSION_* in MaizeModel/maize_c_api.h. A major mismatch raises;
# a minor mismatch warns; older libraries without version exports also warn.
EXPECTED_C_API_VERSION = (1, 1, 0)


class CApiVersionMismatch(RuntimeError):
    """Raised when the loaded shared library's ABI version is incompatible."""


class _Handle(ctypes.Structure):
    pass


class LeafDesc(ctypes.Structure):
    _fields_ = [
        ("id", ctypes.c_int),
        ("distance", ctypes.c_float),
        ("leaf_length", ctypes.c_float),
        ("leaf_width", ctypes.c_float),
        ("azimuth_deg", ctypes.c_float),
        ("leaf_angle", ctypes.c_float),
        ("droopiness", ctypes.c_float),
        ("stem_inclination_deg", ctypes.c_float),
        ("spline_points", ctypes.c_int),
        ("width_taper", ctypes.c_float),
        ("leaf_twist", ctypes.c_float),
        ("leaf_curl", ctypes.c_float),
        ("wave_l_amp", ctypes.c_float),
        ("wave_l_freq", ctypes.c_float),
        ("wave_l_phase", ctypes.c_float),
        ("wave_r_amp", ctypes.c_float),
        ("wave_r_freq", ctypes.c_float),
        ("wave_r_phase", ctypes.c_float),
        ("surface_noise_amp", ctypes.c_float),
        ("surface_noise_freq", ctypes.c_float),
        ("midrib_tip_taper_start", ctypes.c_float),
        ("midrib_texture_strength", ctypes.c_float),
        ("midrib_width", ctypes.c_float),
        ("ligule_wrap_length_scale", ctypes.c_float),
        ("ligule_unfold_sharpness", ctypes.c_float),
        ("sheath_outer_scale", ctypes.c_float),
        ("use_ctrl_overrides", ctypes.c_int),
    ]


class TillerDesc(ctypes.Structure):
    _fields_ = [
        ("type", ctypes.c_char_p),
        ("radius", ctypes.c_float),
        ("alpha_deg", ctypes.c_float),
        ("stem_shrink", ctypes.c_float),
        ("azimuth_deg", ctypes.c_float),
        ("azimuth_noise", ctypes.c_float),
        ("random_seed", ctypes.c_uint),
    ]


class LeafSplineTraits(ctypes.Structure):
    _fields_ = [
        ("connection_x", ctypes.c_float),
        ("connection_y", ctypes.c_float),
        ("connection_z", ctypes.c_float),
        ("tip_x", ctypes.c_float),
        ("tip_y", ctypes.c_float),
        ("tip_z", ctypes.c_float),
        ("base_tangent_x", ctypes.c_float),
        ("base_tangent_y", ctypes.c_float),
        ("base_tangent_z", ctypes.c_float),
        ("leaf_length", ctypes.c_float),
        ("azimuth_deg", ctypes.c_float),
        ("inclination_deg", ctypes.c_float),
    ]


class Maize:
    """High-level wrapper used to build and export maize plant meshes.

    The wrapper owns a single C++ ``MaizeHandle`` for its lifetime. Use it as
    a context manager so the handle is always released::

        with Maize() as gen:
            gen.reset()
            gen.set_species("maize")
            ...
            gen.save_obj("plants/output")

    On construction the loaded shared library's ABI version and struct sizes
    are validated against ``EXPECTED_C_API_VERSION``. A major-version mismatch
    raises :class:`CApiVersionMismatch`; older libraries without version
    exports emit a :class:`UserWarning`.

    Parameters
    ----------
    library_path:
        Optional explicit path to the shared library. If omitted, the wrapper
        searches ``phenoframe/libs/`` first, then ``MaizeModel/`` build outputs.
    """

    def __init__(self, library_path: Optional[str] = None):
        self._dll_dir_handles = []
        self._lib = self._load_library(library_path)
        self._bind_api()
        self._check_api_compatibility()
        self._handle = self._lib.maize_create()
        if not self._handle:
            raise RuntimeError("maize_create() returned null handle")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        self.close()

    def close(self):
        """Release the underlying C++ handle. Safe to call multiple times."""
        if getattr(self, "_handle", None):
            self._lib.maize_destroy(self._handle)
            self._handle = None

    @staticmethod
    def _library_name():
        if platform.system() == "Windows":
            return "maize_c_api.dll"
        if platform.system() == "Darwin":
            return "libmaize_c_api.dylib"
        return "libmaize_c_api.so"

    @staticmethod
    def _candidate_lib_paths(lib_name: str) -> Sequence[Path]:
        package_dir = Path(__file__).resolve().parent
        repo_root = package_dir.parent
        return (
            package_dir / "libs" / lib_name,
            repo_root / "MaizeModel" / lib_name,
            repo_root / "MaizeModel" / "build" / lib_name,
            repo_root / lib_name,
            repo_root / "build" / lib_name,
        )

    def _load_library(self, library_path: Optional[str]) -> ctypes.CDLL:
        if library_path:
            candidates: Sequence[Path] = (Path(library_path),)
        else:
            candidates = self._candidate_lib_paths(self._library_name())

        load_errors = []
        for candidate in candidates:
            if candidate.exists():
                try:
                    self._add_windows_dll_dirs(candidate)
                    return ctypes.cdll.LoadLibrary(str(candidate))
                except OSError as exc:
                    load_errors.append(f"{candidate}: {exc}")
                    continue
        raise FileNotFoundError(
            "Could not find compiled maize_c_api shared library.\n"
            f"Searched: {', '.join(map(str, candidates))}\n"
            + (f"Load errors: {' | '.join(load_errors)}\n" if load_errors else "")
            + "Build it from source via MaizeModel/build_api.sh (Linux/macOS) "
            "or MaizeModel/build_api.bat (Windows)."
        )

    # The DLL search path is a property of the *process*, not of an instance, so
    # directories are registered once and shared. Previously every Maize() added
    # its own cookies and close() never released them; a few hundred instances --
    # ordinary batch work over a plant corpus -- exhausted the process-wide slots,
    # and the next library to call os.add_dll_directory (numpy, pandas, sklearn)
    # failed with WinError 206. Keyed by resolved path so a non-default
    # library_path still registers its own directory exactly once.
    _dll_dirs_added: dict = {}

    def _add_windows_dll_dirs(self, candidate: Path) -> None:
        if platform.system() != "Windows" or not hasattr(os, "add_dll_directory"):
            return

        package_dir = Path(__file__).resolve().parent
        repo_root = package_dir.parent

        dirs = [
            candidate.parent,
            package_dir / "libs",
            repo_root / "MaizeModel",
            repo_root,
            Path(r"C:\msys64\ucrt64\bin"),
            Path(r"C:\msys64\mingw64\bin"),
        ]

        for env_key in ("MINGW_PREFIX", "MSYSTEM_PREFIX"):
            env_value = os.environ.get(env_key)
            if env_value:
                dirs.append(Path(env_value) / "bin")

        seen = set()
        for dll_dir in dirs:
            if not dll_dir.exists():
                continue
            resolved = str(dll_dir.resolve())
            if resolved in seen:
                continue
            seen.add(resolved)
            if resolved in Maize._dll_dirs_added:
                continue
            try:
                cookie = os.add_dll_directory(resolved)
            except OSError:
                continue
            Maize._dll_dirs_added[resolved] = cookie
            self._dll_dir_handles.append(cookie)

    def _bind_api(self):
        lib = self._lib
        lib.maize_create.restype = ctypes.POINTER(_Handle)
        lib.maize_destroy.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_set_species.argtypes = [ctypes.POINTER(_Handle), ctypes.c_char_p]
        lib.maize_set_species.restype = ctypes.c_int
        lib.maize_reset_plant.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_set_global_settings.argtypes = [
            ctypes.POINTER(_Handle),
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_int,
            ctypes.c_float,
        ]
        lib.maize_set_global_settings.restype = ctypes.c_int
        lib.maize_add_tiller.argtypes = [ctypes.POINTER(_Handle), ctypes.POINTER(TillerDesc)]
        lib.maize_add_tiller.restype = ctypes.c_int
        lib.maize_set_tiller.argtypes = [ctypes.POINTER(_Handle), ctypes.c_int, ctypes.POINTER(TillerDesc)]
        lib.maize_set_tiller.restype = ctypes.c_int
        lib.maize_tiller_count.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_tiller_count.restype = ctypes.c_int
        lib.maize_add_leaf.argtypes = [ctypes.POINTER(_Handle), ctypes.c_int, ctypes.POINTER(LeafDesc)]
        lib.maize_add_leaf.restype = ctypes.c_int
        lib.maize_set_leaf.argtypes = [ctypes.POINTER(_Handle), ctypes.c_int, ctypes.c_int, ctypes.POINTER(LeafDesc)]
        lib.maize_set_leaf.restype = ctypes.c_int
        lib.maize_leaf_count.argtypes = [ctypes.POINTER(_Handle), ctypes.c_int]
        lib.maize_leaf_count.restype = ctypes.c_int
        lib.maize_rebuild_geometry.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_rebuild_geometry.restype = ctypes.c_int
        lib.maize_get_triangle_count.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_get_triangle_count.restype = ctypes.c_int
        lib.maize_leaf_spline_trait_count.argtypes = [ctypes.POINTER(_Handle)]
        lib.maize_leaf_spline_trait_count.restype = ctypes.c_int
        lib.maize_get_leaf_spline_traits.argtypes = [ctypes.POINTER(_Handle), ctypes.c_int, ctypes.POINTER(LeafSplineTraits)]
        lib.maize_get_leaf_spline_traits.restype = ctypes.c_int
        lib.maize_save_obj.argtypes = [
            ctypes.POINTER(_Handle),
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
        ]
        lib.maize_save_obj.restype = ctypes.c_int
        lib.maize_load_xml.argtypes = [ctypes.POINTER(_Handle), ctypes.c_char_p]
        lib.maize_load_xml.restype = ctypes.c_int
        lib.maize_save_xml.argtypes = [ctypes.POINTER(_Handle), ctypes.c_char_p]
        lib.maize_save_xml.restype = ctypes.c_int
        lib.maize_default_tiller_desc.argtypes = [ctypes.POINTER(TillerDesc)]
        lib.maize_default_leaf_desc.argtypes = [ctypes.POINTER(LeafDesc)]

        # ABI introspection (added in C-API 1.0.0). Bind only if the symbol
        # exists so that older libraries still load and emit a warning instead
        # of erroring out.
        for name in (
            "maize_c_api_version_major",
            "maize_c_api_version_minor",
            "maize_c_api_version_patch",
            "maize_leaf_desc_size",
            "maize_tiller_desc_size",
            "maize_leaf_spline_traits_size",
        ):
            fn = getattr(lib, name, None)
            if fn is not None:
                fn.argtypes = []
                fn.restype = ctypes.c_int

    def _check_api_compatibility(self) -> None:
        """Validate the loaded library's ABI version + struct sizes against expectations."""
        lib = self._lib
        getter = getattr(lib, "maize_c_api_version_major", None)
        if getter is None:
            warnings.warn(
                "Loaded maize_c_api shared library does not export ABI version "
                "symbols (likely pre-1.0.0). Struct layout cannot be validated; "
                "rebuild the C++ library to silence this warning.",
                stacklevel=3,
            )
            return

        lib_major = lib.maize_c_api_version_major()
        lib_minor = lib.maize_c_api_version_minor()
        lib_patch = lib.maize_c_api_version_patch()
        expected_major, expected_minor, _ = EXPECTED_C_API_VERSION

        if lib_major != expected_major:
            raise CApiVersionMismatch(
                f"maize_c_api major version mismatch: "
                f"library reports {lib_major}.{lib_minor}.{lib_patch}, "
                f"Python wrapper expects {expected_major}.{expected_minor}.x. "
                f"Rebuild the C++ library to match (MaizeModel/build_api.sh)."
            )
        if lib_minor < expected_minor:
            warnings.warn(
                f"maize_c_api minor version {lib_major}.{lib_minor}.{lib_patch} is older than "
                f"the Python wrapper's expected {expected_major}.{expected_minor}.x; "
                f"newer wrapper features may be unavailable. Consider rebuilding the C++ library.",
                stacklevel=3,
            )

        # Struct-size sanity check — catches silent layout changes (e.g. someone
        # adds a field on the C++ side without bumping the ABI version).
        size_checks = [
            ("MaizeLeafDesc", lib.maize_leaf_desc_size(), ctypes.sizeof(LeafDesc)),
            ("MaizeTillerDesc", lib.maize_tiller_desc_size(), ctypes.sizeof(TillerDesc)),
            ("MaizeLeafSplineTraits", lib.maize_leaf_spline_traits_size(), ctypes.sizeof(LeafSplineTraits)),
        ]
        mismatches = [
            f"{name}: C={c_size}B, Python ctypes={py_size}B"
            for name, c_size, py_size in size_checks
            if c_size != py_size
        ]
        if mismatches:
            raise CApiVersionMismatch(
                "maize_c_api struct layout mismatch between Python wrapper and shared library:\n  "
                + "\n  ".join(mismatches)
                + "\nRebuild the C++ library against the current header (MaizeModel/build_api.sh)."
            )

    def reset(self):
        """Reset the underlying plant to an empty state (no tillers, no leaves)."""
        self._lib.maize_reset_plant(self._handle)

    def set_species(self, species: str):
        """Set the plant's species name (free-form identifier).

        Returns ``True`` on success.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_species(self._handle, species.encode("utf-8")))

    def set_global_settings(self, density_v: float, density_u: float, stem_density_scale: float, stem_row_override: int, stem_ribbon_spacing: float):
        """Configure mesh tessellation density and stem ribbon parameters.

        Parameters
        ----------
        density_v, density_u:
            Vertex density along the leaf longitudinal (V) and transverse (U) axes.
        stem_density_scale:
            Multiplier applied to the stem mesh's vertex density.
        stem_row_override:
            If non-zero, forces the stem to use this number of rows of vertices.
        stem_ribbon_spacing:
            Spacing between successive stem ribbons (meters).

        Returns ``True`` on success.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(
            self._lib.maize_set_global_settings(
                self._handle,
                float(density_v),
                float(density_u),
                float(stem_density_scale),
                int(stem_row_override),
                float(stem_ribbon_spacing),
            )
        )

    def default_tiller(self) -> TillerDesc:
        """Return a ``TillerDesc`` populated with the C++ engine's default values."""
        desc = TillerDesc()
        self._lib.maize_default_tiller_desc(ctypes.byref(desc))
        return desc

    def default_leaf(self) -> LeafDesc:
        """Return a ``LeafDesc`` populated with the C++ engine's default values."""
        desc = LeafDesc()
        self._lib.maize_default_leaf_desc(ctypes.byref(desc))
        return desc

    def add_tiller(self, desc: Optional[TillerDesc] = None) -> int:
        """Append a tiller to the plant. Returns the new tiller's zero-based index.

        If ``desc`` is omitted, the engine's default tiller parameters are used.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        desc = self.default_tiller() if desc is None else desc
        return int(self._lib.maize_add_tiller(self._handle, ctypes.byref(desc)))

    def set_tiller(self, index: int, desc: TillerDesc) -> bool:
        """Replace the tiller at ``index`` with ``desc``. Returns ``True`` on success."""
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_tiller(self._handle, int(index), ctypes.byref(desc)))

    def tiller_count(self) -> int:
        """Return the number of tillers currently attached to the plant."""
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_tiller_count(self._handle))

    def add_leaf(self, tiller_index: int, desc: Optional[LeafDesc] = None) -> int:
        """Append a leaf to the tiller at ``tiller_index``. Returns the new leaf index.

        If ``desc`` is omitted, the engine's default leaf parameters are used.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        desc = self.default_leaf() if desc is None else desc
        return int(self._lib.maize_add_leaf(self._handle, int(tiller_index), ctypes.byref(desc)))

    def set_leaf(self, tiller_index: int, leaf_index: int, desc: LeafDesc) -> bool:
        """Replace the leaf at ``(tiller_index, leaf_index)``. Returns ``True`` on success."""
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_leaf(self._handle, int(tiller_index), int(leaf_index), ctypes.byref(desc)))

    def leaf_count(self, tiller_index: int) -> int:
        """Return the number of leaves currently on the tiller at ``tiller_index``."""
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_leaf_count(self._handle, int(tiller_index)))

    def rebuild(self) -> bool:
        """Regenerate the plant mesh from the current descriptor.

        Must be called before :meth:`triangle_count`, :meth:`save_obj`, or
        :meth:`get_leaf_spline_traits` reflect any changes made via
        :meth:`add_leaf`, :meth:`set_leaf`, or descriptor mutators.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_rebuild_geometry(self._handle))

    def triangle_count(self) -> int:
        """Return the total triangle count of the most recently rebuilt mesh."""
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_get_triangle_count(self._handle))

    def leaf_spline_trait_count(self) -> int:
        """Return the number of leaves for which spline traits are available.

        Equals the total leaf count across all tillers after :meth:`rebuild`.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_leaf_spline_trait_count(self._handle))

    def get_leaf_spline_traits(self, leaf_index: int):
        """Return spline-derived traits for the leaf at the global ``leaf_index``.

        The returned dict has the keys

        - ``leaf_index`` (int)
        - ``connection_point`` (``{x, y, z}`` in world space, meters)
        - ``tip_position`` (``{x, y, z}`` in world space, meters)
        - ``base_tangent`` (``{x, y, z}`` direction vector at the leaf base)
        - ``leaf_length`` (float, meters)
        - ``leaf_angle`` (``{azimuth_deg, inclination_deg}``)

        See :func:`phenoframe.compute_traits_from_descriptor` for the equivalent
        pure-Python computation directly from a descriptor XML.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        out = LeafSplineTraits()
        ok = bool(self._lib.maize_get_leaf_spline_traits(self._handle, int(leaf_index), ctypes.byref(out)))
        if not ok:
            raise RuntimeError(f"Failed to get spline traits for leaf index {leaf_index}")
        return {
            "leaf_index": int(leaf_index),
            "connection_point": {"x": float(out.connection_x), "y": float(out.connection_y), "z": float(out.connection_z)},
            "tip_position": {"x": float(out.tip_x), "y": float(out.tip_y), "z": float(out.tip_z)},
            "base_tangent": {"x": float(out.base_tangent_x), "y": float(out.base_tangent_y), "z": float(out.base_tangent_z)},
            "leaf_length": float(out.leaf_length),
            "leaf_angle": {
                "azimuth_deg": float(out.azimuth_deg),
                "inclination_deg": float(out.inclination_deg),
            },
        }

    def save_obj(
        self,
        path: str,
        leaf_texture_path: str = "plants/maize_leaf.png",
        stem_texture_path: str = "plants/maize_stem_texture.png",
        include_stem: bool = True,
        include_leaves: bool = True,
        separate_leaves: bool = False,
    ) -> bool:
        """Export the most recently rebuilt mesh as Wavefront OBJ + MTL.

        Parameters
        ----------
        path:
            Output path. ``.obj`` is appended automatically when missing
            (unless ``separate_leaves=True``, in which case ``path`` is treated
            as a directory).
        leaf_texture_path, stem_texture_path:
            Texture image paths recorded in the MTL file. Pass empty strings
            to omit textures.
        include_stem, include_leaves:
            Toggles for which mesh layers to write.
        separate_leaves:
            When ``True``, each leaf is written to its own OBJ file inside
            ``path/`` (treated as a directory).

        Returns ``True`` on success. The MTL filename is normalized to
        ``maize.mtl`` post-write for consistency.
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        if (not separate_leaves) and (not path.lower().endswith(".obj")):
            path = path + ".obj"
        leaf_ptr = leaf_texture_path.encode("utf-8") if leaf_texture_path else None
        stem_ptr = stem_texture_path.encode("utf-8") if stem_texture_path else None
        ok = bool(
            self._lib.maize_save_obj(
                self._handle,
                path.encode("utf-8"),
                leaf_ptr,
                stem_ptr,
                int(bool(include_stem)),
                int(bool(include_leaves)),
                int(bool(separate_leaves)),
            )
        )
        if ok and not separate_leaves:
            self._normalize_obj_mtllib(path)
        return ok

    @staticmethod
    def _normalize_obj_mtllib(obj_path: str) -> None:
        obj_file = Path(obj_path)
        if not obj_file.exists():
            return

        try:
            lines = obj_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            return

        desired_mtl = "maize.mtl"
        current_mtl = None
        mtllib_index = None

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.lower().startswith("mtllib "):
                mtllib_index = i
                current_mtl = stripped.split(None, 1)[1].strip() if len(stripped.split(None, 1)) > 1 else None
                break

        if current_mtl and current_mtl.lower() != desired_mtl.lower():
            old_mtl = obj_file.parent / current_mtl
            new_mtl = obj_file.parent / desired_mtl
            if old_mtl.exists() and not new_mtl.exists():
                try:
                    shutil.copyfile(old_mtl, new_mtl)
                except OSError:
                    pass

        if mtllib_index is not None:
            lines[mtllib_index] = f"mtllib {desired_mtl}"
        else:
            insert_at = 1 if lines and lines[0].lstrip().startswith("#") else 0
            lines.insert(insert_at, f"mtllib {desired_mtl}")

        try:
            obj_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except OSError:
            pass

    def load_xml(self, path: str, rebuild: bool = True) -> bool:
        """Load a plant descriptor XML and (optionally) rebuild the mesh.

        See ``DESCRIPTOR_FORMAT.md`` for the full schema. Returns ``True`` on
        successful load (and successful rebuild, when ``rebuild=True``).
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        ok = bool(self._lib.maize_load_xml(self._handle, path.encode("utf-8")))
        if ok and rebuild:
            ok = self.rebuild()
        return ok

    def save_xml(self, path: str) -> bool:
        """Serialize the current plant descriptor to XML at ``path``.

        ``.xml`` is appended automatically when missing. Returns ``True`` on
        success. The output round-trips through :meth:`load_xml` losslessly
        for parameters present in the C++ struct (override control points are
        preserved when populated).
        """
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        if not path.lower().endswith(".xml"):
            path = path + ".xml"
        return bool(self._lib.maize_save_xml(self._handle, path.encode("utf-8")))


@dataclass
class Leaf:
    """Pythonic builder for a leaf descriptor.

    Mirrors :class:`LeafDesc` (the underlying ctypes struct) but with default
    values and Python-native types. Convert to a ctypes payload with
    :meth:`to_ctype` before passing to :meth:`Maize.add_leaf`.

    Lengths are in meters; angles are in degrees. See
    ``DESCRIPTOR_FORMAT.md`` for per-field semantics. Note ``splinePoints``
    is clamped to ``[4, 40]`` by the engine on load.
    """

    id: int = 0
    distance: float = 0.0
    leaf_length: float = 0.7
    leaf_width: float = 0.08
    azimuth_deg: float = 0.0
    leaf_angle: float = 45.0
    droopiness: float = 0.5
    stem_inclination_deg: float = 0.0
    spline_points: int = 40
    width_taper: float = 2.5
    leaf_twist: float = 0.0
    leaf_curl: float = 0.0
    wave_l_amp: float = 0.0
    wave_l_freq: float = 0.0
    wave_l_phase: float = 0.0
    wave_r_amp: float = 0.0
    wave_r_freq: float = 0.0
    wave_r_phase: float = 0.0
    surface_noise_amp: float = 0.01
    surface_noise_freq: float = 8.0
    midrib_tip_taper_start: float = 0.75
    midrib_texture_strength: float = 0.35
    midrib_width: float = 0.075
    ligule_wrap_length_scale: float = 6.8
    ligule_unfold_sharpness: float = 2.2
    sheath_outer_scale: float = 1.12
    use_ctrl_overrides: bool = False

    def to_ctype(self) -> LeafDesc:
        """Convert to a ctypes :class:`LeafDesc` payload for the C ABI."""
        return LeafDesc(
            self.id,
            float(self.distance),
            float(self.leaf_length),
            float(self.leaf_width),
            float(self.azimuth_deg),
            float(self.leaf_angle),
            float(self.droopiness),
            float(self.stem_inclination_deg),
            int(self.spline_points),
            float(self.width_taper),
            float(self.leaf_twist),
            float(self.leaf_curl),
            float(self.wave_l_amp),
            float(self.wave_l_freq),
            float(self.wave_l_phase),
            float(self.wave_r_amp),
            float(self.wave_r_freq),
            float(self.wave_r_phase),
            float(self.surface_noise_amp),
            float(self.surface_noise_freq),
            float(self.midrib_tip_taper_start),
            float(self.midrib_texture_strength),
            float(self.midrib_width),
            float(self.ligule_wrap_length_scale),
            float(self.ligule_unfold_sharpness),
            float(self.sheath_outer_scale),
            int(self.use_ctrl_overrides),
        )


@dataclass
class Tiller:
    """Pythonic builder for a tiller (primary stem) descriptor.

    Mirrors :class:`TillerDesc` with default values and Python-native types.
    Convert to a ctypes payload with :meth:`to_ctype` before passing to
    :meth:`Maize.add_tiller`. See ``DESCRIPTOR_FORMAT.md`` for per-field
    semantics.
    """

    type: str = "main"
    radius: float = 0.015
    alpha_deg: float = 0.0
    stem_shrink: float = 0.001
    azimuth_deg: float = 180.0
    azimuth_noise: float = 45.0
    random_seed: int = 1337

    def to_ctype(self) -> TillerDesc:
        """Convert to a ctypes :class:`TillerDesc` payload for the C ABI."""
        return TillerDesc(
            self.type.encode("utf-8"),
            float(self.radius),
            float(self.alpha_deg),
            float(self.stem_shrink),
            float(self.azimuth_deg),
            float(self.azimuth_noise),
            int(self.random_seed),
        )


def build_example(output_obj: str = "plants/example") -> bool:
    """Build a small five-leaf test plant and export it to OBJ.

    Useful as a smoke test that the C++ engine and texture lookups work.
    Returns ``True`` on successful export. Texture paths are relative
    (``plants/maize_leaf.png``, ``plants/maize_stem_texture.png``) and assume
    the working directory contains a ``plants/`` folder.
    """
    with Maize() as gen:
        gen.reset()
        gen.set_species("maize")
        gen.set_global_settings(density_v=200.0, density_u=200.0, stem_density_scale=0.15, stem_row_override=0, stem_ribbon_spacing=0.15)
        t = gen.add_tiller(Tiller().to_ctype())
        leaves = [
            Leaf(id=0, distance=0.24, leaf_length=0.45, leaf_width=0.15, leaf_angle=35, droopiness=-30.0, width_taper=5.0, leaf_curl=0.20),
            Leaf(id=1, distance=0.12, leaf_length=0.55, leaf_width=0.15, leaf_angle=38, droopiness=-20.0, width_taper=5.0, leaf_curl=0.20),
            Leaf(id=2, distance=0.20, leaf_length=0.65, leaf_width=0.15, leaf_angle=42, droopiness=-15.0, width_taper=5.0, leaf_curl=0.20),
            Leaf(id=3, distance=0.30, leaf_length=0.55, leaf_width=0.15, leaf_angle=46, droopiness=-10.0, width_taper=5.0, leaf_curl=0.20),
            Leaf(id=4, distance=0.40, leaf_length=0.50, leaf_width=0.15, leaf_angle=50, droopiness=0, width_taper=5.0, leaf_curl=0.20),
        ]
        for leaf in leaves:
            gen.add_leaf(t, leaf.to_ctype())
        gen.rebuild()
        return gen.save_obj(output_obj, leaf_texture_path="plants/maize_leaf.png", stem_texture_path="plants/maize_stem_texture.png")


def from_xml_to_obj(
    xml_path: str,
    output_obj: str,
    leaf_texture_path: Optional[str] = None,
    stem_texture_path: Optional[str] = None,
    include_stem: bool = True,
    include_leaves: bool = True,
    separate_leaves: bool = False,
    density_v: float = 200.0,
    density_u: float = 200.0,
    stem_density_scale: float = 0.15,
) -> bool:
    """Load a plant descriptor XML, rebuild geometry, and export it as OBJ.

    Convenience wrapper that opens a :class:`Maize`, calls :meth:`load_xml`
    + :meth:`save_obj`, and closes the handle. See :meth:`Maize.save_obj`
    for the meaning of the ``include_*`` and ``separate_leaves`` flags.

    The ``density_v``, ``density_u``, and ``stem_density_scale`` parameters
    control mesh tessellation resolution. Defaults are set to the maximum
    quality (200×200) matching the C++ engine's constructor defaults.

    Returns ``True`` on success.
    """
    with Maize() as gen:
        gen.set_global_settings(
            density_v=density_v,
            density_u=density_u,
            stem_density_scale=stem_density_scale,
            stem_row_override=0,
            stem_ribbon_spacing=0.15,
        )
        if not gen.load_xml(xml_path, rebuild=True):
            return False
        return gen.save_obj(
            output_obj,
            leaf_texture_path,
            stem_texture_path,
            include_stem=include_stem,
            include_leaves=include_leaves,
            separate_leaves=separate_leaves,
        )
