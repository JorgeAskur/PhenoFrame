
"""Python wrapper for the C++ Maize procedural plant generator."""

from __future__ import annotations

import ctypes
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
import platform
from typing import Optional, Sequence


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
    """High-level wrapper used to build and export maize plant meshes."""

    def __init__(self, library_path: Optional[str] = None):
        self._dll_dir_handles = []
        self._lib = self._load_library(library_path)
        self._bind_api()
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
        if getattr(self, "_handle", None):
            self._lib.maize_destroy(self._handle)
            self._handle = None

    # --- high-level API ---------------------------------------------------

    @staticmethod
    def _library_name():
        if platform.system() == "Windows":
            return "maize_c_api.dll"
        if platform.system() == "Darwin":
            return "libmaize_c_api.dylib"
        return "libmaize_c_api.so"

    def _load_library(self, library_path: Optional[str]) -> ctypes.CDLL:
        candidates: Sequence[Path] = []
        repo_root = Path(__file__).resolve().parent
        if library_path:
            candidates = (Path(library_path),)
        else:
            lib_name = self._library_name()
            candidates = (
                repo_root / "MaizeModel" / lib_name,
                repo_root / lib_name,
                repo_root / "MaizeModel" / "build" / lib_name,
                repo_root / "build" / lib_name,
            )
        load_errors = []
        for candidate in candidates:
            if candidate.exists():
                try:
                    self._add_windows_dll_dirs(repo_root, candidate)
                    return ctypes.cdll.LoadLibrary(str(candidate))
                except OSError as exc:
                    load_errors.append(f"{candidate}: {exc}")
                    continue
        raise FileNotFoundError(
            "Could not find compiled maize_c_api shared library.\n"
            f"Searched: {', '.join(map(str, candidates))}\n"
            + (
                f"Load errors: {' | '.join(load_errors)}\n"
                if load_errors
                else ""
            ),
            "Build it first (example: g++ -shared -O2 -std=c++17 -DMAIZE_C_API_BUILD "
            "MaizeModel/maize_c_api.cpp MaizeModel/Maize.cpp MaizeModel/Descriptor.cpp "
            "MaizeModel/vect3d.cpp MaizeModel/tinyxml2.cpp -IMaizeModel -o MaizeModel\\maize_c_api.dll)."
        )

    def _add_windows_dll_dirs(self, repo_root: Path, candidate: Path) -> None:
        if platform.system() != "Windows" or not hasattr(os, "add_dll_directory"):
            return

        dirs = [
            candidate.parent,
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
            try:
                self._dll_dir_handles.append(os.add_dll_directory(resolved))
            except OSError:
                pass

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

    def reset(self):
        self._lib.maize_reset_plant(self._handle)

    def set_species(self, species: str):
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_species(self._handle, species.encode("utf-8")))

    def set_global_settings(self, density_v: float, density_u: float, stem_density_scale: float, stem_row_override: int, stem_ribbon_spacing: float):
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
        desc = TillerDesc()
        self._lib.maize_default_tiller_desc(ctypes.byref(desc))
        return desc

    def default_leaf(self) -> LeafDesc:
        desc = LeafDesc()
        self._lib.maize_default_leaf_desc(ctypes.byref(desc))
        return desc

    def add_tiller(self, desc: Optional[TillerDesc] = None) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        desc = self.default_tiller() if desc is None else desc
        return int(self._lib.maize_add_tiller(self._handle, ctypes.byref(desc)))

    def set_tiller(self, index: int, desc: TillerDesc) -> bool:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_tiller(self._handle, int(index), ctypes.byref(desc)))

    def tiller_count(self) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_tiller_count(self._handle))

    def add_leaf(self, tiller_index: int, desc: Optional[LeafDesc] = None) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        desc = self.default_leaf() if desc is None else desc
        return int(self._lib.maize_add_leaf(self._handle, int(tiller_index), ctypes.byref(desc)))

    def set_leaf(self, tiller_index: int, leaf_index: int, desc: LeafDesc) -> bool:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_set_leaf(self._handle, int(tiller_index), int(leaf_index), ctypes.byref(desc)))

    def leaf_count(self, tiller_index: int) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_leaf_count(self._handle, int(tiller_index)))

    def rebuild(self) -> bool:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return bool(self._lib.maize_rebuild_geometry(self._handle))

    def triangle_count(self) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_get_triangle_count(self._handle))

    def leaf_spline_trait_count(self) -> int:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        return int(self._lib.maize_leaf_spline_trait_count(self._handle))

    def get_leaf_spline_traits(self, leaf_index: int):
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
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        ok = bool(self._lib.maize_load_xml(self._handle, path.encode("utf-8")))
        if ok and rebuild:
            ok = self.rebuild()
        return ok

    def save_xml(self, path: str) -> bool:
        if not self._handle:
            raise RuntimeError("Wrapper is closed")
        if not path.lower().endswith(".xml"):
            path = path + ".xml"
        return bool(self._lib.maize_save_xml(self._handle, path.encode("utf-8")))


@dataclass
class Leaf:
    id: int = 0
    distance: float = 0.0
    leaf_length: float = 0.7
    leaf_width: float = 0.08
    azimuth_deg: float = 0.0
    leaf_angle: float = 45.0
    droopiness: float = 0.5
    stem_inclination_deg: float = 0.0
    spline_points: int = 4
    width_taper: float = 1.0
    leaf_twist: float = 0.0
    leaf_curl: float = 0.0
    wave_l_amp: float = 0.0
    wave_l_freq: float = 0.0
    wave_l_phase: float = 0.0
    wave_r_amp: float = 0.0
    wave_r_freq: float = 0.0
    wave_r_phase: float = 0.0
    use_ctrl_overrides: bool = False

    def to_ctype(self) -> LeafDesc:
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
            int(self.use_ctrl_overrides),
        )


@dataclass
class Tiller:
    type: str = "main"
    radius: float = 0.015
    alpha_deg: float = 0.0
    stem_shrink: float = 0.001
    azimuth_deg: float = 180.0
    azimuth_noise: float = 45.0
    random_seed: int = 1337

    def to_ctype(self) -> TillerDesc:
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
    """Build a small test plant and export it to OBJ."""
    with Maize() as gen:
        gen.reset()
        gen.set_species("maize")
        gen.set_global_settings(density_v=120.0, density_u=120.0, stem_density_scale=0.15, stem_row_override=0, stem_ribbon_spacing=0.15)
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
) -> bool:
    """Load a plant XML descriptor and export it as OBJ."""
    with Maize() as gen:
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
