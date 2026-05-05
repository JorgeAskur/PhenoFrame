"""
Compute maize leaf traits from a descriptor XML by rebuilding center splines in Python.

Trait definitions:
1) Leaf-stem connection point: center spline point at u=0.
2) Leaf tip position: center spline point at u=1.
3) Leaf length: arc length of the center spline polyline.
4) Leaf angle:
   - Base tangent t = C(u1) - C(0), using the first two center spline points.
   - Azimuth (deg) = atan2(tz, tx), normalized to [0, 360).
   - Inclination (deg) = atan2(ty, sqrt(tx^2 + tz^2)).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
import xml.etree.ElementTree as ET


Vec3 = Tuple[float, float, float]


def _v_add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _v_sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _v_scale(v: Vec3, s: float) -> Vec3:
    return (v[0] * s, v[1] * s, v[2] * s)


def _v_dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _v_cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _v_len(v: Vec3) -> float:
    return math.sqrt(_v_dot(v, v))


def _v_norm(v: Vec3) -> Vec3:
    n = _v_len(v)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    return _v_scale(v, 1.0 / n)


def _v_lerp(a: Vec3, b: Vec3, t: float) -> Vec3:
    return _v_add(_v_scale(a, 1.0 - t), _v_scale(b, t))


def _clamp(x: float, a: float, b: float) -> float:
    return max(a, min(b, x))


def _radians(deg: float) -> float:
    return deg * math.pi / 180.0


def _dist(a: Vec3, b: Vec3) -> float:
    return _v_len(_v_sub(a, b))


def _hash_to_float_minus1_1(seed: int) -> float:
    seed = seed & 0xFFFFFFFF
    seed = ((seed ^ 61) ^ (seed >> 16)) & 0xFFFFFFFF
    seed = (seed + (seed << 3)) & 0xFFFFFFFF
    seed = (seed ^ (seed >> 4)) & 0xFFFFFFFF
    seed = (seed * 0x27D4EB2D) & 0xFFFFFFFF
    seed = (seed ^ (seed >> 15)) & 0xFFFFFFFF
    return (float(seed & 0xFFFFFF) / (0xFFFFFF / 2.0)) - 1.0


def _open_uniform_knots(n_ctrl: int, degree: int) -> List[float]:
    n = n_ctrl - 1
    m = n + degree + 1
    knots = [0.0] * (m + 1)
    for j in range(degree + 1):
        knots[j] = 0.0
    for j in range(m - degree, m + 1):
        knots[j] = 1.0
    for j in range(degree + 1, m - degree):
        knots[j] = float(j - degree) / float(m - 2 * degree)
    return knots


def _find_span(n: int, degree: int, u: float, knots: Sequence[float]) -> int:
    if u >= knots[n + 1]:
        return n
    low = degree
    high = n + 1
    while high - low > 1:
        mid = (low + high) // 2
        if u < knots[mid]:
            high = mid
        else:
            low = mid
    return low


def _de_boor(u: float, ctrl: Sequence[Vec3]) -> Vec3:
    if len(ctrl) < 2:
        return ctrl[0] if ctrl else (0.0, 0.0, 0.0)
    if len(ctrl) < 4:
        uu = _clamp(u, 0.0, 1.0) * (len(ctrl) - 1)
        i = int(_clamp(uu, 0.0, len(ctrl) - 2))
        t = uu - i
        return _v_lerp(ctrl[i], ctrl[i + 1], t)

    degree = 3
    n = len(ctrl) - 1
    uu = _clamp(u, 0.0, 1.0)
    knots = _open_uniform_knots(len(ctrl), degree)
    span = _find_span(n, degree, uu, knots)
    d = [ctrl[span - degree + j] for j in range(degree + 1)]

    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            denom = knots[span + 1 + j - r] - knots[span - degree + j]
            a = (uu - knots[span - degree + j]) / denom if abs(denom) > 1e-8 else 0.0
            d[j] = _v_add(_v_scale(d[j - 1], 1.0 - a), _v_scale(d[j], a))
    return d[degree]


def _float_attr(elem: ET.Element, key: str, default: float) -> float:
    raw = elem.get(key)
    return float(raw) if raw is not None else default


def _int_attr(elem: ET.Element, key: str, default: int) -> int:
    raw = elem.get(key)
    return int(raw) if raw is not None else default


def _parse_spline_attributes(parent: ET.Element, name: str) -> List[Vec3]:
    spline_elem = parent.find(name)
    if spline_elem is None:
        return []
    count = _int_attr(spline_elem, "count", 0)
    points: List[Vec3] = []
    for i in range(count):
        x = _float_attr(spline_elem, f"p{i}x", 0.0)
        y = _float_attr(spline_elem, f"p{i}y", 0.0)
        z = _float_attr(spline_elem, f"p{i}z", 0.0)
        points.append((x, y, z))
    return points


@dataclass
class LeafDesc:
    id: int = 0
    distance: float = 0.0
    leaf_length: float = 0.7
    leaf_width: float = 0.08
    leaf_angle: float = 45.0
    droopiness: float = 0.5
    stem_inclination_deg: float = 0.0
    spline_points: int = 4
    use_ctrl_overrides: bool = False
    ctrl_center_override: List[Vec3] = field(default_factory=list)


@dataclass
class TillerDesc:
    type: str = "main"
    radius: float = 0.02
    alpha_deg: float = 0.0
    stem_shrink: float = 0.001
    azimuth_deg: float = 180.0
    azimuth_noise: float = 45.0
    random_seed: int = 1337
    use_stem_ctrl_overrides: bool = False
    stem_ctrl_override: List[Vec3] = field(default_factory=list)
    leaves: List[LeafDesc] = field(default_factory=list)


def _parse_descriptor(xml_path: Path) -> List[TillerDesc]:
    tree = ET.parse(xml_path)
    root = tree.getroot()
    species = root.find("species")
    if species is None:
        raise ValueError("Invalid descriptor XML: missing <species>")

    tillers: List[TillerDesc] = []
    for t in species.findall("Tiller"):
        td = TillerDesc(
            type=t.get("type", "main"),
            radius=_float_attr(t, "radius", 0.02),
            alpha_deg=_float_attr(t, "alpha", 0.0),
            stem_shrink=_float_attr(t, "stemShrink", 0.001),
            azimuth_deg=_float_attr(t, "beta", 180.0),
            azimuth_noise=45.0,
            random_seed=_int_attr(t, "randomSeed", 1337),
        )

        if _int_attr(t, "useStemCtrlOverrides", 0) != 0:
            td.stem_ctrl_override = _parse_spline_attributes(t, "stemCtrl")
            td.use_stem_ctrl_overrides = len(td.stem_ctrl_override) >= 2

        leaves_elem = t.find("leaves")
        if leaves_elem is not None:
            for lf in leaves_elem.findall("leaf"):
                leaf = LeafDesc(
                    id=_int_attr(lf, "id", 0),
                    distance=_float_attr(lf, "distance", 0.0),
                    leaf_length=_float_attr(lf, "leafLength", 0.7),
                    leaf_width=_float_attr(lf, "leafWidth", 0.08),
                    leaf_angle=_float_attr(lf, "leafAngle", 45.0),
                    droopiness=_float_attr(lf, "droopiness", 0.5),
                    stem_inclination_deg=_float_attr(lf, "stemInclinationDeg", 0.0),
                    spline_points=_int_attr(lf, "splinePoints", 4),
                )
                leaf.spline_points = max(4, min(leaf.spline_points, 40))
                if _int_attr(lf, "useCtrlOverrides", 0) != 0:
                    leaf.ctrl_center_override = _parse_spline_attributes(lf, "ctrlCenter")
                    leaf.use_ctrl_overrides = len(leaf.ctrl_center_override) > 0
                td.leaves.append(leaf)

        tillers.append(td)

    return tillers


def _stem_state(td: TillerDesc) -> Tuple[List[Vec3], List[Vec3], List[float], List[float], float]:
    stem_ctrl: List[Vec3] = [(0.0, 0.0, 0.0)]
    leaf_stem_pos: List[Vec3] = []
    leaf_stem_len: List[float] = []
    leaf_inclinations: List[float] = []

    cumulative_incl = 0.0
    cumulative_len = 0.0

    if td.use_stem_ctrl_overrides and len(td.stem_ctrl_override) >= 2:
        stem_ctrl = list(td.stem_ctrl_override)

        total_len = 0.0
        for i in range(1, len(stem_ctrl)):
            total_len += _dist(stem_ctrl[i], stem_ctrl[i - 1])
        cumulative_len = total_len

        for i, leaf in enumerate(td.leaves):
            ctrl_idx = min(i + 1, len(stem_ctrl) - 1)
            leaf_stem_pos.append(stem_ctrl[ctrl_idx])

            sub_len = 0.0
            for s in range(1, ctrl_idx + 1):
                sub_len += _dist(stem_ctrl[s], stem_ctrl[s - 1])
            leaf_stem_len.append(sub_len)

            cumulative_incl += leaf.stem_inclination_deg
            leaf_inclinations.append(cumulative_incl)
    else:
        stem_pos: Vec3 = (0.0, 0.0, 0.0)
        for leaf in td.leaves:
            d = max(0.0, leaf.distance)
            cumulative_len += d
            cumulative_incl += leaf.stem_inclination_deg
            inc = _radians(cumulative_incl)
            direction = (math.sin(inc), math.cos(inc), 0.0)
            stem_pos = _v_add(stem_pos, _v_scale(direction, d))
            stem_ctrl.append(stem_pos)
            leaf_stem_pos.append(stem_pos)
            leaf_stem_len.append(cumulative_len)
            leaf_inclinations.append(cumulative_incl)

    stem_length = max(0.001, cumulative_len)
    return stem_ctrl, leaf_stem_pos, leaf_stem_len, leaf_inclinations, stem_length


def _stem_pos_at(stem_ctrl: Sequence[Vec3], t: float) -> Vec3:
    return _de_boor(_clamp(t, 0.0, 1.0), stem_ctrl)


def _stem_tangent_at(stem_ctrl: Sequence[Vec3], t: float) -> Vec3:
    dt = 1.0 / max(2, len(stem_ctrl))
    p0 = _stem_pos_at(stem_ctrl, _clamp(t - dt, 0.0, 1.0))
    p1 = _stem_pos_at(stem_ctrl, _clamp(t + dt, 0.0, 1.0))
    tan = _v_sub(p1, p0)
    if _v_len(tan) < 1e-6:
        return (0.0, 1.0, 0.0)
    return _v_norm(tan)


def _build_leaf_center_local(leaf: LeafDesc, stem_radius_at_exit: float) -> List[Vec3]:
    if leaf.use_ctrl_overrides and leaf.ctrl_center_override:
        return list(leaf.ctrl_center_override)

    if leaf.spline_points < 2:
        return [(0.0, stem_radius_at_exit, 0.0)]

    start_dir = (0.0, 1.0, 0.0)
    angle_rad = _radians(leaf.leaf_angle)
    target_dir = (0.0, math.cos(angle_rad), math.sin(angle_rad))

    out_center: List[Vec3] = [(0.0, stem_radius_at_exit, 0.0)]
    current_pos: Vec3 = (0.0, stem_radius_at_exit, 0.0)

    segment_len = leaf.leaf_length / float(leaf.spline_points - 1)
    gravity_factor = leaf.droopiness * 5.0
    t_p1 = 1.0 / float(max(1, leaf.spline_points - 1))

    blend_end = (stem_radius_at_exit * 2.0) / max(leaf.leaf_length, 1e-6)
    blend_end = _clamp(blend_end, 0.001, 1.0)

    base_dir_first = _v_norm(target_dir)
    rot_amt_first = _radians(gravity_factor * (t_p1 * t_p1))
    cy, cz = base_dir_first[1], base_dir_first[2]
    first_dir = (base_dir_first[0], cy * math.cos(rot_amt_first) - cz * math.sin(rot_amt_first), cy * math.sin(rot_amt_first) + cz * math.cos(rot_amt_first))
    first_dir = _v_norm(first_dir)
    first_pos = _v_add(current_pos, _v_scale(first_dir, segment_len))

    if leaf.spline_points > 2:
        for k in range(1, 3):
            t_seg = float(k) / 3.0
            out_center.append(_v_lerp(current_pos, first_pos, t_seg))

    current_dir = first_dir if _v_len(first_dir) > 1e-8 else start_dir
    current_pos = first_pos
    out_center.append(current_pos)

    for i in range(2, leaf.spline_points):
        t = float(i) / float(leaf.spline_points - 1)
        t_blend = _clamp(t / blend_end, 0.0, 1.0)
        _ = t_blend * t_blend * (3.0 - 2.0 * t_blend)  # kept for parity with C++ flow

        base_dir = _v_norm(target_dir)
        rot_amt = _radians(gravity_factor * (t * t))
        cy, cz = base_dir[1], base_dir[2]
        target_step_dir = (base_dir[0], cy * math.cos(rot_amt) - cz * math.sin(rot_amt), cy * math.sin(rot_amt) + cz * math.cos(rot_amt))
        target_step_dir = _v_norm(target_step_dir)

        current_dir = _v_norm(_v_lerp(current_dir, target_step_dir, 0.5))
        current_pos = _v_add(current_pos, _v_scale(current_dir, segment_len))
        out_center.append(current_pos)

    return out_center


def _tilt_xy(p: Vec3, alpha_deg: float) -> Vec3:
    tilt = _radians(alpha_deg)
    ct = math.cos(tilt)
    st = math.sin(tilt)
    return (p[0] * ct - p[1] * st, p[0] * st + p[1] * ct, p[2])


def _world_leaf_center(td: TillerDesc, leaf_index: int, stem_ctrl: Sequence[Vec3], leaf_stem_pos: Sequence[Vec3], leaf_stem_len: Sequence[float], stem_length: float) -> List[Vec3]:
    leaf = td.leaves[leaf_index]
    noise = _hash_to_float_minus1_1(td.random_seed + leaf_index)
    leaf_azimuth = td.azimuth_deg * leaf_index + noise * td.azimuth_noise
    phi = _radians(leaf_azimuth)
    rhat = (math.cos(phi), 0.0, math.sin(phi))

    t_leaf = (leaf_stem_len[leaf_index] / stem_length) if stem_length > 1e-6 else 0.0
    t_leaf = _clamp(t_leaf, 0.0, 1.0)
    global_radius = max(0.001, td.radius - td.stem_shrink * t_leaf)

    local_center = _build_leaf_center_local(leaf, global_radius)
    sheath_exit = leaf_stem_pos[leaf_index]

    blade_y = _stem_tangent_at(stem_ctrl, t_leaf)
    blade_z = rhat
    blade_x = _v_cross(blade_y, blade_z)
    if _v_len(blade_x) < 1e-6:
        blade_z = (0.0, 0.0, 1.0)
        blade_x = _v_cross(blade_y, blade_z)
    blade_x = _v_norm(blade_x)
    blade_z = _v_norm(_v_cross(blade_x, blade_y))

    world: List[Vec3] = []
    for p in local_center:
        w = _v_add(
            sheath_exit,
            _v_add(
                _v_scale(blade_x, p[0]),
                _v_add(_v_scale(blade_y, p[2]), _v_scale(blade_z, p[1])),
            ),
        )
        world.append(_tilt_xy(w, td.alpha_deg))
    return world


def _leaf_traits_from_center_spline(center: Sequence[Vec3]) -> Dict[str, object]:
    if not center:
        raise RuntimeError("Empty center spline")

    connection = center[0]
    tip = center[-1]
    length = 0.0
    for i in range(1, len(center)):
        length += _dist(center[i], center[i - 1])

    if len(center) >= 2:
        tangent = _v_sub(center[1], center[0])
    else:
        tangent = (0.0, 1.0, 0.0)

    tx, ty, tz = tangent
    azimuth = math.degrees(math.atan2(tz, tx))
    if azimuth < 0.0:
        azimuth += 360.0
    inclination = math.degrees(math.atan2(ty, math.sqrt(tx * tx + tz * tz)))

    return {
        "connection_point": {"x": connection[0], "y": connection[1], "z": connection[2]},
        "tip_position": {"x": tip[0], "y": tip[1], "z": tip[2]},
        "base_tangent": {"x": tx, "y": ty, "z": tz},
        "leaf_length": length,
        "leaf_angle": {"azimuth_deg": azimuth, "inclination_deg": inclination},
        "center_spline": center,
    }


def compute_traits_from_descriptor(xml_path: Path) -> List[Dict[str, object]]:
    tillers = _parse_descriptor(xml_path)
    traits: List[Dict[str, object]] = []
    global_leaf_index = 0

    for ti, td in enumerate(tillers):
        stem_ctrl, leaf_stem_pos, leaf_stem_len, _, stem_length = _stem_state(td)
        for li, leaf in enumerate(td.leaves):
            center_world = _world_leaf_center(td, li, stem_ctrl, leaf_stem_pos, leaf_stem_len, stem_length)
            t = _leaf_traits_from_center_spline(center_world)
            t["leaf_index"] = global_leaf_index
            t["tiller_index"] = ti
            t["leaf_id"] = leaf.id
            traits.append(t)
            global_leaf_index += 1

    return traits


def _fmt(v: float) -> str:
    return f"{v:.9g}"


def write_traits_xml(xml_path: Path, source_descriptor: Path, traits: Sequence[Dict[str, object]]) -> None:
    root = ET.Element("traits")
    root.set("source_descriptor", str(source_descriptor))

    for t in traits:
        leaf_elem = ET.SubElement(root, "leaf")
        leaf_elem.set("leaf_id", str(t["leaf_id"]))
        leaf_elem.set("tiller_id", str(t["tiller_index"]))

        angle = t["leaf_angle"]
        angle_elem = ET.SubElement(leaf_elem, "leaf_angle")
        angle_elem.set("azimuth_deg", _fmt(angle["azimuth_deg"]))
        angle_elem.set("inclination_deg", _fmt(angle["inclination_deg"]))

        length_elem = ET.SubElement(leaf_elem, "leaf_length")
        length_elem.set("value", _fmt(t["leaf_length"]))

        cp = t["connection_point"]
        cp_elem = ET.SubElement(leaf_elem, "connection_point")
        cp_elem.set("x", _fmt(cp["x"]))
        cp_elem.set("y", _fmt(cp["y"]))
        cp_elem.set("z", _fmt(cp["z"]))

        tip = t["tip_position"]
        tip_elem = ET.SubElement(leaf_elem, "tip_position")
        tip_elem.set("x", _fmt(tip["x"]))
        tip_elem.set("y", _fmt(tip["y"]))
        tip_elem.set("z", _fmt(tip["z"]))

    xml_path.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(root)
    try:
        ET.indent(tree, space="  ")
    except AttributeError:
        pass
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute spline-based maize leaf traits from descriptor XML.")
    parser.add_argument("--xml", default="plant_0.xml", help="Path to descriptor XML")
    parser.add_argument("--out-xml", default="Traits/traits.xml", help="Output traits XML path")
    args = parser.parse_args()

    xml_path = Path(args.xml)
    if not xml_path.exists():
        fallback = Path("plants") / args.xml
        if fallback.exists():
            xml_path = fallback
        else:
            raise FileNotFoundError(f"XML file not found: {args.xml}")

    traits = compute_traits_from_descriptor(xml_path)
    out_xml = Path(args.out_xml)
    write_traits_xml(out_xml, xml_path.resolve(), traits)

    print(f"Computed traits for {len(traits)} leaves.")
    print(f"Saved XML: {out_xml}")
    for t in traits:
        ang = t["leaf_angle"]
        tip = t["tip_position"]
        print(
            f"Leaf {t['leaf_index']:02d}: "
            f"len={t['leaf_length']:.4f}, "
            f"az={ang['azimuth_deg']:.2f} deg, "
            f"inc={ang['inclination_deg']:.2f} deg, "
            f"tip=({tip['x']:.4f}, {tip['y']:.4f}, {tip['z']:.4f})"
        )


if __name__ == "__main__":
    main()
