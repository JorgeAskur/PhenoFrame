# PyMaize Plant Descriptor Format

Version: 1.0 (matches `MAIZE_C_API_VERSION_MAJOR = 1`)

This document specifies the XML descriptor format consumed by both the C++
geometry engine (`maize_load_xml` in `MaizeModel/maize_c_api.h`) and the
pure-Python trait pipeline (`pymaize.compute_traits_from_descriptor`).

A formal XSD is provided alongside this document at
[`pymaize/schemas/descriptor.xsd`](pymaize/schemas/descriptor.xsd).

## Document outline

```
<plant>
  <species name="...">
    <Tiller ...>
      <stemCtrl count="N" p0x="..." p0y="..." p0z="..." ... />   (optional)
      <leaves number="M">
        <leaf ...>
          <ctrlCenter count="K" p0x="..." p0y="..." p0z="..." ... />   (optional)
        </leaf>
        ...
      </leaves>
    </Tiller>
    ...
  </species>
</plant>
```

Multiple tillers per species are supported by the data model; the current
geometry engine renders the first tiller only (single-tiller limitation —
see README).

## Coordinate frame convention

All coordinates are world-space, Y-up:

- **+Y** is the stem growth direction (vertical for an unrotated plant).
- **X / Z** form the horizontal ground plane.
- Units are meters throughout (length, distance, radius).
- Angles are in **degrees**, except where the C++ source explicitly converts
  to radians internally.

## Element reference

### `<plant>`

Root element. Has no attributes. Must contain exactly one `<species>` child.

### `<species>`

| Attribute | Type | Required | Notes |
|---|---|---|---|
| `name` | string | yes | Free-form identifier (e.g. `"Maize_Procedural"`). |

Contains zero or more `<Tiller>` children.

### `<Tiller>`

A single primary stem axis.

| Attribute | Type | Default | Meaning |
|---|---|---|---|
| `type` | string | `"main"` | Free-form (`"main"`, `"side"`, etc.). |
| `radius` | float | required | Stem radius at the base, meters. |
| `alpha` | float | `0` | Whole-tiller tilt around the world Z axis, degrees. |
| `beta` | float | `180` | Per-leaf azimuth step around the stem, degrees (`180` → alternate phyllotaxy). |
| `stemShrink` | float | `0.001` | Linear taper of stem radius from base to tip, meters. |
| `randomSeed` | unsigned int | `1337` | Seed for the deterministic per-leaf azimuth noise. |
| `useStemCtrlOverrides` | int (0/1) | `0` | If `1`, supply explicit stem control points via `<stemCtrl>`. |

May contain at most one `<stemCtrl>` and one `<leaves>` child.

### `<leaves>`

Container for `<leaf>` elements.

| Attribute | Type | Required | Notes |
|---|---|---|---|
| `number` | int | optional | Cosmetic count; the parser uses the actual `<leaf>` element count, not this attribute. |

### `<leaf>`

A single blade attached to the stem.

| Attribute | Type | Default | Meaning |
|---|---|---|---|
| `id` | int | required | Stable leaf identifier (preserved through round-trips). |
| `distance` | float | required | Internode length below this leaf, meters. |
| `leafLength` | float | required | Total blade length along the midrib, meters. Ignored when `useCtrlOverrides=1`. |
| `leafWidth` | float | required | Maximum blade width, meters. |
| `leafAngle` | float | required | Insertion angle of the blade relative to the stem, degrees. Ignored when `useCtrlOverrides=1`. |
| `droopiness` | float | `0.0` | Gravity-pull factor; negative values curl the tip downward. Ignored when `useCtrlOverrides=1`. |
| `stemInclinationDeg` | float | `0.0` | Per-leaf bend of the stem segment that immediately precedes this leaf, degrees. |
| `splinePoints` | int | `4` | Number of control points used to build the procedural midrib. **Clamped to [4, 40] on load.** The Python parser emits a `UserWarning` when clamping. |
| `widthTaper` | float | `1.0` | Width-profile shape exponent (higher = sharper tip). |
| `leafTwist` | float | `0.0` | Twist around the midrib, degrees from base to tip. |
| `leafCurl` | float | `0.0` | Lateral curl of the blade. |
| `waveLAmp` / `waveLFreq` / `waveLPhase` | float | `0.0` | Sinusoidal wave on the left edge. |
| `waveRAmp` / `waveRFreq` / `waveRPhase` | float | `0.0` | Sinusoidal wave on the right edge. |
| `useCtrlOverrides` | int (0/1) | `0` | If `1`, supply explicit midrib control points via `<ctrlCenter>` (and the procedural parameters above are bypassed). |

May contain at most one `<ctrlCenter>` child.

### `<stemCtrl>` and `<ctrlCenter>` (spline control payloads)

Both elements share the same shape: a `count` attribute giving the number of
control points, followed by `count` triples of attributes
`p{i}x` / `p{i}y` / `p{i}z` for `i` in `[0, count-1]`.

| Attribute | Type | Required | Notes |
|---|---|---|---|
| `count` | int | yes | Number of control points (≥ 2 for `<stemCtrl>`, ≥ 1 for `<ctrlCenter>`; the engine uses cubic B-spline evaluation when count ≥ 4). |
| `p{i}x`, `p{i}y`, `p{i}z` | float | yes for each `i` | Cartesian components of control point `i`. |

Coordinate space:

- `<stemCtrl>` points are in **world space**.
- `<ctrlCenter>` points are in the **leaf's local frame**, where the local
  +Y axis points outward from the stem (radial), and the local +Z axis is
  along the stem tangent at the sheath exit. The leaf's azimuthal
  rotation around the stem is applied when transforming into world space.

The engine evaluates these as cubic B-splines on an open-uniform knot vector
when `count >= 4`; for smaller counts it falls back to piecewise linear
interpolation.

## Validation

To validate a descriptor against the schema with `xmllint`:

```bash
xmllint --noout --schema pymaize/schemas/descriptor.xsd path/to/plant.xml
```

The Python pipeline does not enforce the schema at runtime — it tolerates
missing optional attributes and falls back to documented defaults. The schema
is provided as a contract for tooling and external producers.

## Sample

A minimal valid descriptor with a single tiller and a single leaf:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<plant>
  <species name="Maize_Procedural">
    <Tiller type="main" radius="0.02" alpha="0" beta="180"
            stemShrink="0.001" randomSeed="1337">
      <leaves number="1">
        <leaf id="0" distance="0.20" leafLength="0.45" leafWidth="0.07"
              leafAngle="35" droopiness="-10" stemInclinationDeg="0"
              splinePoints="8" widthTaper="5.0" leafCurl="0.2"/>
      </leaves>
    </Tiller>
  </species>
</plant>
```

A larger eight-leaf example is shipped at
[`plants/plant_0.xml`](plants/plant_0.xml).
