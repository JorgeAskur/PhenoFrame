# GIC maize longitudinal voxel reconstructions and optimized skeletons (2024)

This dataset contains longitudinal 3D reconstructions of 80 physical maize
plants representing 20 genotypes and four biological replicates. Directory
names encode genotype, replicate, acquisition timestamp, and reconstruction
identifier. The release contains 4,171 scan directories; 3,679 have optimized
skeletons and reference leaf-angle files. The PhenoFrame paper's paired
geometric and procedural analyses use 3,317 scans that passed all analysis
filters. Repeated scans are not independent plants.

## Contents

Each scan directory may contain:

- `voxels.txt`: occupied voxel coordinates. The first line is the number of
  voxels, followed by integer `x y z` coordinates, one voxel per line.
- `optim_skeleton.txt`: optimized curve-skeleton coordinates. The first line
  is the number of points, followed by integer `x y z` coordinates.
- `angles.txt`: one tab-separated row of repeating
  `(normalized insertion height, theta, phi)` triples, ordered from lower to
  upper leaves. Angles are in degrees.
- `error.txt`: reconstruction status or fitting-error value retained from the
  original processing output.

`manifest.csv` lists every scan and records parsed identifiers, file presence,
sizes, and the counts declared in the text files. The archive intentionally
excludes `voxel_cubes.obj`, `voxel_centers.obj`, `plant_mesh.obj`, skeleton
preview PNGs, and other visualization derivatives. Those files account for
approximately 118.6 GiB and can be regenerated from the voxel and skeleton
coordinates.

## Loading with PhenoFrame

The PhenoFrame experiment loader expects a split directory layout:

```
<root>/reconstructed/<scan-id>/voxels.txt
<root>/skeletons/<scan-id>/optim_skeleton.txt
<root>/skeletons/<scan-id>/angles.txt
```

The archive itself uses one directory per scan. Avoid copying the dataset:
create `reconstructed` and `skeletons` directory links that both point to the
extracted scan directory. Then set `PHENOFRAME_VOXEL_PATH=<root>`. See the
PhenoFrame repository's `paper/DATA_PROVENANCE.md` and
`paper/scripts/run_gic_notebook.py` for commands.

## Scope and provenance

These are processed voxel reconstructions and optimized skeletons, not the raw
camera images. The archive supports reconstruction-to-descriptor, trait,
geometry, and heritability analyses. Reproducing the camera-image-to-voxel
stage requires the source image archive and acquisition calibration, which are
not included here.

The corresponding code, frozen analysis inputs, and paper reproduction scripts
are in the PhenoFrame repository:
https://github.com/JorgeAskur/PhenoFrame

## License

Data license: Creative Commons Attribution 4.0 International (CC BY 4.0).
Software in the PhenoFrame repository is distributed separately under the MIT
License.
