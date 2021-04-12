## Introduction

This guide describes a workflow for setting up and executing OpenFOAM studies and uses the simulation of water flow through a header to illustrate the approach. The workflow is comprised of a number of steps discussed individually. The purpose of the study is to evaluate the pressure, the velocity and the temperature patterns developing in the fluid domain for different conditions.

The study consists of two cases `symmetric-velocity-temperature` and `asymmetric-velocity-temperature`. The cases use the same **Geometry** and **Mesh**, but employ different Initial Condition (IC) and Boundary Conditions (BC) for the inlet patches, as presented in the **Model** section.

## Geometry

The geometry generation must create a smooth, clean and watertight geometry. A watertight geometry means a close body with no holes or overlapping surfaces. The mesh quality and, hence the solution quality, depend on the geometry. The geometry should adequately represent the fluid domain. Geometry defeaturing is used to simplify the geometry to retain only the features important for the fluid phenomena under investigation. The solid modeling application used for the geometry generation is [Onshape](https://www.onshape.com). 

The fluid domain is comprised of two inlet pipes, a common header and an outlet pipe. Each inlet pipe is comprised of a cylindrical vertical section and a cylindrical horizontal section. The horizontal section of an inlet pipe is connected to the header at a 90 degrees angle. The header consists of a horizontal pipe section. The outlet pipe consists of a horizontal section connected to the middle of the header at a 90 degrees angle.

The fluid domain geometry is created as a single watertight solid by extruding and sweeping geometry sketches. To facilitate meshing and model development, the solid geometry faces representing patches subsequently used to define fluid boundary conditions (e.g., inlet, outlet) are deleted. After the faces are deleted, the solid geometry becomes a surface geometry. The inlet and outlet patches are recreated as individual surfaces by revolving geometry sketches.

The final fluid domain geometry is comprised of 4 surface objects. The surface objects are exported from Onshape as individual stl files (text format) using the high resolution option. The 4 surface objects, `inlet1.stl`, `inlet2.stl`, `outlet.stl` and `shell.stl`, are saved to the `constant/triSurface` folder. The fluid domain geometry is shown below.

Fluid Domain Geometry and Dimensions [m]              |
:----------------------------------------------------:|
<img src="img/geometry.png"> |

## Mesh

The mesh generation uses `blockMesh` and `snappyHexMesh`.

### blockMesh

`blockMesh` is used to generate the background mesh for `snappyHexMesh`. A quality background mesh is generated when the following criteria are met:
- The background mesh must consist purely of hexahedral cells,
- The cell aspect ratio (i.e. the ratio of the longest to the shortest side of a cell) should be close to 1, at least near the stl surface,
- The cell size shall be sufficiently small to adequately resolve the smallest geometry feature,
- There must be at least one intersection of a cell edge with the stl surface.

`blockMesh` uses the `system/blockMeshDict` dictionary file to generate the mesh. The mesh size is set by variables `xmin`, `xmax`, `ymin`, `ymax`, `zmin` and `zmax` which define a block enclosing the fluid domain geometry. The number of cells along `X`, `Y` and `Z` dimensions is set by variables `xcells`, `ycells` and `zcells` for which selections are made to meet the criteria above.

In addition to the `system/blockMeshDict` dictionary file, `blockMesh` also requires the `system/controlDict` dictionary file. `blockMesh` is executed in the case folder using:

```
blockMesh | tee log.blockMesh
```

`blockMesh` execution generates the background hexahedral mesh in `constant/polyMesh` folder and the log file `log.blockMesh` in the case folder. The mesh information is contained in the `constant/polyMesh/boundary`, `constant/polyMesh/faces`, `constant/polyMesh/neighbour`, `constant/polyMesh/owner` and `constant/polyMesh/points` files.

The `blockMesh` mesh quality is assessed using:

```
checkMesh -allGeometry -allTopology | tee log.checkMesh.block
```

### snappyHexMesh

The fluid domain mesh generation with `snappyHexMesh` requires `system/fvSchemes` and `system/fvSolution` dictionary files, which are discussed in the **Solver** section. The mesh generation process is comprised of a number of steps, described below. 

#### Feature Edge Extraction

system/surfaceFeaturesDict

#### Mesh Parameters

snappyHexMeshDict

#### Mesh Quality Parameters

system/meshQualityDict

#### Mesh Partition

system/decomposeParDict

#### Mesh Execution

snappyHexMesh execution

#### Mesh Quality Evaluation

Check snappyHexMesh mesh

#### Mesh Reconstruction

Reconstruct snappyHexMesh mesh

#### Mesh Visualization

Visualize snappyHexMesh mesh with Paraview


## Model


## Solver


## Execution


## Results

