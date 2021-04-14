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

`blockMesh` uses the `system/blockMeshDict` dictionary file to generate the mesh. The mesh size is set by variables `xmin`, `xmax`, `ymin`, `ymax`, `zmin` and `zmax` which define a block enclosing the fluid domain geometry. The number of cells along `X`, `Y` and `Z` dimensions is set by variables `xcells`, `ycells` and `zcells` selected such that the criteria above are met. The boundary of the mesh, given in the list named `boundary`, consists of 5 patches `frontandback`, `inlet`, `outlet`, `lowerwall` and, `upperwall`. These patches will be discarded during the mesh generation with `snappyHexMesh` and replaced by the geometry information defined in the `system/snappyHexMeshDict` dictionary file.

In addition to the `system/blockMeshDict` dictionary file, `blockMesh` also requires the `system/controlDict` dictionary file. `blockMesh` is executed in the case folder using:

```
blockMesh | tee log.01.blockMesh
```

`blockMesh` execution generates the background hexahedral mesh in `constant/polyMesh` folder and the log file `log.blockMesh` in the case folder. The mesh information is contained in the `constant/polyMesh/boundary`, `constant/polyMesh/faces`, `constant/polyMesh/neighbour`, `constant/polyMesh/owner` and `constant/polyMesh/points` files.

The `blockMesh` mesh quality is assessed using:

```
checkMesh -allGeometry -allTopology | tee log.02.checkMesh.block
```

### snappyHexMesh

The fluid domain mesh generation with `snappyHexMesh` requires `system/fvSchemes` and `system/fvSolution` dictionary files, which are discussed in the **Solver** section. The mesh generation process is comprised of a number of steps, described below. 

#### Geometry Edge Extraction

The `surfaceFeatures` utility is used to extract the geometry edges and allow for better meshing with `snappyHexMesh` on these edges. The `system/surfaceFeaturesDict` dictionary file lists the geometry surface (stl) files in the `surfaces` list. The edges marked for extraction are those whose adjacent surfaces normal are at an angle less than the angle specified by `includedAngle` in the `system/surfaceFeaturesDict` dictionary file. The `surfaceFeatures` utility is executed in the case folder using:

```
surfaceFeatures | tee log.03.surfaceFeatures
```

The `surfaceFeatures` utility execution creates the following files for each stl file:
- A `.eMesh` file in the `constant/triSurface` folder,
- A `.extendedFeatureEdgeMesh` and a `_edgeMesh.obj` in the `constant/extendedFeatureEdgeMesh` folder.

#### Mesh Parameters

The `system/snappyHexMeshDict` dictionary file contains the `snappyHexMesh` parameters. The execution of `snappyHexMesh` consists of three steps `castellating`, `snapping` and `layering` which can be enabled or disabled as needed.

The `geometry` dictionary lists the geometry surfaces (stl) files along with their type and user defined name. The user defined name will be used in the **Model** section for the definition of IC and BC. A refinement region, named `refinementBox` is also defined.

The `castellatedMeshControls` dictionary controls the parameters for mesh refining in the `castellating` step.

- The `features` dictionary is used for edge refinement of the `*.eMesh` edges to the desired refinement level.

- The `refinementSurfaces` dictionary is used for surface  based  refinement.  Every  surface  is  specified  with  two  levels. The first level is the minimum level that every cell intersecting the surface gets refined up to. The second level is the maximum level of refinement. The `patchInfo` dictionary sets the patch type for each surface as required by the boundary condition types associated with each patch in the IC and BC dictionary files discussed in the **Model** section.

- The `resolveFeatureAngle` setting allows the edges, whose adjacent surfaces normal are at an angle higher than the value set, to be resolved. A lower value for `resolveFeatureAngle` results in a better resolution at sharp edges.

- The `refinementRegions` dictionary contains the volume  based  refinement settings for the `shell` region defined in the `geometry` dictionary. The first number of the `levels` setting represents the distance from the geometry within which all cells are refined while the second number represents the level of refinement.

- The `locationInMesh` setting identifies a location in the final mesh (inside the fluid domain) from which `snappyHexMesh` will mark and keep all connected cells.

The `snapControls` dictionary controls the parameters for mesh refining in the `snapping` step which adapts the castellated mesh to the geometry.

The `addLayersControls` dictionary controls the parameters inserting prismatic cell layers on `shell` surface. The number of layers for the `shell` surface is set by `nSurfaceLayers` to 3. Because the `relativeSizes` is set to `false`, the thickness of the first layer is set by `firstLayerThickness` to 0.004 m. The minimum thickness of any layer is set by `minThickness` to 0.004 m. The increase in size from one layer to the next is set by `expansionRatio` to 1.2.

#### Mesh Quality Parameters

The mesh quality for `snappyHexMesh` is controlled by the entries in the `meshQualityControls` dictionary in the `system/snappyHexMeshDict` dictionary file. The `meshQualityControls` dictionary uses an `include` statement to include the mesh quality settings contained in the `system/snappyHexMeshDict`.

#### Mesh Decomposition

The `decomposePar` utility is used to decompose the mesh into sub-domains, allocated to separate processors, to allow the mesh or solver execution in parallel. The `system/decomposeParDict` dictionary file contains the `decomposePar` utility parameters. The mesh is decomposed into 8 sub-domains by setting `numberOfSubdomains` to 8. The  `method` parameter is used to set the decomposition method to `simple`. For the `simple` method, the `n` parameter in the `simpleCoeffs` dictionary decomposes the mesh into 2 sub-domains along the x, y and z directions, respectively.

The `decomposePar` utility is executed in the case folder using:

```
decomposePar | tee log.04.decomposePar
```

The `decomposePar` utility execution creates a `processorX/constant/polyMesh` folder for each processor containing the sub-domain mesh files assigned to that processor.

#### Mesh Execution

The `snappyHexMesh` is executed in parallel in the case folder using:

```
mpirun -np 8 snappyHexMesh -overwrite -parallel | tee log.05.snappyHexMesh
```

#### Mesh Quality Evaluation

The `snappyHexMesh` mesh quality is assessed using:

```
mpirun -np 8 checkMesh -latestTime -allGeometry -allTopology -parallel | tee log.06.checkMesh.snappy
```

#### Mesh Reconstruction

The `reconstructParMesh` utility reads the individual processor mesh file and updates the mesh files in the `constant/polyMesh` folder. The `reconstructParMesh` utility is executed in the case folder using:

```
reconstructParMesh -latestTime -constant | tee log.07.reconstructParMesh
```

After the `snappyHexMesh` mesh is reconstructed, the individual processor mesh files can be removed using:

```
rm -rf processor* > /dev/null 2>&1
```

#### Mesh Optimization

The `renumberMesh` utility is used to reduce bandwidth and speed up computation on the generated mesh. The `renumberMesh` utility is executed in the case folder using:

```
renumberMesh -overwrite | tee log.08.renumberMesh
```

#### Mesh Visualization

The `snappyHexMesh` mesh can be visualized with Paraview using:

```
paraFoam
```

## Model


## Solver


## Execution


## Results

