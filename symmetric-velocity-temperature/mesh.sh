#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Execute blockMesh
start_time="$(date -u +%s.%N)"
blockMesh | tee log.01.blockMesh
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute blockMesh: $elapsed seconds" > log.mesh

# Execute check mesh for blockMesh
start_time="$(date -u +%s.%N)"
checkMesh -allGeometry -allTopology | tee log.02.checkMesh.block
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute check mesh for blockMesh: $elapsed seconds" >> log.mesh

# Execute geometry edge extraction
start_time="$(date -u +%s.%N)"
surfaceFeatures | tee log.03.surfaceFeatures
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute edge extraction: $elapsed seconds" >> log.mesh

# Execute mesh partition
start_time="$(date -u +%s.%N)"
decomposePar | tee log.04.decomposePar
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute mesh partition: $elapsed seconds" >> log.mesh

# Execute snappyHexMesh in parallel (mesh partition required)
start_time="$(date -u +%s.%N)"
mpirun -np 8 snappyHexMesh -overwrite -parallel | tee log.05.snappyHexMesh
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute snappyHexMesh: $elapsed seconds" >> log.mesh

# Execute check mesh for snappyHexMesh in parallel (mesh partition required)
start_time="$(date -u +%s.%N)"
mpirun -np 8 checkMesh -latestTime -allGeometry -allTopology -parallel | tee log.06.checkMesh.snappy
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute check mesh for snappyHexMesh: $elapsed seconds" >> log.mesh

# Reconstruct snappyHexMesh mesh
start_time="$(date -u +%s.%N)"
reconstructParMesh -latestTime -constant | tee log.07.reconstructParMesh
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Reconstruct snappyHexMesh mesh: $elapsed seconds" >> log.mesh

# Remove mesh partition folders
rm -rf processor* > /dev/null 2>&1

# Renumber mesh to make the sparse matrix more diagonal
start_time="$(date -u +%s.%N)"
renumberMesh -overwrite | tee log.08.renumberMesh
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute solver partition: $elapsed seconds" > log.mesh

# Optional - Create VTK files to check mesh and problem setup
# start_time="$(date -u +%s.%N)"
# foamToVTK -latestTime | tee log.09.foamToVTK
# end_time="$(date -u +%s.%N)"
# elapsed="$(bc <<<"$end_time-$start_time")"
# echo "Create VTK files: $elapsed seconds" >> log.mesh

