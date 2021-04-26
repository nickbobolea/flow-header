#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Copy 0 boundary conditions
cp -r 0.orig 0

# Execute solver partition
start_time="$(date -u +%s.%N)"
decomposePar | tee log.09.decomposeParSolver
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute solver partition: $elapsed seconds" > log.solve

# Execute solver
start_time="$(date -u +%s.%N)"
# Parallel (solver partition required)
mpirun -np 8 buoyantPimpleFoam -parallel | tee log.10.buoyantPimpleFoam
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Execute solver: $elapsed seconds" >> log.solve

# Optional - Solver monitoring
# foamMonitor -l postProcessing/residuals/0/*.dat
# foamMonitor postProcessing/outlet_temperature/0/*.dat
# foamMonitor postProcessing/inlet1_massflow/0/*.dat
# foamMonitor postProcessing/inlet2_massflow/0/*.dat
# foamMonitor postProcessing/outlet_massflow/0/*.dat
# foamMonitor postProcessing/pressure_drop/0/*.dat

# Reconstruct solver solution
start_time="$(date -u +%s.%N)"
reconstructPar | tee log.11.reconstructPar
# reconstructPar -fields '(epsilon k p phi p_rgh T U yPlus)' | tee log.11.reconstructPar
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Reconstruct solver solution: $elapsed seconds" >> log.solve

# Optional - Post-process results
# buoyantPimpleFoam -postProcess -latestTime -dict system/controlDict | tee log.13.postProcess
# mpirun -np 8 buoyantPimpleFoam -parallel -postProcess -latestTime -dict system/controlDict | tee log.13.postProcess

#------------------------------------------------------------------------------
