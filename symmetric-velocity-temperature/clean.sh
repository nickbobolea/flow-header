#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

# Remove "0" folder (ICs and BCs)
rm -rf 0 > /dev/null 2>&1

# Clean "constant" folder
# Remove blockMesh files
rm -rf constant/polyMesh > /dev/null 2>&1
# Remove edge extraction files
rm -f constant/triSurface/*.eMesh > /dev/null 2>&1
rm -rf constant/extendedFeatureEdgeMesh > /dev/null 2>&1

# Remove mesh / solver execution folders
# Remove partition folders
rm -rf processor* > /dev/null 2>&1
# Remove solver folders
rm -rf 0 0.[0-9]* [1-9]*
# Remove postProcessing
rm -rf postProcessing > /dev/null 2>&1
# Remove log files
rm -f log.* > /dev/null 2>&1

#------------------------------------------------------------------------------
