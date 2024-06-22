# NACASMesher

NACASMesher is a fork of the project curiosityFluidsAirfoilMesher.

<img src="https://github.com/LHopitalAstralis/NACASMesher/blob/master/Grid%20Example.png?raw=true" width="500">

## Description

The project consists in a script for generating blockMeshDict files, for OpenFOAM, of NACA 4-digit series airfoils, based on the required parameters.

This project adds, internally to the code, the generation of points for a given selected NACA 4-digit series airfoil.

This generation of points avoids the need, present in the original project, for the importing of external .dat files.

The generation of the mesh is modified from the one in the curiosityFluidAirfoilMesher to a generation based on the project by Thien Phan. 
Available at: https://www.phanquocthien.org/mesh-geometry/blockmesh/airfoil 



