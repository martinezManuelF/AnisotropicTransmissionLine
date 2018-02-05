# AnisotropicTransmissionLine
A MATLAB script that simulates transmission lines embedded in anisotropic media.
Makes use of Yee Staggered grid schemes to compute derivatives using the finite difference method.

Functions included:
* yeeder.m: Derivative operator calculators on a Yee staggered grid.
  -Supports any Nx-by-Ny grid.
  -Supports Dirichlet and Periodic Boundary Conditions.
* RotMat.m: Rotation matrix calculations for 3-by-3 tensors
  -Supports rotations along x,y, and z axes.
* AnisoTL:  Anisotropic transmission line calculator.
  -Supports arbitrary number of signal lines for symmetric or asymmetric transmission lines.
  
Simulation files included:
* AnisoMicroStrip.m:  Builds and simulates a micro strip transmission line embedded in anisotropic media.
* Coplanar.m:         Builds and simulates a coplanar transmission line.
* MicroStrip.m:       Builds and simulates a micro strip transmission line.
* ParallelPlate.m:    Builds and simulates a parallel plate transmission line.
* RG59Coax.m:         Builds and simulates an RG-59 Coaxial transmission line.

Sample simulations:
![RG-59 Coaxial Cable](https://github.com/TasartirAmras/AnisotropicTransmissionLine/blob/master/Graphics/CoaxE.png)

Created as a part of EE 5322--21st Century Electromagnetics at the University of Texas at El Paso.
