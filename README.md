# HausdorffIsingModel
#### A simulation of the Ising model as a function of the lattice Hausdorff dimension.
------

The Ising model is a simple example of a physical system displaying dimensionality- dependent phase transitions. However, the notion of dimensionality considered in the context of this model tends to follow the Lebesgue picture, as opposed to the more general Hausdorff picture. In particular, the Ising model displays nontrivial dimensionality dependence on fractal lattices. For instance, the nearest-neighbor model on a Sierpinski gasket lattice features no critical behavior above zero temperature, while that on a Sierpinski carpet lattice does. This indicates that there may be a fractional lower critical dimension for the Ising model on Cantor set-like lattices.

Ultimately, the point of intrigue in a critical dimension is to understand: 

1. What constraints limit the physical behavior of a system when it is embedded in a given dimension, and 
2. Whether there is a natural or "equilibrium" interpretation to dimensionality which explains why we can describe real physical phenomena with simple models. 
  
At a large scale, for instance, we might address why and how an equilibrium description of the universe could have three spacial dimensions.

#####Contents Summary:

- `condor`: contains templates for job submission to the USCMS T2 at FNAL 
- `scripts`: `.py`, `.sh` and other useful scripts for various utilities
- `src`: C++ source code for the analysis
- `steerScript.sh`: launch this script to see the steps in the analysis

#####Example lattice plots:

![1D Ising Model](https://github.com/eacoleman/HausdorffIsingModel/blob/master/doc/hausdorff-1.png "1D Ising Model")
![2D Ising Model](https://github.com/eacoleman/HausdorffIsingModel/blob/master/doc/hausdorff-2.png "2D Ising Model")
![3D Ising Model](https://github.com/eacoleman/HausdorffIsingModel/blob/master/doc/hausdorff-3.png "3D Ising Model")
![4D Ising Model](https://github.com/eacoleman/HausdorffIsingModel/blob/master/doc/hausdorff-4.png "4D Ising Model")
