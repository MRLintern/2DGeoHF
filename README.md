## 2D Heat Diffusion Through Various Geological Materials

### TODO

* Basic one has been provided.
* Theres issues with new lines/indentation that need resolving.
* Issues with the _Configuration Class_; to be resolved.
* Separate `modules/classes` into different files; not a priority to make the software run/work. Just more OOD principles.
* A script for plotting data/results; Probablt `Python based` using `Matplotlib`. 

### Introduction

* This `Fortran Solver` models `Heat Diffusion` through a `2D Domain`, which models 3 different `Geological Materials`.

* The only physical property we need is the ___Coefficient of Thermal Diffusivity___. 

* The rest of the work is purely Mathematical and Computational. 

### Geological Materials

* Thermal Diffusivities for the different rock types can be found in attached `PDF` file called ___Thermal Properties of Rocks___. 

* We consider materials from the: `Sedimentary`, `Metamorphic` and `Igneous Rock Groups`.

#### Sedimentary Rock

* _Sandstone_. Coefficient of Thermal Diffusivity: __1.3e-4 (m*m/s)__.

#### Metamorphic Rock

* _Quartzite_. Coefficient of Thermal Diffusivity: __2.6e-4 (m*m/s)__. 

#### Igneous Rock

* _Basalt_. Coefficient of Thermal Diffusivity: __9e-5 (m*m/s)__.


### Mathematical Model

* `Finite Discretisation` is used to `discretise` the domain. The Solver uses `Space & Time Marching` to solve the `2D Heat-Diffusion Equation`.
* `OpenMP` has been used to `parallelise` the `Marching-Based Numerical Method`.

### Requirements

* Developed in `Linux Ubuntu 20.04`.
* `gfortran compiler`. Version `9.4.0` was used.
* A text editor for making changes.
* `OpenMP` for `Multithreading`.
* `GNU Make`
* This isn't a tutorial. There's an expectation that you know some Physics and Appled Mathematics; e.g. `Thermodynamics`, `Numerical Linear Algebra` and `Numerical Partial Differential Equations`. 

### Getting the Software

* `$ git clone https://github.com/MRLintern/2DGeoHF.git`
* `$ cd 2DGeoHF` for the purposes of building and running the software.

### Resources

* [Temperature dependence of the thermal diffusivity of sandstone](https://www.sciencedirect.com/science/article/pii/S0920410516312712#sec)
* [A 2-D Barakat-Clark finite difference numerical method and its comparison for water quality simulation modeling](https://environmentalsystemsresearch.springeropen.com/articles/10.1186/2193-2697-2-11)


