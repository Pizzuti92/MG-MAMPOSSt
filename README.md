<a name="top"></a>

# MG-MAMPOSSt, a code to test gravity with the mass profiles of galaxy clusters 

[(Main page)](https://github.com/Pizzuti92/MG-MAMPOSSt)

## Introduction

MG-MAMPOSSt is a FORTRAN code that extended the MAMPOSSt algorithm of G. Mamon, A. Biviano and G. Boué - 
which performs Bayesian fits of models of mass and velocity anisotropy profiles to the distribution of tracers in projected phase space -
to handle modified gravity models and constrain their parameter space. The new version implements two distinct types of gravity modifications, 
namely general chameleon (including $f(\mathcal{R})$ models) and beyond Horndeski gravity (Vainshtein screening). MG-MAMPOSSt efficently explores the parameter space either by computing the likelihood over a multi-dimensional grid of points or by performing a simple Metropolis-Hastings MCMC. 
The code requires a Fortran90 compiler or higher and makes use of the Python3 [getdist package](https://github.com/cmbant/getdist) of Antony Lewis to plot the marginalized distributions in the MCMC mode.

### Documentation

A full description of the code functionalities, input parameters and output files is given in the pdf user manual available [here](https://arxiv.org/abs/2201.07194).


---

## Install

To install and run MG-MAMPOSSt on a Linux terminal, download the .zip file which contains all the dependencies. Note that the main source code gomamposstoptS.f requires additional routines which are stored in the various folders shipped within the code. Some routines are taken from free FORTRAN libraries available on the web. Credits to the developers are given in the header of the source files of these routines.

### Install with CMake

The configuration and installation requires at least CMake 3.17.1. 
Execute the following commands:
```bash
$ mkdir build
$ cd build/
$ cmake ..
$ cmake --build .
$ sudo cmake --install .
$ cd ..
```
To run and test the code, execute:
```bash
$ gomamposst < gomamposst_x.inp
```
which produces the main outputs, or
```bash
$ ./script/script_runmam.sh  
```
which also generates additional plots if the MCMC mode is selected (see below).
Note that to run the above script, permissions should be changed to make it executable. Otherwise, one can simply use the ``` sh ``` environment.


### Input parameters
Input parameters are stored in the file Input_pars/pars_all_N_O_spec_DS.
Additional options are in Options.txt where one can customized the mode 
All the required information about the structure of these files can be found in the [pdf user manual](https://arxiv.org/abs/2201.07194)


### Output
Output files can be found in the Output folder. Main is MaxLik_test.dat, 
which stores the tabulated $-\ln\mathcal{L}$, where $mathcal{L}$ is the 
likelihood/posterior, as a function of the model parameters.



Go to [Top](#top)

