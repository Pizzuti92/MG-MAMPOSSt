<a name="top"></a>

# MG-MAMPOSSt, a code to test gravity with the mass profiles of galaxy clusters 

[(Main page)](https://github.com/Pizzuti92/MG-MAMPOSSt)

## Introduction

MG-MAMPOSSt is a FORTRAN code that extended the MAMPOSSt algorithm of G. Mamon, A. Biviano and G. Boué - 
which performs Bayesian fits of models of mass and velocity anisotropy profiles to the distribution of tracers in projected phase space -
to handle modified gravity models and constrain their parameter space. The new version implements two distinct types of gravity modifications, 
namely general chameleon (including $f(\mathcal{R})$ models) and beyond Horndeski gravity (Vainshtein screening). MG-MAMPOSSt efficently explores the parameter space either by computing the likelihood over a multi-dimensional grid of points or by performing a simple Metropolis-Hastings MCMC. 
The code requires a Fortran90 compiler or higher and makes use of the Python3 [@getdist package](https://github.com/cmbant/getdist) of Antony Lewis to plot the marginalized distributions in the MCMC mode.

### Documentation

A full description of the code functionalities, input parameters and output files is given in [the user manual](https://arxiv.org/abs/2201.07194)


---

## Install

To install and run MG-MAMPOSSt on a Linux terminal, download the .zip file which contains all the dependencies. Notice that the main source code gomamposstoptS.f requires additional routines which are stored in the various folders shipped within the code. Some routines are taken from free libraries available on the web. Credits to the developers are given as comments in the source files of these routines.

### Install with CMake

The configuration and installation requires at least CMake 3.17.1.  Unzip the MG-MAMPOSSt folder and navigate to the working directory, e.g.
```bash
$ unzip MG-MAMPOSSt-main.zip
$ cd MG-MAMPOSSt-main/
```
Proceed to configure and install the code through the following commands:
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
which only produces the main outputs, or
```bash
$ ./script/script_runmam.sh  
```
which also generates additional plots if the MCMC mode is selected (see below).
Note that to run the above script, permissions should be changed to make it executable. Otherwise, one can simply use the ```bash sh ``` environment.

### Install and run with bash scripts

An alternative installation can be made by using bash script to generate the libraries and compile the source code.
Execute the following scripts in this order:
```bash
$ ./script/script_Lib.sh 
$ ./script/script_compile.sh
$ ./script/script_runmam.sh  
```
The first script creates the libraries and should be run only once. The second script compiles the main code while the third execute the MG-MAMPOSSt method further producing the additional outputs.

As before, one should change either the permissions of the scripts or run the ```bash sh ``` command.

Type -h for help.

Default compiler is "f95", default directory is "$PWD". If you want to change compiler or set your own path, use the -f and -d options:
```bash
$ ./script/script_Lib.sh -f <compiler_name> -d <directory_path> 
$ ./script/script_compile.sh -f <compiler_name> 
$ ./script/script_runmam.sh  
```

be sure to put a backslash at the end of the path, e.g.:
```bash
./script/script_Lib.sh -f <compiler_name> -d /home/pizzuti/MG-MAMPOSSt/
```

### Input parameters
Input parameters are stored in the file Input_pars/pars_all_N_O_spec_DS.
Additional options are in Options.txt where one can customized the mode 
All the required information about the structure of these files can be found in the [pdf user manual](https://arxiv.org/abs/2201.07194)


### Output
Output files can be found in the Output folder. Main is MaxLik_test.dat, 
which stores the tabulated $-\ln\mathcal{L}$, where $mathcal{L}$ is the 
likelihood/posterior, as a function of the model parameters.



Go to [Top](#top)

