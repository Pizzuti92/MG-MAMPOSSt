<a name="top"></a>

# MG-MAMPOSSt 

[(Main page)](https://github.com/Pizzuti92/MG-MAMPOSSt)

### MG-MAMPOSSt, a code to test gravity with the mass profiles of galaxy clusters

MG-MAMPOSSt is a FORTRAN code that extended the MAMPOSSt algorithm of G. Mamon, A. Biviano and G. Boué, 
that performs Bayesian fits of models of mass and velocity anisotropy profiles to the distribution of tracers in projected phase space, 
to handle modified gravity models and constrain its parameters. The new version implements two distinct types of gravity modifications, 
namely general chameleon (including $f(\mathcal{R})$ models) and beyond Horndeski gravity (Vainshtein screening).


---

### Install

To install and run MG-MAMPOSSt on your Linux terminal, download the .zip file which contains all the dependencies. Notice that the main source code gomamposstoptS.f requires additional routines which are stored in the various folders shipped within the code. Some routines are taken from free libraries available on the web. Credits to the developers are given as comments in the source files of these routines.

Execute the following scripts in this order:
```bash
./script/script_Lib.sh 
./script/script_compile.sh
./script/script_runmam.sh  
```
The first script creates the libraries and should be run only once. The second script compiles the main code while the third execute the MG-MAMPOSSt method further producing the additional outputs.

Note that you should change the permissions of the scripts to make them executable. Otherwise, one can simply use the ```bash sh ``` environment.

Type -h for help.

Default compiler is "f95", default directory is "$PWD". If you want to change compiler or set your own path, use the -f and -d options:
```bash
./script/script_Lib.sh -f <compiler_name> -d <directory_path> 
./script/script_compile.sh -f <compiler_name> 
./script/script_runmam.sh  
```

be sure to put a backslash at the end of the path, e.g.:
```bash
./script/script_Lib.sh -f <compiler_name> -d /home/pizzuti/MG-MAMPOSSt/
```

### Input parameters
Input parameters are stored in the file Input_pars/pars_all_N_O_spec_DS.
Additional options are in Options.txt.
All the required information about the structure of these files can be found in the pdf user manual available  
[(here)](https://arxiv.org/abs/2201.07194)

### Output
Output files can be found in the Output folder. Main is MaxLik_test.dat, 
which stores the tabulated $-\ln\mathcal{L}$, where $mathcal{L}$ is the 
likelihood/posterior, as a function of the model parameters.



Go to [Top](#top)

