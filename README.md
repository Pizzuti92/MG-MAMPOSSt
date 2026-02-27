<a name="top"></a>

# MG-MAMPOSSt, a code to test gravity with the mass profiles of galaxy clusters 

- [Go to Main page](https://github.com/Pizzuti92/MG-MAMPOSSt)
- [Introduction](#Introduction)
- [Documentation](#Documentation)
- [Install](#Install)
- [Overview of usage](#Overview-of-usage)
- [Test and Tutorial](#Test-and-Tutorial)
- [Extension and Update](#Extension-and-Update)
## Introduction

MG-MAMPOSSt is a FORTRAN code that extended the MAMPOSSt algorithm of G. Mamon, A. Biviano and G. Boué - 
which performs Bayesian fits of models of mass and velocity anisotropy profiles to the distribution of tracers in projected phase space -
to handle modified gravity models and constrain their parameter space. The new version implements several types of gravity modifications, 
namely general chameleon (including $`f(\mathcal{R})`$ models), beyond Horndeski gravity (Vainshtein screening) Refracted Gravity, Clustering Dark Energy. MG-MAMPOSSt efficently explores the parameter space 
either by computing the likelihood over a multi-dimensional grid of points or by performing a simple Metropolis-Hastings MCMC. When Brightest Cluster Galaxy (BCG) data and gas/galaxies mass profiles are available,
the algorithm can be run in a multi-component fashon to extract the profile of the dark matter in $\Lambda$CDM and non-standard scenarios.
The code requires a Fortran90 compiler or higher and makes use of the Python3 [getdist package](https://github.com/cmbant/getdist) of Antony Lewis to plot the marginalized distributions in the MCMC mode.

### Documentation

MG-MAMPOSSt has been developed by L. Pizzuti, I.D. Saltas G. Mamon L. Amendola and A. Biviano from the MAMPOSSt version of A. Biviano, with a precious contribution by S. Sartor and S. Sanseverinati. 
This is version 2.0, which implements some additional features with respect to version 1.0, including:
- Multi-Component run
- new MG models
- Gaussian priors for MCMC
- possibility to introduce a lensing chain and perform a more robust joint analysis.
- new models for velocity anisotropy
- a new free parameter for the velocity anisotropy (the scale radius $r_\beta$).
- a single input file to work with.


To test some of the subroutines, refer to the source code testMAM.f, which can be compiled and executed as
```bash
 f95 -o testMAM.e testMAM.f -L build/ -lMAM -L build/GamI/ -lGamI -L build/Newuoa/ -lNewuoa -L build/JJin/ -lJJin -L build/Utili/ -lUtili -L build/Powell/ -lPowell
 ./testMAM.e
``` 
Details of the original MAMPOSSt method can be found in [Mamon et al., 2013](https://academic.oup.com/mnras/article/429/4/3079/1012297?login=false). The updated version of the (GR) MAMPOSSt code can be downloaded [here](https://gitlab.com/gmamon/MAMPOSSt). 


---

## Install

To install and run MG-MAMPOSSt, download the .zip file which contains all the dependencies. Note that the main source code gomamposstoptS.f requires additional routines which are stored in the various folders shipped within the code. Some routines are taken from free FORTRAN libraries available on the web. Credits to the developers are given in the header of the source files of these routines.

### Install with CMake

The configuration and installation requires at least CMake 3.17.1. 
In the working directory, execute the following commands:
```bash
$ mkdir build
$ cd build/
$ cmake ..
$ cmake --build .
$ sudo cmake --install .
$ cd ..
```
To install it locally, type ```$ cmake ---DCMAKE_INSTALL_PREFIX:PATH=/Your/path_to/folder/MG_MAMPOSSt/sw ..````

To run and test the code, execute:
```bash
$ gomamposstopt < gomamposst_x.inp
```
which produces the main outputs, or
```bash
$ ./script/script_runmam.sh  
```
which also generates additional plots if the MCMC mode is selected (see below).
Note that to run the above script, permissions should be changed to make it executable. Otherwise, one can simply use the ``` sh ``` environment.

## Overview of usage

Here we summarize all the necessary information to perform a complete run of the code.

### Dataset

The directory data/ stores the datafiles of projected phase spaces (p.p.s) that serve as input to the MG-MAMPOSSt procedure. The files are structured as a table where the number of rows coincides with the number of data points. The first column is the projected radius in units of ```kpc```, the second and thirds columns represent the l.o.s. velocities and the associated errors in units of ```km/s```. Note that data points should be given in the rest frame of the cluster.
The first two lines are considered as comment lines when MG-MAMPOSSt read the data.
If a fourth column is given, the numbers should correspond to the weights associated to each measurement. Set the parameter **weights=1** in the input parameter file (see [next section](#input-parameters)), otherwise the weights will be ignored.

The test data file provided corresponds to a projected phase space of the Coma galaxy cluster obtained using publicy available DESI DR1 data. 
A mock velocity dispersion file for a BCG and a mock lensing chain are further provided for testing the different code options.


# MG-MAMPOSSt

## Input Parameters

Input parameters are stored in:

```
input_pars/pars_test.txt
```

All parameters must follow the format:

```
label = value
```

- **Labels are mandatory**
- **Values** can be omitted only if a default is implemented
- Guess values must always be provided when the parameter is explored

Parameters can be:

- `integer*4`
- `real*8`
- `string`

They are divided into the following groups.

---

# 1. Number of Iterations / Steps (integers)

These control whether a parameter is fixed or explored.

## Grid mode (`nmcmc = 0`)
- Value = number of grid points.
- If an even number is provided, the code may internally use `n+1` to symmetrically sample around the guess value.

## MCMC mode (`nmcmc = 1`)
- If steps ≠ 0 → parameter is free.
- If steps = 0 → parameter is fixed to its guess value.
- Special values (-1, -2) activate specific physical assumptions.

---

## Structural Parameters

| Parameter | Meaning |
|------------|----------|
| `nr200` | Steps for virial radius $r_{200}$ |
| `nrc` | Steps for tracer scale radius $r_c$ |
| `nrs` | Steps for mass scale radius $r_s$ |

Special cases:

```
nrc = -1   → Light Follows Mass (r_c = r_s)
nrc = -2   → Fit N(R) first, then fix in MG-MAMPOSSt

nrs = -1   → Mass Follows Light (r_s = r_c)
nrs = -2   → LCDM concentration relation (e.g. Macciò+2008)
```

---

## Velocity Anisotropy Parameters

| Parameter | Meaning |
|------------|----------|
| `nbeta` | Steps for primary anisotropy parameter |
| `nbeta2` | Steps for second anisotropy parameter (used only for `gT` or `gOM`) |
| `nrbeta` | Steps for anisotropy scale radius `r_beta` |

Special cases:

```
nbeta = -1   → ML profile with r_beta = r_s
nbeta = -2   → Hansen–Moore relation
```

---

## Modified Gravity Parameters

| Parameter | Meaning |
|------------|----------|
| `nA1` | Steps for MG parameter $A_1$ |
| `nA2` | Steps for MG parameter $A_2$ |
| `nA3` | Steps for MG parameter $A_3$ |

Special case:

```
nA2 = -1   → Enforces f(R) sub-case (chameleon)
```

---

## BCG Parameters

| Parameter | Meaning |
|------------|----------|
| `njaf` | Steps for Jaffe radius |
| `nxmas` | Steps for stellar M/L |

---

# 2. Guess Values (reals)

Used as initial guesses for grid or MCMC exploration.

If steps = 0 → guess remains fixed.

---

## Structural Parameters

| Parameter | Units | Meaning |
|------------|-------|----------|
| `r200g` | Mpc | Virial radius |
| `rcg` | Mpc | Tracer scale radius |
| `rsg` | Mpc | Mass scale radius |

---

## Velocity Anisotropy

| Parameter | Meaning |
|------------|----------|
| `betag` | Primary anisotropy parameter |
| `beta2g` | Secondary anisotropy parameter |
| `rbetag` | Anisotropy scale radius (Mpc) |

Note:
- Some models use the scaled parameter $(1 - \beta)^(-1/2)$
- ML/OM-like models use Mpc units for scale radius

---

## Modified Gravity/ additional parameters

| Parameter | Meaning |
|------------|----------|
| `A1g` | MG parameter A1 |
| `A2g` | MG parameter A2 |
| `A3g` | MG parameter A3 |

Units and meanings depend on selected `M(r)` model.
- in General Chameleon: A1 is the background field $\phi_\infty$, A2 is the coupling constant $Q$, A3 is not used.
if the model is ```mtotal_GC```, then A3 is the $\gamma$ inner slope of the gNFW density profile of DM.
- if the model is ```gNFW``` A3 is the $\gamma$ inner slope of the gNFW density profile of DM.
- if the model is Beyond-Horndeski/Vainsthein screening: A1 is the parameter $Y_1$ and A2 is the parameter $Y_2$ (used only if lensing distribution is included!). A3 is not used
- if the model is Refracted gravity, A1 is the permittivity in vacuum $\epsilon_0$, A2 is the logarithm of the effective density in ($M_\odot/Mpc^3$) and A3 is the steepness of the transition $Q$.
- in Clustering DE, A1 is the speed of sound, A2 the ratio between DE and DM perturbation, and A3 is the exponent of the gNFW model.



---

## BCG

| Parameter | Units | Meaning |
|------------|-------|----------|
| `rjafg` | Mpc | Jaffe radius |
| `xmstarg` | M☉/L☉ | Stellar M/L |

---

# 3. Model Options (reals / integers / strings)

---

## Cosmology

| Parameter | Default | Description |
|------------|----------|-------------|
| `H0` | 70 | km/s/Mpc |
| `za` | 0 | Cluster redshift |
| `Olam` | 0.7 | ΩΛ |
| `Omegam` | 0.3 | Ωm |

---

## Radial Cuts

| Parameter | Meaning |
|------------|----------|
| `Rlow` | Inner projected radius (Mpc) |
| `Rup` | Outer projected radius (Mpc) |

---

## Tracer Number Density Profile

| Parameter | Options |
|------------|----------|
| `N(R)` | `pNFW`, `pHer`, `beta` |
| `al` | Exponent (only if `beta` model is selected) |

---

## Mass / Potential Model (`M(r)`)

### GR models

```
NFW, Her, PIEMD, Bur, SoftIS, Eis5, gNFW
```

### Modified gravity models

```
mNFW_BH, mNFW_GC, mBur_GC,
mgNFW_GC, mIso_GC, mTotal_GC,
RG, BS
```

Additional parameters:

| Parameter | Meaning |
|------------|----------|
| `rcut` | Truncation radius (used only for PIEMD) |
| `eim` | Exponent of Einasto model (if selected), exponent of gNFW and Hernquist if mgNFW_GC or mNFW_GC are selected |

---

## Velocity Anisotropy Model (`Beta(r)`)

Available models:

```
C, ML, OM, WJ, T, O, gOM, gT
```

Additional option:

| Parameter | Meaning |
|------------|----------|
| `FreeBeta` | 0/1 — allow anisotropy scale radius to vary |

---

## Performance / Optimization

| Parameter | Meaning |
|------------|----------|
| `FASTMODE` | 0/1 — likelihood interpolation mode |
| `OPT` | Optimization algorithm: -1 (none), 0 (BOBYQA), 1 (NEWUOA), 2 (POWELL) |


---

## BCG and Baryons

| Parameter | Meaning |
|------------|----------|
| `BCG` | Include Jaffe BCG |
| `BCG_Data` | Include external BCG VDP data |
| `Baryons` | Include gas + galaxies as fixed components |

Gas parameters:

```
rho0gas
rsgas
srhoga
srsga
```

Galaxy parameters:

```
rhostar
rstar
srhostar
srstar
```

---

# 4. Parameter Limits (reals)

Used when parameter limits are enabled in the run options.

| Parameter | Meaning |
|------------|----------|
| `r2low`, `r2up` | Bounds for r200 |
| `rclow`, `rcup` | Bounds for r_c |
| `rslow`, `rsup` | Bounds for r_s |
| `blow`, `bup` | Bounds for primary anisotropy |
| `b2low`, `b2up` | Bounds for secondary anisotropy |
| `A1low`, `A1up` | Bounds for A1 |
| `A2low`, `A2up` | Bounds for A2 |
| `A3low`, `A3up` | Bounds for A3 |
| `rbetalow`, `rbetaup` | Bounds for anisotropy scale radius |
| `rjalow`, `rjaup` | Bounds for Jaffe radius |
| `xmalow`, `xmaup` | Bounds for BCG stellar M/L |

---

# 5. Exploration Switches

| Parameter | Meaning |
|------------|----------|
| `nmcmc` | 0 = grid, 1 = MCMC |
| `Nsample` | Number of MCMC samples |
| `GaussianPrior` | 0 flat prior, 1 Gaussian prior |
| `FlatPriorBeta` | Keep flat prior for β |
| `FlatPriorBeta0` | Keep flat prior for β₀ |

---

# 6. Lensing Options (optional)

| Parameter | Meaning |
|------------|----------|
| `nlens` | Include lensing likelihood |
| `r200t`, `rst` | Peak values for lensing priors |
| `delta1`, `delta2`, `delta3` | Lensing hyperparameters |
| `lensing_file` | Include lensing chains from file |
| `nLogP` | Explore $\phi$ in log-space (chameleon mode) |

---

For a complete example configuration, see:

```
input_pars/pars_test.txt
```


### Working Options

The file Options.txt is deprecated.


### Output
Output files can be found in the Output folder. The names and the location of those files can be changed by modyfing the file **gomamposst_x.inp**. The main output file is **MaxLik.dat**, 
which stores the tabulated ln(L) where L is the  likelihood/posterior, as a function of the model parameters. The other accessory files are described in the documentation.


## Test and Tutorial

### Basic smoke test
In order to check that the installation process was successfull and to produce the test example, type:
```bash
$ ./script/script_runmam.sh  
```
without changing any parameter in the input files. This way, the code should execute a 100000-points MCMC exploration of the parameter space in Vainsthein screening gravity  (**M(r)=mNFW_BH**) sampling the full kinematics+lensing likelihood by using the test data-set included in the data/ folder. The script runs in Fast Mode **FASTMODE=1** and should provide the complete output posterior in Output/MaxLik.dat within less than half an hour if the execution is performed over an average laptop. 
At the end of the run, if successufull the following text will be displayed:
```bash
 Best-fit from optimization
   r_200    =  1.406
   r_tracer =  0.330 (c_tracer=    4.26)
   r_mass   =  0.369 (c_mass=      3.81)
   Anisotropy parameter =   1.5691
   First MG parameter =  0.0828
   Second MG parameter =  0.0100
   Second Anisotropy parameter =  1.0000
 Likelihood =   5911.38999
 
   After MG-MAMPOSSt:

     build output files of binned N(R), VDP
     of best-fit N(R), VDP and of
     input N(R), VDP solutions ('true' values)

  Binned N(R) computed
  Using           13  bins for the VDP(R)
  Evaluating expected VDP for
  Max Lik solution values:
   r_200    =  1.390
   r_tracer =  0.330 (c_tracer=    4.21)
   r_mass   =  0.346 (c_mass=      4.01)
   Anisotropy parameter =   1.5457
   Second anisotropy parameter =   1.0000
   First MG parameter =  -0.0991
   Second MG parameter =  -0.0725

  sigma_r evaluated
 ```
 The plot of the marginalized distributions of the free parameters generated by the test run should correspond to what is shown by [this figure](https://github.com/Pizzuti92/MG-MAMPOSSt/blob/main/test/plot_example_test1.png). Note that the final best fit(s) can slightly chnages depending on the points sampled in the MCMC/grid.
 
### A simple tutorial of usage
In the following, we present a simple tutorial to guide the user over a step-by-step execution of the code. 
#### File name selection
First, open the gomamposst_x.inp file, which should appear as:
```bash
data/datphys_test.dat
input_pars/pars_test.txt
Output/rnvn.dat
Output/NRpar.dat
Output/NRbin.dat
Output/NRfit.dat
Output/MaxLik.dat
Output/svbin.dat
Output/svfit.dat
```
The first line is the input data-set - Be sure that the format of the data-set agrees with the documentation (the first two lines are considered as comments) - while the second line corresponds to the file storing the input parameters. All the other lines indicate the output file names. In this tutorial we will perform a test run in Chameleon Screening (**M(r)=mNFW_GC**); as such, we change the output likelihood file as
```bash
Output/MaxLik_ChameleonTest.dat
```
#### Parameter definition
Now, we define the details of the run by working on the parameters in pars_test.txt. In particular, we choose to work with the two free parameters defining the modified gravity model, which are the background value of the chameleon field **\phi_\infty**, given in unit of 1e-5 and the coupling constant **Q**. This means that the only numbers different from zero in the first block should be the fifth and the sixth:
```bash
**************  number of  steps in free parameters ********************

nr200 = 0     
               ! number of steps for r200 fit

nrc = 0        
               ! number of steps for rc fit, scale radius of N(R)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces LfM, c_nu=c]
                            !   [if = -2 fit N(R) outside MAMPOSSt]
nrs =  0       
			   ! number of steps for rs fit, scale radius of M(r)
                            !   [if = 0 takes guess value]
                            !   [if = -1 forces MfL, c=c_nu]
                            !   [if = -2 forces LCDM, c=c(M)]

nbeta = 0     
			   ! number of steps for anisotropy parameter
                            !   [if = -1 forces a_ML=r_s]
                            !   [if = -2 forces Hansen+Moore]
nA1 = 50 
			   ! number of steps in the first MG parameter 

nA2 = 30        
			   ! number of steps in the second MG parameter
			   ! If equal to -1 force the case of chameleon f(R) gravity
			   
nbeta2 = 0     ! number of steps for the second anisotropy parameter                             
************************************************************************
```
Note that, in the case of a fixed-grid exploration, this corresponds to a 51 X 31 points grid. As for the second group, we select the following guess values of the input parameters:
```bash
************* free parameters initial guess values *********************

r200g = 1.41   
			   ! mass profile r200 initial guess (Mpc)

rcg = 0.33    
			   ! N(R) scale radius initial guess (Mpc)

rsg = 0.33    
			   ! mass profile scale radius initial guess (Mpc)

betag = 1.41   
			   ! Anisotropy initial guess, beta_C, a_ML, a_OM, a_W, beta_inf 

A1g = 1.2      
			   ! first MG parameter initial guess

A2g = 0.4      
			   ! second MG parameter initial guess

beta2g = 1.0   ! Second Anisotropy parameter initial guess beta0 for 
               ! gOM and gT models
************************************************************************
```
Except for the modified gravity parameters, the others are fixed to the true values of the halo from  which the test data-set is generated. Now, let's give a look to the third group of parameters:
```bash
*************  Model options *******************************************

H0 = 70.0      
               ! Hubble constant at z=0              

za = 0.0      
               ! average redshift of the cluster (needed to evaluate Hz)                                        

Olam = 0.7     
               ! Omega lambda density parameter

Omegam = 0.3   
               ! Omega matter density parameter

Rlow = 0.05    
               ! Inner projected radius for sample selection (Mpc)

Rup = 1.41     
			   ! Outer projected radius for sample selection (Mpc)


N(R) = pNFW
               ! model for number density profile N(R)
               ! Here we list the labels for N(R) and then the corresponding profile:
               ! pNFW=projected NFW / pHer=projected Hernquist / beta=beta-model (1/2/3)
    
al = -0.00
               ! N(R) negative exponent (only if knfit=3)

M(r) = mNFW_GC	  
              ! mass profile model. 
              ! Here we list the labels for M(r) and then the corresponding profile:
              !  
              ! NFW=NFW/ Her=Hernquist/ PIEMD=PIEMD/ Bur=Burkert/ SoftIS=SoftIS/
              ! Eis5=Einasto_m=5/ mNFW_LH = mod_NFW linear Horndeski/
              ! mNFW_BH = mod_NFW beyond Horndeski/mNFW_GC= mod_NFW general chameleon/
              ! gNFW = generalized NFW
              

Beta(r) = T   
              ! Anisotropy model.
              ! Here we list the labels for beta(r) and then the corresponding profile: 
              ! C=Constant/  ML=Mamon&Lokas/ OM=Osipkov-Merritt/  WJ=simplified Wojtak/
              ! T=simplified Tiret/ O= (Opposite) modified Tiret 
              
rcut= 1.0     
              ! PIEMD model rcut in Mpc (only used if M(r)=PIEMD) 

FASTMODE= 0  
              ! Run MG-MAMPOSSt in fast mode interpolating over data points

OPT= 2        
              ! optimization algorithm: 0/1/2=bobyqa/newuao/powell
              ! -1 skip optimization
              
screen= 0   
              !     for linear f(R) kmp.eq.7, one can decide to set an instantaneous
	      !     transition between screeening and linear regime, by using the 
	      !     analytical approximation of Lombriser+12.                
              !-1/0/1/2=noscreen (general Hordenski)/noscreen f(R)/screen(instantaneous transition)
              !/screen (arctan transition)                    
			 
************************************************************************
```
The first four values represent the background cosmology (note that the test phase space has been generated at redshift 0). Then, we set the range of tracers to be considered in the MG-MAMPOSSt analysis [0.05 Mpc,1.41 Mpc] and we choose the  model **pNFW** for the number density profile of the tracers. The next parameter **-0.0** is ignored in this case.

The effective gravitational potential (mass profile) in Chameleon gravity is selected by the number **9** and the model of the velocity anisotrpy profile is chosen to be a Tiret. The following parameter 1.0 is not considered unless MG-MAMPOSSt runs with the model **PIEMD** for the mass profile. 
If we want to execute the code in fast mode, set the next paramater equal to 1, otherwise (as this is the case), type '0'. 
Finally, we choose the optimization algorithm **2** (Powell); the last parameter is ignored unless the mass profile model is **mNFW_LH**.

The last step is to choose the parameter limits for your exploration (only if **kpro=1**, see below).
```bash
********* parameter limits *********************************************

r2low= 0.4    
              !r200 lower bound

r2up= 5.0     
              !r200 upper bound

rclow= 0.05   
              !rc lower bound

rcup= 3.9     
              !rc upper bound

rslow= 0.04   
              !rs lower bound


rsup= 3.9     
              !rs upper bound

blow= 0.5    
              !beta lower bound

bup= 7.1      
              !beta upper bound

A1low = 0.006  
              !first MG parameter lower bound

A1up =   100.0
			  !first MG parameter upper bound

A2low=  0.02 
			  !second MG parameter lower bound

A2up= 100.0    
			  !second MG parameter upper bound

b2low = 0.5
			  !second beta parameter lower bound

b2up = 7.1 
			  !second beta parameter upper bound


************************************************************************
```

#### Select running options
Close the "pars_test.txt" file and open "Options.txt" to customize the execution of MG-MAMPOSSt. To choose the fixed-grid parameter space exploration, type 
```bash
nmcmc=0
```
In this case, the following parameter **Nsample** is ignored. We further exclude the additional lensing distribution by setting
```bash
nlens=0
```
To select the customized limits of the parameter space (the last group of inputs in "Pars_all_N_O_spec_DS"), set 
```bash
kpro=1
```
We leave the rest of the file unchanged. 

#### Execution
In the parent working directory, execute
```bash
$ ./script/script_runmam.sh -t  
```
or
```bash
$ sh script/script_runmam.sh -t  
```
where the additional option ```-t``` prints the execution time. If everything worked fine, at the end of the run one should obtain the following printed message:
```bash 
 Best-fit from optimization (kinematic only)
   r_200    =  1.410
   r_tracer =  0.330 (c_tracer=    4.27)
   r_mass   =  0.330 (c_mass=      4.27)
   Anisotropy parameter =   1.4100
   First MG parameter =  0.0087
   Second MG parameter =  0.0386
   Second Anisotropy parameter =  1.0000
 Likelihood =   5908.30636


  After MG-MAMPOSSt:

     build output files of binned N(R), VDP
     of best-fit N(R), VDP and of
     input N(R), VDP solutions ('true' values)

  Binned N(R) computed
  Using           13  bins for the VDP(R)
  Evaluating expected VDP for
  Max Lik solution values:
   r_200    =  1.410
   r_tracer =  0.330 (c_tracer=    4.27)
   r_mass   =  0.330 (c_mass=      4.27)
   Anisotropy parameter =   1.4100
   Second anisotropy parameter =   1.0000
   First MG parameter =   0.1197
   Second MG parameter =   0.0249

  sigma_r evaluated
******************************
time of execution:  12  seconds
******************************
No plot for the grid case
```
Note that this time of execution refers to a laptop ASUS Intel(R) Core(TM) i7-8565U CPU 1.99 GHz. The output "MaxLik_ChameleonTest.dat" file should be the same found in the "test" folder in this repository. 

In order to perform the same run in MCMC mode type ```nmcmc=1``` in "Options.txt". We select ```Nsample=100000``` points in the chain.

Running ```$ sh script/script_runmam.sh -t``` in this case should produce the [marginalized distributions](https://github.com/Pizzuti92/MG-MAMPOSSt/blob/main/test/plot_test_MCMC_Chameleon.png) in the "test" folder. The corresponding (log) likelihood is stored in "MaxLike_ChameleonTestMCMC.dat" in the same folder; the execution on a laptop with the features specified above should take roughly 10 minutes.   


## Extension and update

Since the expression of the gravitational potential is involved only in the solution of the Jeans' equation, new mass models or modified gravity/dark energy parametrizations can be introduced by modifying the routines where the potential enters, which are the functions fa(tlog) and sr2int(alr) in MAM.f. Add a new entry in the condition chain, identified by a new integer number **kmp** which should be string for **M(r)**. 






Go to [Top Page](#top)

