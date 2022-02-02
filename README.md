<a name="top"></a>

# MG-MAMPOSSt, a code to test gravity with the mass profiles of galaxy clusters 

- [Go to Main page](https://github.com/Pizzuti92/MG-MAMPOSSt)
- [Introduction](#Introduction)
- [Documentation](#Documentation)
- [Install](#Install)
- [Overview of usage](#Overview-of-usage)

## Introduction

MG-MAMPOSSt is a FORTRAN code that extended the MAMPOSSt algorithm of G. Mamon, A. Biviano and G. Boué - 
which performs Bayesian fits of models of mass and velocity anisotropy profiles to the distribution of tracers in projected phase space -
to handle modified gravity models and constrain their parameter space. The new version implements two distinct types of gravity modifications, 
namely general chameleon (including $`f(\mathcal{R})`$ models) and beyond Horndeski gravity (Vainshtein screening). MG-MAMPOSSt efficently explores the parameter space either by computing the likelihood over a multi-dimensional grid of points or by performing a simple Metropolis-Hastings MCMC. 
The code requires a Fortran90 compiler or higher and makes use of the Python3 [getdist package](https://github.com/cmbant/getdist) of Antony Lewis to plot the marginalized distributions in the MCMC mode.

### Documentation

MG-MAMPOSSt has been developed by L. Pizzuti, I.D. Saltas G. Mamon L. Amendola and A. Biviano from the MAMPOSSt version of A. Biviano. A precious contribution has been given by J. Dvoracek and S. Sartor. A full description of the code functionalities, input parameters and output files is given in the pdf user manual available [here](https://arxiv.org/abs/2201.07194). Details of the original MAMPOSSt method can be found in [Mamon et al., 2013](https://academic.oup.com/mnras/article/429/4/3079/1012297?login=false). The updated version of the (GR) MAMPOSSt code can be downloaded [here](https://gitlab.com/gmamon/MAMPOSSt) 


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

## Overview of usage

Here we summarize all the necessary information to perform a complete run of the code.

### Dataset

The directory data/ stores the datafiles of projected phase spaces (p.p.s) that serve as input to the MG-MAMPOSSt procedure. The files are structured as a table where the number of rows coincides with the number of data points. The first column is the projected radius in units of ```kpc```, the second and thirds columns represent the l.o.s. velocities and the associated errors in units of ```km/s```. Note that data points should be given in the rest frame of the cluster.
The first two lines are considered as comment lines when MG-MAMPOSSt read the data.

### Input parameters
Input parameters are stored in the file Input_pars/pars_all_N_O_spec_DS. Different parameters go on different rows and must be written starting by the first column. They are all mandatory and are divided in four main groups:

* **Number of iterations (lines 1-6, integers):** for the grid search mode, they represent the number of points in each free parameter over which the likelihood is computed. If set to 0 or 1, the parameter is fixed to its guess value, except for specific cases. When MG-MAMPOSSt is in MCMC mode, if **number of iterations** is different from zero, then the corresponding parameters are optimized within the chain. Otherwise the parameters are fixed to the guess value. 
   - **nr200**: Number of iteration for the virial radius.

   - **nrnu**: Number of iteration for the tracers scale radius r_nu. For **nrnu=-1**  assumes that r_nu is equal to the scale radius of the mass profile r_s in the fit (option: "Light Follows Mass"). If **nrnu=-2** the code first fits the number density profile alone and the best fit found is then used as fixed value in the MG-MAMPOSSt procedure. 
   - **nrs**: Number of iteration for the mass profile scale radius.  For **nrs-1**  assumes that r_s is equal to the scale radius of the number density profile r_nu in the fit (option: "Mass Follows Light"). If **nrs=-2** the mass scale radius is computed by using the theoretical relation of Macciò et al., 2008 (LambdaCDM option).
   - **nbeta**: Number of iteration for the anisotropy parameter. If **nbeta=-1** the anisotropy profile is forced to be a Mamon&Lokas profile with beta=r_beta=r_s. If **nbeta=-2** the Hansen-Moore model is assumed (beta(r) related to the matter density rho_m(r)).
   - **nA1**: Number of iteration for the first MG parameter.
   - **nA2**: Number of iteration for the second MG parameter. For the case of general chameleon gravity, A_2 corresponds to the coupling constant **Q**. If **nA2=-1** it forces the case of $f(R)$ gravity (Q=1/\sqrt{6}).
   
* **Guess values (lines 8-13, reals):** they serve as an initial guess for the MG-MAMPOSSt fit both in grid and MCMC mode. If the corresponding number of iterations is set to zero, or one when in the grid-search option, the parameter guess is kept fixed within the code.
   * **r_200**: guess starting value of the characteristic "virial" radius for the mass profile, measured in units of Mpc. 

    * **r_nu**: guess starting vaule for the scale radius of the number density profile,  units of Mpc. 
    
    * **r_s**: guess starting value for the scale radius of the mass profile, in  units of Mpc. 

    * **beta**: starting guess for the velocity anisotropy profile parameter. If the selected model is "constant", "Tiret" or "Modified Tiret", the parameter is dimensionless, otherwise it should be given in units of Mpc.
    
    * **A_1**: Initial guess for the first modified gravity parameter. For the case of gNFW model **kmp=10**, **A_1=gamma** is the free exponent which characterizes the profile.
    
    * **A_2**: Initial guess for the second modified gravity parameter.
  
* **Model options (lines 15-27, reals/integers):** this family of values allow to select cosmological environment such as the value of the  Hubble parameter and the average redshift of the cluster, the mass or number density profiles, the velocity anisotropy profiles, as well as the optimization algorithm choices.
   * **H0 (real):** The value of the Hubble parameter evaluated at redshift z=0, measured in units of km/s/Mpc.
   * **z (real):** Average redshift of the cluster's center of mass.
   * **Omega_Lambda (real):** Value of the Omega_Lambda density parameter today.
   * **Omega_m (real):** Value of the Omega_m density parameter today.
   * **R_low (real):** Inner projected radius, defining the minimum projected radius from the cluster center (given in Mpc) at which galaxies should be considered in the MG-MAMPOSSt fit.
    
   * **R_up (real):** Outer projected radius, defining the maximum projected radius from the cluster center (given in Mpc) at which galaxies should be considered in the MG-MAMPOSSt fit.
   *  **kn (integer):** Number density profile model. It selects the model for the tracer's projected number density profile. In the current version three possible choices are allowed: projected NFW (kn=1), projected Hernquist (kn=2) and beta-model (kn=3).
    
   * **nb (real):** Exponent for beta-model. The beta-model choice requires an additional input exponent which should be negative. The input parameter is ignored otherwise.
    
   * **kmp (integer):** Mass profile/gravitational potential model. It selects the allowed parametrisations of the mass profile to be used in the Jeans' equation. For **kmp=1...6** and **kmp=10** GR is assumed, while if **kmp=7,8,9** the mass profile is the effective dynamical (Navarro-Frenk-White) mass obtained in linear Horndeski, Vainshtein screening and general chameleon screening respectively.  Details on each single model are given in the user manual.
    
   * **kani (integer):** Anisotropy profile model. It selects the  
    velocity anisotropy beta(r) in the Jeans' equation. The currently implemented profiles  are: constant anisotropy  **kani=0** , Mamon&Lokas **kani=1** , Osipkov-Merritt **kani=2**, simplified Wojtak **kani=3** , Tiret **kani=4** , modified Tiret **kani=5** .
    
    * **rcut (real):** Truncation radius. Value of the truncation radius needed in the pseudo-isothermal elliptical mass distribution (PIEMD), which correspond to kmp=3. It is ignored for other mass profile models.
    
    * **kbsp (integer): Fast mode**. If equal to 1, the likelihood is estimated by using a grid of values (default 60 points) in the phase space (R,v_z) of the galaxies and then bispline-interpolating over the data points. 
    
    * **kopt (integer): Optimization choice.** If required, an optimization algorithm is launched before the grid/MCMC likelihood computation to find the minimum of $-\ln \mathcal{L}$. Eventually, the resulting best fit parameters are used as guess values for the subsequent parameter space exploration. Currently, three choices are available:  BOBYQA **kopt=0**, NEWUOA **kopt=1** or POWELL **kopt=2** We point out that the POWELL algorithm does not work efficiently in MG mode (kmp=7,8,9) due to round-off errors in the integration routines. If **kopt=-1** the optimization is skipped. 
    
    * **kscr (integer): Screening mode** (available only for kmp=7). In linear Horndeski, one can choose to adopt the $f(R)$ sub-case, **kscr=0**, where A2 is fixed to 1/sqrt{6}. In this framework, there is the possibility to include a model-dependent screened $f(R)$ model with Hu&Sawicki functional form,implemented by assuming a simple analytical approximation. The transition between the screened and linear regime can be instantaneous (**kscr=1**), or smoothed with an additional parameter controlling the sharpness (**kscr=2**). For **kscr=-1**, the general linear Horndeski with two free parameters is selected.
     
* **Parameter limits (lines 29-40, reals):** limits in the parameter space exploration. It works only if the option **kpro=1** in the file Option.txt
   *  **r200_low:** Lower limit for r200 (in Mpc).
   *  **r200_up:** Upper limit for r200 (in Mpc).
   *  **rnu_low:** Lower limit for rnu (in Mpc).
   *  **rnu_up:** Upper limit for rnu (in Mpc).  
   *  **rs_low:** Lower limit for rs (in Mpc).
   *  **rs_up:** Upper limit for rs (in Mpc).  
   *  **beta_low:** Lower limit for for the anisotropy parameter beta. For kani=1 the unit is Mpc.
   *  **beta_up:** Upper limit for for the anisotropy parameter beta. For kani=1 the unit is Mpc.
   *  **A1_low:** Lower limit for the first modified gravity parameter A1. 
   *  **A1_up:** Upper limit for the first modified gravity parameter A1. 
   *  **A2_low:** Lower limit for the second modified gravity parameter A2. 
   *  **A2_up:** Upper limit for the second modified gravity parameter A2. 

### Working Opions

The file Option.txt contains various options and switches for the new features in MG-MAMPOSSt. These are mostly related to the numerical analysis and evaluation of the posterior likelihood. Notice that, the input parameters can be **binary integers** (with values 0 or 1), **integers*4** or **reals*8**. 
All the parameters must be given in a format **"label = "value"**. The "label"s are mandatory while the "value"s, if not given, are set by default. 
* **nmcmc (binary)**. Select between grid-search mode (= 0), and MCMC sampling (= 1). The default value is 1.
* **Nsample (integer)**. Number of points in the MCMC run. Default is 400000.
* **nlens (binary)**. Lensing information: If equal to 1, it adds to the kinematics likelihood, a probability distribution simulating additional information such as provided by a lensing mass profile reconstruction. For each set of values of the parameters, the joint (log) likelihood is then computed as

    **L(joint) =L(dyn)+L(lens),**

    where L(dyn) and L(lens) are the log-likelihoods of the MG-MAMPOSSt procedure and the simulated lensing distribution respectively. 
    For linear Horndeski and Chameleon screening, where photon propagation is not affected by the new degrees of freedom, the lensing likelihood has the form of a bivariate Gaussian distribution L(lens)=P(lens)(rs,r200), specified by central values, standard deviations and correlation (see below). In Vainsthein screening, where lensing mass profile is explicitly modified by the MG parameters, the likelihood is computed by simulating a full tangential shear profile, as explained in [Pizzuti et al., 2021](https://inspirehep.net/literature/1834246). The default value is 0.
    
 * **r200t (real)**. "True" value of the cluster's virial radius (in unit of Mpc) around which the lensing distribution is centered.
 * **rst (real)**.   "True" value of the cluster's scale radius (in unit of Mpc) around which the lensing distribution is centered.
 * **delta1 (real)**. For Vainshtein screening, it represents the intrinsic ellipticity of galaxies in the (weak) lensing simulations. For Chameleon screening, it is the relative uncertainty on the virial radius in the Gaussian distribution sigma(r_{200})/r_{200}.
 * **delta2 (real)**. For Vainshtein screening, it represents the large scale structure contribution to the errors in the (weak) lensing simulations. For Chameleon screening, it is the relative uncertainty on the scale radius in the Gaussian distribution \sigma(r_s)/r_s.
 * **delta3 (integer/real)**. For Vainshtein screening, it represents the number of galaxies per arcminute square in the (weak) lensing simulations. For Chameleon screening, it is the correlation in the Gaussian distribution.
 * **kpro (binary)**. If it is equal to 1, the parameter space exploration is made over a given interval (Xlow,Xup), where X indicates a generic free parameter. Default is 1. 
 * **Nclust (integer)**. MG-MAMPOSSt allows for efficient statistical forecast analyses of the constraints on the implemented MG models. In particular, it is possible to input Nclust realizations of phase spaces at the same time to compute directly the joint likelihood for a given set of parameters, obtained from the combination of the likelihood from each single data-set.  These data files should be located in /data folder and named as datphys_"i".dat, where "i" labels the number of the file in ascending order, starting from 2 (e.g. datphys_2.dat, datphys_3.dat, ...). The file format is the same as the main input datphys.dat. Default is 1. Note that in order to obtain meaningful results using this option, all the data files should be a realization of the same cluster (i.e. characterized by the same values of all parameters).
 * **nsingle (binary)**. If it is equal to 1, the Nclust-clusters likelihood is computed by simply multiplying by Nclust the likelihood from a single data-set. Useful to forecast a fast estimation of the limiting behaviour of the constraints when increasing the number of clusters. Default is 0.
   
All the required information about the structure of these files can be found in the [pdf user manual](https://arxiv.org/abs/2201.07194)


### Output
Output files can be found in the Output folder. Main is MaxLik_test.dat, 
which stores the tabulated $`-\ln\mathcal{L}`$, where $`mathcal{L}`$ is the 
likelihood/posterior, as a function of the model parameters.



Go to [Top](#top)

