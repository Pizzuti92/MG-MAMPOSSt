---
title: 'MG-MAMPOSSt, a Fortran code to test gravity at galaxy-cluster scales'
tags:
  - Fortran
  - astronomy
  - cosmology
  -  gravitation
  -  mass-modeling 
  -  orbital shapes
  -  kinematics
  - galaxy clusters
  - modified gravity
authors:
  - name: Lorenzo Pizzuti^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0001-5654-7580
    affiliation: "1" 
  - name: Ippocratis D. Saltas 
    affiliation: 2
  - name: Andrea Biviano
    affiliation: "3,4"
  - name: Gary Mamon
    affiliation: 5
  - name: Luca Amendola
    affiliation: 6
affiliations:
 - name: Osservatorio Astronomico della Regione Autonoma Valle d'Aosta,  Loc. Lignan 39, I-11020, Nus, Italy
   index: 1
 - name: CEICO, Institute of Physics of the Czech Academy of Sciences, Na Slovance 2, 182 21 Praha 8, Czechia
   index: 2
 - name: INAF, Osservatorio Astronomico di Trieste, via Tiepolo 11, 34143 Trieste, Italy
   index: 3
 - name:  IFPU, Institute for Fundamental Physics of the Universe, via Beirut 2, 34014 Trieste, Italy
   index: 4
 - name:  Institut d'Astrophysique de Paris (UMR 7095, CNRS and Sorbonne Université), 98 bis Bd Arago, F-75014 Paris, France
   index: 5
 - name:  Institute of Theoretical Physics, Philosophenweg 16, Heidelberg University, 69120, Heidelberg, Germany
   index: 6
date: 05 February 2022
bibliography: paper.bib

output:
  pdf_document:
    citation_package: natbib
  bookdown::pdf_book:
    citation_package: biblatex
biblio-style: "mnras"
link-citations: true

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
\textsc{MG-MAMPOSSt} is a license-free \textsc{Fortran90} code to perform tests of General Relativity (GR) through the analyses of kinematical data of galaxy clusters. The code solves the Jeans equation, relying on the \textsc{MAMPOSSt} method of [@Mamon01]. It extends the latter method through new parametrisations of the gravitational potential for general families of gravity theories beyond GR aimed to explain the late-time accelerated expansion of the universe (@riess98, @Perlmutter99).Through appropriate input of projected positions and line-of-sight velocities of cluster's member galaxies, \textsc{MG-MAMPOSSt} reconstructs the cluster mass profile and the velocity anisotropy profile in modified gravity, jointly constraining the kinematics (mass and anisotropy profile) and modified gravity parameters. The code is further supplemented with a new capability to produce weak lensing forecasts for joint kinematic+lensing analyses, offering a valuable tool for studying the nature of gravity at cluster's scales.

# Statement of need

In the last two decades, the interest amongst the communities of cosmology and astrophysics in testing dark energy theories at large scales has been grown up significantly. Among the vast range of probes of modified gravity and dark energy models, galaxy clusters offer a powerful laboratory at scales where possible departures from General Relativity should become observable (e.g. @Cataneo19review and references therein). Cluster mass profiles (@Wilcox15, @Sakstein:2016ggl, @Pizzuti:2016ouw) and the cluster abundance (e.g @Lombriser:2010mp, @Cataneo:2016iav) have been extensively used in the literature to place stringent bounds on some popular classes of extensions of GR. In this context, we developed \textsc{MG-MAMPOSSt}, a \textsc{Fortran90} code capable of performing tests of gravity models with the kinematics of member galaxies in clusters. The code is based upon the original \textsc{MAMPOSSt} method, developed by G. Mamon, A. Biviano and G. Boué (@Mamon01, hereafter MAM13). A public version of \textsc{MAMPOSSt} by G. Mamon can be found at [https://gitlab.com/gmamon/MAMPOSSt](https://gitlab.com/gmamon/MAMPOSSt). Whereas the original code relies on the assumption of a standard Newtonian gravitational potential, \textsc{MG-MAMPOSSt} implements general and viable models of gravity beyond GR, with the aim of placing constraints on their theory space and investigating the essential statistical correlations between model parameters. In addition, the code is capable of producing complementary weak-lensing forecasts for joint kinematics+lensing analyses. \textsc{MG-MAMPOSSt} has been first presented and applied to real galaxy cluster data in @Pizzuti:2017diz, then extended and upgraded by @Pizzuti2021.

\textsc{MAMPOSSt} (Modelling Anisotropy and Mass Profile of Spherical Observed Systems) determines mass profiles of galaxy clusters (or in general, spherical systems in dynamical equilibrium) by analysing the internal kinematics of the cluster members. Given an input of projected positions and line-of-sight (l.o.s) velocities of the member galaxies, and under the assumptions of spherical symmetry and dynamical relaxation, the code solves the Jeans equation to reconstruct the gravitational potential, the velocity anisotropy profile - which measures the difference among the velocity dispersion components in each direction - and (optionally) the projected number density profile. 
\textsc{MG-MAMPOSSt} extends the method to the case where the gravitational potential is explicitly modified by the presence of an additional scalar degree of freedom. The non-standard gravitational potential is confronted against real or synthetic data provided as input to the code to infer the value of the free parameters describing the gravity model. 
In addition to the modified mass profile, the code requires a parametric modelling for the profiles of velocity anisotropy and number density of the member galaxies, with several available choices which can be tuned as input. The code's main output is a tabulated likelihood/posterior as a function of all the free model parameters of gravity and other input physics. The code takes $\sim 1$ second to find the best fit values of the parameters for the selected models, and few hours to perform a complete run of $\sim 10^5$ sampling of the posterior through a  simple (although efficent) Monte Carlo-Markov Chain (MCMC) exploration of the parameter space.

# Mathematics

\textsc{MG-MAMPOSSt} operates relying on the input of the projected phase space $(R,v_\text{z})$ of the cluster member galaxies. Here, $R$ is the projected distance from the cluster center at which a galaxy is seen by the observer, and  $v_\text{z}$ the velocity measured along the l.o.s. in the rest frame of the cluster. The output likelihood is computed by comparing data of galaxies in projected phase space $(R,v_\text{z})$ to the theoretical radial velocity dispersion of the cluster member galaxies, obtained for a given set of models and parameters as the solution of the spherical Jeans' equation (see e.g. @MamLok05),
\begin{equation}
\label{eq:sigmajeans}
\sigma^2_r(r)=\frac{1}{\nu(r)}\int_r^{\infty}{\exp\left[2\int_r^s{\frac{\beta(t)}{t}\text{d}t}\right]\nu(s)\frac{\text{d}\Phi}{\text{d}s}\text{d}s},
\end{equation}
projected in the phase space. The above equation captures the main necessary input required for the \textsc{MG-MAMPOSSt} method: the gradient of the gravitational potential, which in turn contains information about the model of gravity, as well as the velocity anisotropy profile ($\beta(r)$) and the projected number density of tracers ($\nu(r)$).

The current version of the code can handle up to a six-dimensional parameter space, with two parameters defining the mass profile, one parameter for the velocity anisotropy profile, one for the number density profile, and finally, two parameters related to the modified gravity framework. Each parameter can be either treated as free in the fitting procedure, or it can be assigned pre-defined values. 

In the original code of MAM13 there are several possible choices for the modelling of the dark matter mass profile in the $\Lambda$CDM scenario. \textsc{MG-MAMPOSSt} adds new parametrizations to handle popular modifications of gravity. At the moment, all the implemented non-standard profiles rely on the Navarro-Frenk-White (NFW hereafter, @navarro97) mass density profile to model the matter density distribution, which has been shown to provide a good description for simulated and observed galaxy clusters, both in GR and in modified gravity (see e.g. @Umetsu20, @Peirani17, @Wilcox:2016guw). Nevertheless, other mass models are going to be included in the code with upcoming versions. Since the mass profile, or equivalently, the gravitational potential, enters only in the expression of the radial velocity dispersion, \autoref{eq:sigmajeans}, the implementation of new models can be performed directly by the user, modifying the subroutines in the auxiliary file \texttt{MAM.f} where the above equation is involved in a straightforward way. In particular, the functions called \texttt{sr2int(alr)} and \texttt{fa(tlog)} are the only parts of the \textsc{MG-MAMPOSSt} code where parametrizations for the mass profile appear. 

The code is equipped with the two most popular and observationally viable  classes of dark energy models beyond GR based on a single, extra scalar field. These correspond to the so--called chameleon models and Beyond Horndeski/DHOST models. Common ground between the two families of models is the presence of the extra dynamical scalar degree of freedom ($\phi$) which introduces a new gravitational force modifying the (gradient of the) gravitational potential as (e.g. @Saltas:2019ius):
\begin{equation}
\frac{\text{d}\Phi}{\text{d}r}=\frac{G}{r^2}\left[M(r,\bold{\theta}_{DM})+f(r,\bold{\theta}_{DM},\bold{\theta}_\text{MG})\right] .
\end{equation}
In the above equation, $M(r,\bold{\theta}_{DM})$ is the total mass profile at radius $r$, as a function of the parameter vector $\bold{\theta}_{DM}$, while $f(r,\bold{\theta}_{DM},\bold{\theta}_\text{MG})$ is the contribution of the fifth force, which depends on the parametrisation of the mass density and on the parameters defining the modified gravity models $\bold{\theta}_\text{MG}$.
Both families are characterized by two free parameters determining the action of the fifth force, which can be constrained with \textsc{MG-MAMPOSSt}. For a detailed exposition of these models and the associated equations we refer to our main paper [@Pizzuti2021], as well as to the original papers where the models were  first introduced (@Kobayashi:2014ida, @Crisostomi:2017lbg, @Dima:2017pwp).

As shown in [@Pizzuti2021], internal kinematics alone is generally not enough to provide stringent bounds on the modified gravity parameters, due to the strong degeneracy between model parameters. For this reason, \textsc{MG-MAMPOSSt} gives the possibility to include a simulated lensing information to your kinematics analysis in modified gravity, a feature which is particularly useful for forecasting the constraining power of the method in view of upcoming imaging and spectroscopic surveys such as Euclid or LSST. 

# Functionality and Design

A complete run of the \textsc{MG-MAMPOSSt} code is based upon several correlated files which store the input/output information. 
In particular:

-  \texttt{gomamposstopt}_\texttt{x.inp} contains the names and locations of the input data file and of the input parameter file, as well as the names and locations of the output files. Each of them can be customised by the user

-  \texttt{data/datphys.dat} is the input data file, structured as a table where the number of rows coincides with the number of data points. The first column is the projected radius in units of $\text{kpc}$, the second and thirds columns represent the l.o.s. velocities and the associated errors in units of km/s.

-  \texttt{input\textunderscore pars/pars\textunderscore all\textunderscore N\textunderscore O\textunderscore spec\textunderscore DS} is the input parameters file, where one can select the number of free parameters and their guess values, the models of the various kinematic components (gravitational potential, number density profile and velocity anisotropy profile) and other relevant physical quantities for the \textsc{MG-MAMPOSSt} analysis.

-  \texttt{Options.txt} contains additional options required by \textsc{MG-MAMPOSSt}, e.g. how to explore the parameter space (fixed grid of values or MCMC) and the details of the lensing simulation in modified gravity.

- \texttt{Output/MaxLik.dat} is the main output. It is organized as a table where each row indicates the values of the parameters for a given point in the six-dimensional parameter space and the corresponding value of the logarithm of the Likelihood/Posterior.

The \textsc{MG-MAMPOSSt} run further produces additional output files, stored in the \texttt{Output} folder and, optionally, a plot of the marginalized posteriors for the free parameters when the MCMC exploration mode is selected. The plots generation requires the Python [getdist package](https://github.com/cmbant/getdist) [@Lewis:2019xzd]. In \autoref{fig:example} an example of a typical output plot is shown for a MCMC sampling in the case of DHOST model of gravity with five free parameters. The run has been performed by using the sample data-set shipped together with the code; the test of execution described in the documentation should produce in the same Figure. For a complete description of the code basic usage and functionalities, see the code's [documentation](https://github.com/Pizzuti92/MG-MAMPOSSt/blob/main/Documentation.pdf) or refer to [@Pizzuti22man]. 

![Example of the marginalized distribution for the free parameters of the DHOST modified gravity model implemented in \textsc{MG-MAMPOSSt}, obtained by using the test data-set provided in the code repository. The MCMC sampling has been performed over $10^5$ points in the parameter space. Dark blue and light blue areas corresponds to $1\sigma$ and $2 \sigma$ confidence region, respectively. $r_s$ and $r_{200}$ are the mass profile parameters, $\mathcal{A}_\infty$ is the velocity anisotropy profile parameter and $Y_1$, $Y_2$ are the parameters defining the modified gravity model analysed  .\label{fig:example}](https://github.com/Pizzuti92/MG-MAMPOSSt/blob/main/test/plot_example_test1.png)


# Acknowledgements

LP is partially supported by a 2019 "Research and Education" grant from Fondazione CRT. The OAVdA is managed by the Fondazione Cle\'ment Fillietroz-ONLUS, which is supported by the Regional Government of the Aosta Valley, the Town Municipality of Nus and the "Unite\' des Communes valdotaines Mont-E\'milius.
I. D. Saltas is supported by the Grant Agency of the Czech Republic (GAČR), under the grant number 21-16583M. The authors further acknowledge all the developers of the free FORTRAN routines used in \textsc{MG-MAMPOSSt}. Credits are given in the header of each routine.

# References
