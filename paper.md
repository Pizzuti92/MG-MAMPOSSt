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






# References
