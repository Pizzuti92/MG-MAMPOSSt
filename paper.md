---
title: 'MG-MAMPOSSt, a Fortran code to test gravity at galaxy-cluster scales'
tags:
  - Fortran
  - cosmology
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

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
\textsc{MG-MAMPOSSt} is a license-free \textsc{Fortran95} code that performs tests of General Relativity (GR) through the analyses of kinematical data of galaxy clusters based on the Jeans' equation. The code has been developed starting from the \textsc{MAMPOSSt} method of [@Mamon01], and extends it through new parametrisations of the gravitational potential for general families of gravity theories beyond GR aimed to explain the late-time accelerated expansion of the universe. By using input of projected positions and line-of-sight velocities of cluster's member galaxies, \textsc{MG-MAMPOSSt} reconstructs the cluster mass profile and the velocity anisotropy profile in modified gravity, jointly constraining the kinematics (mass and anisotropy profile) and modified gravity parameters. The code is further supplemented with a new capability to produce weak lensing forecasts for joint kinematic+lensing analyses, offering a valuable tool for studying the nature of gravity at cluster's scales.

# Statement of need

In the last two decades, the interest amongst the communities of cosmology and astrophysics in testing the nature of gravity and dark energy theories at large scales, aiming at explaining the origin of the late-time accelerated expansion of the universe ([@riess98],[@Perlmutter99]), has been growing up. Among the vast range of probes of modified gravity and dark energy models, galaxy clusters offer an interesting field of investigation at those scales where possible departures from General Relativity should become observable (e.g. [@Cataneo19review] and references therein). Cluster mass profiles ([@Wilcox:2015kna],[@Sakstein:2016ggl],[@Pizzuti:2016ouw]) and cluster abundance (e.g [@Lombriser:2010mp],[@Cataneo:2016iav]) have been extensively used in the literature to put stringent bounds on some popular classes of non-standard theories. In this context, we developed \textsc{MG-MAMPOSSt}, a FORTRAN95 code capable of performing tests of gravity models with the kinematics of member galaxies in clusters. The code is based upon the original \textsc{MAMPOSSt} method, developed by G. Mamon, A. Biviano and G. Boué ([@Mamon01], hereafter MAM13). A public version of \textsc{MAMPOSSt} by G. Mamon can be found at [https://gitlab.com/gmamon/MAMPOSSt](https://gitlab.com/gmamon/MAMPOSSt). Whereas the original code relies on the assumption of a standard Newtonian gravitational potential, \textsc{MG-MAMPOSSt} implements general and viable models of gravity beyond General Relativity (GR), with the aim of placing constraints on their theory space at galaxy-cluster scales, as well as of investigating the essential statistical degeneracy between model parameters. In addition, the code is capable of producing complementary weak-lensing forecasts towards joint kinematics+lensing analyses. \textsc{MG-MAMPOSSt} has been first presented and applied to real galaxy cluster data in  [@Pizzuti:2017diz], then extended and upgraded by [@Pizzuti2021].

\textsc{MAMPOSSt} (Modelling Anisotropy and Mass Profile of Spherical Observed Systems) determines mass profiles of galaxy clusters (or in general, spherical systems in dynamical equilibrium) by analysing the internal kinematics of the cluster members. Given an input of projected positions and line-of-sight (l.o.s) velocities of the member galaxies, and under the assumptions of spherical symmetry and dynamical relaxation, the code solves the Jeans equation to reconstruct the gravitational potential, the velocity anisotropy profile - which measures the difference among the velocity dispersion components in each direction - and (optionally) the projected number density profile. 
\textsc{MG-MAMPOSSt} extends the method to gravity scenarios beyond GR, where the gravitational potential is explicitly modified by the presence of an additional scalar degree of freedom, resulting in effective mass profiles different from GR. The effective mass profile is confronted against real or synthetic data provided as input to the code to infer the value of the free parameters describing the gravity model. 
In addition to the modified mass profile, the code requires a parametric modelling for the profiles of velocity anisotropy and number density of the member galaxies, with several available choices which can be tuned as input. The code's main output is a tabulated likelihood/posterior as a function of all the free model parameters of gravity and other input physics. The code takes $\sim 1$ second to find the best fit values of the parameters for the selected models, and few hours to perform a complete run of $\sim 10^5$ sampling of the posterior through a  simple (although efficent) Monte Carlo-Markov Chain (MCMC) exploration of the parameter space.

# Mathematics

\textsc{MG-MAMPOSSt} operates relying on the input of the projected phase space $(R,v_\text{z})$ of the cluster member galaxies. Here, $R$ is the projected distance from the cluster center at which a galaxy is seen by the observer, and  $v_\text{z}$ the velocity measured along the l.o.s. in the rest frame of the cluster. As mentioned before, the current version of the code assumes parametric expressions for the various kinematical quantities. Moreover, the 3-dimensional velocity distribution is taken to be a Gaussian. The latter assumption has been well-tested through cosmological simulations, as explained in the original MAMPOSSt paper (MAM13). 

The output likelihood is computed by comparing data of galaxies in projected phase space $(R,v_\text{z})$ to the theoretical radial velocity dispersion of the cluster member galaxies, obtained for a given set of models and parameters as a solution of the spherical Jeans' equation (see e.g. [@MamLok05]),
\begin{equation}
\label{eq:sigmajeans}
\sigma^2_r(r)=\frac{1}{\nu(r)}\int_r^{\infty}{\exp\left[2\int_r^s{\frac{\beta(t)}{t}\text{d}t}\right]\nu(s)\frac{\text{d}\Phi}{\text{d}s}\text{d}s},
\end{equation}
projected in the phase space. The above equation captures the main necessary input required for the \textsc{MG-MAMPOSSt} method: the gradient of the gravitational potential, which in turn contains information about the model of gravity, as well as the velocity anisotropy profile ($\beta(r)$) and the projected number density of tracers ($\nu(r)$).

The current version of the code can handle up to a six-dimensional parameter space, with two parameters defining the mass profile, one parameter for the velocity anisotropy profile, one for the number density profile, and finally, two parameters related to the modified gravity framework. Each parameter can be either treated as free in the fitting procedure, or it can be assigned pre-defined values. 

In the original code of MAM13 there are several possible choices for the modelling of the dark matter mass profile in the $\Lambda$CDM scenario. \textsc{MG-MAMPOSSt} adds new parametrizations to handle popular modifications of gravity. At the moment, all the implemented non-standard profiles rely on the Navarro-Frenk-White (NFW, [@navarro97]) mass density profile to model the matter density distribution, which has been shown to provide a good description for simulated and observed galaxy clusters, both in GR and in modified gravity (see e.g. [@Umetsu20], [@Peirani17],[@Wilcox:2016guw]). Nevertheless, other mass models are going to be included in the code with upcoming versions. Since the mass profile, or equivalently, the gravitational potentials, enters only in the expression of the radial velocity dispersion, \autoref{eq:sigmajeans}, the implementation of new models can be performed directly by the user, modifying the subroutines where the above equation is involved in the main source code. In particular, the functions \texttt{sr2int(alr)} and \texttt{fa(tlog)} are the only parts of the \textsc{MG-MAMPOSSt} code where parametrizations for the mass profile appear. 

The code is equipped with the two most popular and observationally viable  classes of dark energy models beyond GR based on a single, extra scalar field. These correspond to the so--called chameleon models and Beyond Horndeski/DHOST models. Common ground between the two families of models is the presence of the extra dynamical scalar degree of freedom ($\phi$) which introduces a new gravitational force. However, the structure of the fifth force  in each of them exhibits different characteristics; in particular, the screening mechanism, introduced to recover standard gravity at small scales and high density regions, acts in a different way producing very peculiar imprint on the gravitational potential. Both families are characterized by two free parameters, determining the action of the fifth force, which can be constrained with \textsc{MG-MAMPOSSt}. For a detailed exposition of these models and the associated equations we refer to our main paper [@Pizzuti2021], as well as to the original papers where the models were  first introduced [@Kobayashi:2014ida],[@Crisostomi:2017lbg],[@Dima:2017pwp].

As shown in [@Pizzuti2021], internal kinematics alone is generally not enough to provide stringent bounds on the modified gravity parameters, due to the strong degeneracy between model parameters. For this reason, \textsc{MG-MAMPOSSt} gives the possibility to include a simulated lensing information to your kinematics analysis in modified gravity, a feature which is particularly useful for forecasting the constraining power of the method. 


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

LP is partially supported by a 2019 "Research and Education" grant from Fondazione CRT. The OAVdA is managed by the Fondazione Cle\'ment Fillietroz-ONLUS, which is supported by the Regional Government of the Aosta Valley, the Town Municipality of Nus and the "Unite\' des Communes valdotaines Mont-E\'milius.
I. D. Saltas is supported by the Grant Agency of the Czech Republic (GAČR), under the grant number 21-16583M. The authors further acknowledge all the developers of the free FORTRAN routines used in \textsc{MG-MAMPOSSt}. Credits are given in the header of each routine.

# References
