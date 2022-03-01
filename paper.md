---
title: 'MG-MAMPOSSt, a code to test gravity at galaxy-cluster scales'
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
MG-MAMPOSSt is a license-free Fortran95 code that performs tests of General Relativity (GR) through the analyses of kinematical data of galaxy clusters based on the Jeans' equation. The code has been developed starting from the MAMPOSSt method of [@Mamon01], and extends it through new parametrisations of the gravitational potential for general families of gravity theories beyond GR aimed to explain the late-time accelerated expansion of the universe. By using input of projected positions and line-of-sight velocities of cluster's member galaxies, MG-MAMPOSSt reconstructs the cluster mass profile and the velocity anisotropy profile in modified gravity, jointly constraining the kinematics (mass and anisotropy profile) and modified gravity parameters. The code is further supplemented with a new capability to produce weak lensing forecasts for joint kinematic+lensing analyses, offering a valuable tool for studying the nature of gravity at cluster's scales.

# Statement of need

In the last two decades, the interest amongst the communities of cosmology and astrophysics in testing the nature of gravity and dark energy theories at large scales, aiming at explaining the origin of the late-time accelerated expansion of the universe ([@riess98],[@Perlmutter99]), has been growing up. Among the vast range of probes of modified gravity and dark energy models, galaxy clusters offer an interesting field of investigation at those scales where possible departures from General Relativity should become observable (e.g. [@Cataneo19review] and references therein). Cluster mass profiles ([@Wilcox:2015kna],[@Sakstein:2016ggl],[@Pizzuti:2016ouw]) and cluster abundance (e.g [@Lombriser:2010mp],[@Cataneo:2016iav]) have been extensively used in the literature to put stringent bounds on some popular classes of non-standard theories. In this context, we developed MG-MAMPOSSt, a FORTRAN95 code capable of performing tests of gravity models with the kinematics of member galaxies in clusters. The code is based upon the original MAMPOSSt method, developed by G. Mamon, A. Biviano and G. Boué ([@Mamon01], hereafter MAM13). A public version of MAMPOSSt by G. Mamon can be found at [https://gitlab.com/gmamon/MAMPOSSt](https://gitlab.com/gmamon/MAMPOSSt). Whereas the original MAMPOSSt code relies on the assumption of a standard Newtonian gravitational potential, MG-MAMPOSSt implements general and viable models of gravity beyond General Relativity (GR), with the aim of placing constraints on their theory space at galaxy-cluster scales, as well as of investigating the essential statistical degeneracy between model parameters. In addition, the code is capable of producing complementary weak-lensing forecasts towards joint kinematics+lensing analyses. 

MAMPOSSt (Modelling Anisotropy and Mass Profile of Spherical Observed Systems) determines mass profiles of galaxy clusters (or in general, spherical systems in dynamical equilibrium) by analysing the internal kinematics of the cluster members. MAMPOSSt has also been used for elliptical galaxies traced by globular clusters and  dwarf spheroidals traced by their stars [@Mamon+15], and has been recently extended into MAMPOSSt-PM to handle proper motions in star clusters ([@Mamon&Vitral22], see [@Vitral&Mamon21]).} Given an input of projected positions and line-of-sight (l.o.s) velocities of the member galaxies, and under the assumptions of spherical symmetry and dynamical relaxation, the code solves the Jeans equation to reconstruct the gravitational potential, the velocity anisotropy profile, and (optionally) the projected number density profile. 
MG-MAMPOSSt extends the method to gravity scenarios beyond GR, where the gravitational potential is explicitly modified by the presence of an additional scalar degree of freedom, resulting in effective mass profiles different from GR. The resulting effective mass profile can be confronted against real or synthetic data provided as input to the code. In addition, the code requires a parametric modelling for the profiles of velocity anisotropy, mass and number density, with several available choices which can be tuned as input. The code's main output is a tabulated likelihood as a function of all the free model parameters of gravity and other input physics.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
