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


# References
