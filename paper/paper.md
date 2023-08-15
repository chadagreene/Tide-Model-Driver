---
title: 'Tide Model Driver for MATLAB'
tags:
  - MATLAB
  - tides
  - oceanography
  - TMD
authors:
  - name: Chad A. Greene
    orcid: 0000-0001-6710-6297 
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Svetlana Erofeeva
    orcid: 0000-0002-4489-7505 
    affiliation: 2
  - name: Laurie Padman
    orcid: 0000-0003-2010-642X
    affiliation: 3
  - name: Susan L. Howard
    orcid: 0000-0002-9183-0178 
    affiliation: 3
  - name: Tyler Sutterley
    orcid: 0000-0002-6964-1194  
    affiliation: 4
  - name: Gary Egbert
    affiliation: 2
affiliations:
 - name: Institute for Geophysics, John A. and Katherine G. Jackson School of Geosciences, University of Texas at Austin, Austin, TX, USA
   index: 1 
 - name: College of Earth, Ocean, and Atmospheric Sciences, Oregon State University, Corvallis, OR, USA
   index: 2
 - name: Earth and Space Research, Corvallis, OR, USA
   index: 3
 - name: University of Washington Applied Physics Laboratory Polar Science Center, Seattle, WA, 98122, USA
   index: 4
date: 4 August 2023
bibliography: paper.bib

---

# Summary

Astronomically-forced tides influence ocean surface height and currents on timescales of minutes to years. Tides contribute to ocean mixing [@munk:1998] and mean flows [@loder:1980], ice sheet and sea ice dynamics, and melting of ice shelves and marine terminating glaciers [@padman:2018]. Tidal signals must be accurately removed when calculating long-term trends in ocean surface height from tidally aliased satellite altimetry measurements [@smith:2000], and the tidal component of ocean circulation must be removed for analyses of ship-based current measurements [@carrillo:2005]. Several models have been made publicly available to predict tides on global [@stammer:2014] or regional [e.g., @padman:2004] scales. Each model contains only information about the complex coefficients of tidal constituents, and therefore requires software to extract and manipulate the data into a meaningful form.  

The Tide Model Driver for MATLAB version 3.0 (TMD3.0) allows users to access and calculate tidal predictions from model coefficients at arbitrary locations and times, and this version of the software represents decades of development by the tide community. The underlying equations for TMD were originally written in Fortran by Richard Ray, and were converted into MATLAB functions by Oregon State University and Earth and Space Research in 2005 [@padman:2022]. The MATLAB version of TMD developed a global user base well before the advent of GitHub or modern documentation standards, but no major updates have been implemented since its inception. The updated toolbox presented here has been restructured for computational efficiency and ease of use, but relies on the same mathematical equations that have been employed since its earliest implementation. The documentation has been greatly expanded to include clear descriptions of syntax along with many thoroughly explained, replicable examples that use real-world ocean data provided with the software. For TMD3.0, we also introduce a new, consolidated NetCDF tide model data format that is compact, user friendly, and can be used for any barotropic ocean or load tide model. TMD3.0 is issued as a free, open-source software package that is available to all. 

# Acknowledgements

A Fortran implementation of the Tide Model Driver software was originally written by Richard Ray, and we appreciate his feedback on the present software. 

# References
