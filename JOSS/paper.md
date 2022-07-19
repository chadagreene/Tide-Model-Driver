---
title: 'Tide Model Driver for MATLAB'
tags:
  - MATLAB
  - tides
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
  - name: Susan Howard
    orcid: 0000-0002-9183-0178 
    affiliation: 3
  - name: Tyler Sutterley
    orcid: 0000-0002-6964-1194  
    affiliation: 4
  - name: Gary Egbert
    affiliation: 2
affiliations:
 - name: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA
   index: 1
 - name: College of Earth, Ocean, and Atmospheric Sciences, Oregon State University, Corvallis, OR, USA
   index: 2
 - name: Earth and Space Research, Corvallis, OR, USA
   index: 3
 - name: University of Washington Applied Physics Laboratory Polar Science Center, Seattle, WA, 98122, USA
   index: 4
date: 19 July 2022
bibliography: paper.bib

---

# Summary

The Tide Model Driver for MATLAB version 3.0 (TMD3.0) predicts tidal heights or tidal current transports by summing the coefficients of tidal constituents at specified times and locations. The underlying equations for the Tide Model Driver were originally written in Fortran, and were translated into MATLAB nearly two decades ago. This update uses the same underlying equations, but the package has been restructured for computational efficiency and ease of use. For TMD3.0, we introduce a new, consolidated NetCDF tide model data format that is compact, user friendly, and can be used for any barotropic ocean or load tide model. 


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

A Fortran implementation of the Tide Model Driver software was originally written by Richard Ray, and we appreciate his feedback on the present software. 

# References