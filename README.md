# Robust current regulator design for LCL-VSI based on D-decomposition and Parameter Space approach

## Overview 
This repository contains the companion MATLAB code for the conference paper:

I. Tyuryukanov, M. Popov, "D-Decomposition Based Robust Discrete-Time Current Regulator for Grid-Connected VSI" in Proc. 29th IEEE International Symposium on Industrial Electronics (ISIE 2020),
2020, pp. 1–8.

The code allows to generate robust low-order current regulators for LCL filter interfaced voltage source inverters (VSI).
The controller synthesis is based on D-Decomposition and parameter space approach. 
The repository contains multiple files that implement the two abovementioned design techniques, as well as some extra auxiliary MATLAB files. 

The file **start.m** is the program entry calling all other necessary function. 
The file **start.m** can be modified by the user familiar with MATLAB to achieve different results.
The file **test.m** contains some generic tests that are non-essential.

Although the code in the repository was designed with due care to avoid errors, unforeseen bugs and errors can't be excluded.

The implemented robust controller design approach relies on designing a controller satisfying multiple representatives of the uncertain plant.
Therefore, the resulting controller has to be tested against a sufficiently large set of uncertain plant representatives to verify its robustness properties. 
A simple (but not unique) way to accomplish this is to plot the closed-loop eigenvalues for a large number of uncertain plant representatives (see **evals_plot.m**).

- - - -

## Basic Requirements
For the reasons of numerical accuracy and robustness, the operations with polynomials in D-decomposition rely on multiprecision computations.
This implementation relies on the Advanpix Multiprecion Computing Toolbox (https://www.advanpix.com/) for computing with arbitrary precision. 

In addition, the controller synthesis program relies on the following third-party scripts (included into the repository):

1. X. Beudaert. Range intersection. 2012. url: http://www.mathworks.com/matlabcentral/fileexchange/31753 (visited on 05/14/2014).
2. D. Engwirda and S. Paris. Fast points-in-polygon test. 2011. url: http://www.mathworks.com/matlabcentral/fileexchange/15410-inpoly-mex-file (visited on 05/14/2014).
3. A. Murta. GPC: General Polygon Clipper library. n.d. url: http://www.cs.man.ac.uk/~toby/gpc/ (visited on 05/14/2014).
4. S. Hölz. Polygon Clipper. 2006. url: https://www.mathworks.com/matlabcentral/fileexchange/8818-polygon-clipper (visited on 17/05/2020).
5. D. Schwarz. Fast and Robust Curve Intersections. 2010. url: http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections (visited on 05/14/2014).
6. M. Yoshpe and A. Weinstein. Distance from a point to polygon. 2008. url: http://www.mathworks.com/matlabcentral/fileexchange/19398  (visited on 05/14/2014)
7. P. J. Acklam. Chebpoly & Chebpoly2. 2004. 
8. K. J�nasson. RGB. 2009

All credits for these scripts belong to their respective creators. 
- - - -

## Usage
Run **start.m** from MATLAB.
To understand the internal workings of **start.m**, step through it with the interactive MATLAB debugger.
To change the case study, change the relevant lines of **start.m**.

- - - -

## License
BSD 3-clause license.

- - - -

## Authors
Ilya Tyuryukanov
