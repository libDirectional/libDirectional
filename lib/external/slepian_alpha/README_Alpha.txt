Welcome to the README for SLEPIAN_Alpha, a set of routines for the
computation of spherical harmonics, Slepian functions, and transforms.
This package is hosted by CSDMS at
http://csdms.colorado.edu/wiki/Model:SLEPIAN_Alpha

SLEPIAN_Alpha mostly consists of routines that pertain to Simons et
al., 2006 paper in SIAM Review. Please cite this reference and
perhaps the doi of the software itself if you find this package of
use:

Frederik J. Simons, F. A. Dahlen, and Mark A. Wieczorek
Spatiospectral Concentration on a Sphere
SIAM Rev., 48(3), 2006, 504-536.
http://dx.doi.org/10.1137/S0036144504445765


Where do I start?

We make two suggestions where users can begin to use this package.
First, this package contains a set of scripts which reproduce the
figures from Simons et al., 2006. The names of these scripts and
examples of the figures they produce can be found at
http://geoweb.princeton.edu/people/simons/software.html#SIAMReview
Reproducing these figures yourself will help ensure you have set up
the code properly, and serve as examples to help users perform their
own analysis. Second, we suggest that users run the demos for the
function LOCALIZATION. These demos illustrate various functions of the
code and precompute some useful data files. Several other functions
also have demos users may find useful.

Note: The demos above and the package in general expects a couple of
environmental variables to be set (e.g. your storage location:
$IFILES) and a directory tree for data storage to already be created
(i.e. subdirectories in $IFILES). The programs will initially display
errors instructing the user to set up these shell variables and
folders. Also, if you have an open pool of Matlab workers, a few
functions such as KERNELCP will take advantage of the parallel
computing resources and run much faster than otherwise.



Other important stuff

This software is distributed under the GNU Public License v2, which can be
found at http://geoweb.princeton.edu/people/simons/license.html  and also
copied below.

Copyright (C) 2014. Developer can be contacted by email. 

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version. 

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. 

You can receive a copy of the GNU General Public License by writing to the
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA. 
