Welcome to the README for SLEPIAN_Bravo, a set of routines for the
computation of spherical harmonics, Slepian functions, and transforms.
This package is hosted by CSDMS at
http://csdms.colorado.edu/wiki/Model:SLEPIAN_Bravo

SLEPIAN_Bravo mostly consists of routines that pertain to Simons and
Dahlen, 2006 paper in Geoph. J. Int. Please cite this reference and
perhaps the doi of the software itself if you find this package of
use:

Frederik J. Simons and F. A. Dahlen
Spherical Slepian functions and the polar gap in geodesy
Geoph. J. Int., 2006, 166 (3), 1039-1061, 
doi:10.1111/j.1365-246X.2006.03065.x


Where do I start?

We make two suggestions where users can begin to use this package.
First, this package contains a set of scripts which reproduce the
figures from Simons and Dahlen, 2006. The names of these scripts and
examples of the figures they produce can be found at
http://geoweb.princeton.edu/people/simons/software.html#GJI2006
Reproducing these figures yourself will help ensure you have set up
the code properly, and serve as examples to help users perform their
own analysis. Second, we suggest that users run the demos for the
functions XYZ2SLEP and PLM2XYZ. These demos illustrate various functions of the
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
