TessellateS3
============

TessellateS3 is an efficient implementation of the tetrahedral/octahedral tessellation algorithm for the four-dimensional unit sphere as a Mex File for MATLAB. It is based on the algorithm by Schaefer et al. and the C implementation by Jared Glover contained in [libBingham](https://github.com/SebastianRiedel/bingham). All parts of libBingham that are not required for the S3 tessellation have been removed.

Requirements
------------

  * a reasonably recent version of MATLAB
  * a C++ Compiler (tested with Visual Studio 2015)

Features
--------

  * compute a tessellation of S3, the four-dimensional unit sphere
  * easy to use interface for MATLAB
  * variable number of subdivision levels

Example usage
-------------
First, compile the mex file by calling

	>> compileTessellate

Then, you can run the tessellation algorithm using

	>> x = tessellate_S3(n)

Here, n corresponds to the number of points desired. The algorithm will return the next larger possible tessellation. Possible tessellations have 16, 128, 1024, 8192, etc. points. You can try the included unit test using

	>> run(TessellateS3Test)

References
----------
	Schaefer, S. and Hakenberg, J. and Warren, J.,
	Smooth Subdivision of Tetrahedral Meshes
	Proceedings of the 2004 Eurographics/ACM SIGGRAPH Symposium on Geometry Processing
	2004

	Jared Glover and Leslie Pack Kaelbling
	Tracking 3-D Rotations with the Quaternion Bingham Filter
	MIT-CSAIL-TR-2013-005
	2013

	Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
	Discretization of SO(3) Using Recursive Tesseract Subdivision
	Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017), Daegu, Korea, November 2017.

License
-------

TessellateS3 is licensed under the GPLv3 license.

Contact
-------

Author: Gerhard Kurz

Mail: gerhard.kurz (at) kit (dot) edu

Web: [http://isas.uka.de/User:Kurz](http://isas.uka.de/User:Kurz)

