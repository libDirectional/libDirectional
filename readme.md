libDirectional
==============

libDirectional is a library for directional statistics as well as recursive estimation on directional manifolds. The supported manifolds include 

  * unit circle
  * unit hypersphere
  * complex unit hypersphere  
  * torus
  * SE(2)

For each of these manifolds, the library contains certain probability distributions as well as recursive filtering algorithms.

Installation
------------

Requirements:

  * Matlab 2022a or later (earlier versions of libDirectional are compatible with older versions of Matlab)
  * a suitable compiler (Visual Studio 2015, Visual Studio 2017, MinGW64, gcc 4.7 or later, XCode)

To use libDirectional, add the entire lib-folder including subdirectories to Matlab's search path. Then, change to the lib-folder and run the `compileAll.m` script. This script should compile all mex-files used by libDirectional. 

If you experience any issues, run `mex -setup` and `mex -setup C++` to ensure that you have selected the correct compiler. In case you have trouble compiling with gcc, make sure that you are using a version that is officially supported by MATLAB.

The following toolboxes are recommended for libDirectional:

  * image_toolbox
  * optimization_toolbox
  * statistics_toolbox
  * symbolic_toolbox

The statistics toolbox and the optimization toolbox are fairly widely used, but the other toolboxes are only required for certain very specific features.

Example: Plotting Probability Density Functions
-----------------------------------------------

For example, we can generate a two-dimensional plot of the pdf of a wrapped normal distribution with parameters mu = 2 and sigma = 1.3 simply by typing the following two commands.

	>> wn = WNDistribution(2, 1.3);
	>> wn.plot2d();

We can then set the labels and axis using the following code: 

	>> setupAxisCircular('x');
	>> xlabel('x'); ylabel('f(x)');

Similarly, we can create plots of other distributions. A three-dimensional plot of the pdf of a von Mises distribution with parameters mu = 6 and kappa=0.5 could be generated using the following code.

	>> vm = VMDistribution(6, 0.5);
	>> vm.plot3d('color','red');
	>> hold on; vm.plotCircle('color','black'); hold off;
	>> xlabel('cos(x)'); ylabel('sin(x)'); zlabel('f(x)');

Example: Numerical and Analytical Calculation
---------------------------------------------

Let us again consider the wrapped normal distribution defined in the previous example. Suppose we want to calculate the first trigonometric moment, i.e., E(exp(ix)), of this distribution. For this purpose, we simply call the corresponding function:

	>> wn.trigonometricMoment(1)

This produces the output

	ans =
	  -0.1788 + 0.3906i

In the case of the wrapped normal distribution, `trigonometricMoment` is a function inside the class `WNDistribution` that implements an analytic calculation of the trigonometric moment. If no analytic solution was implemented, the function `trigonometricMoment` in the base class `AbstractCircularDistribution` would automatically fall back to an algorithm based on numerical integration. Even though an analytical solution is available for the wrapped normal distribution, we can still call the numerical algorithm as follows.

	>> wn.trigonometricMomentNumerical(1)

We obtain the result

	ans =
	  -0.1788 + 0.3906i

This can, for example, be used to compare the numerical and analytical results in order to ensure correctness of the analytical implementation. In this case, both results match up to the displayed number of digits, but in certain cases, analytical and numerical solutions may differ more significantly.

Example: Nonlinear Circular Filtering
-------------------------------------

Let us consider a system with circular state x_k in [0, 2pi) and dynamics

	x_{k+1} = a(x_k) + w_k 
	a_k(x_k) = x + 0.1 cos(x_k) mod 2 pi 

where w_k is distributed according to WN(x; 0, 0.4) If we assume that the current state is distributed according to WN(x; 2, 0.5), we can perform the prediction step with the WN-assumed filter using the following commands.

	>> filter = WNFilter();
	>> filter.setState(WNDistribution(2,0.5));
	>> a = @(x) mod(x + 0.1*cos(x),2*pi);
	>> filter.predictNonlinear(a, WNDistribution(0, 0.4));
	>> filter.getEstimate()

This produces the output

	ans = 
	  WNDistribution with properties:
	       mu: 1.9623
	    sigma: 0.6092

As you can see, the predicted density is returned as a wrapped normal distribution. Now we consider the measurement model
 
	z_k = h_k(x_k) + v_k
 
with

	h_k: [0,2 pi) -> R,  h_k(x) = sin(x)
 
where v_k is additive noise distributed according to N(x; 0, 0.7). As you can see, we have a circular state, but a real-valued measurement, in this case. However, a circular measurement (or a measurement on a completely different manifold) would be possible as well. If we obtain a measurement, say z = 0.3, we can perform the measurement update as follows.

	>> h = @(x) sin(x);
	>> measurementNoise = GaussianDistribution(0, 0.7);
	>> likelihood = LikelihoodFactory.additiveNoiseLikelihood(h, measurementNoise);
	>> filter.updateNonlinearProgressive(likelihood, 0.3)
	>> filter.getEstimate()

This produces the output

	ans = 
	  WNDistribution with properties:
	       mu: 2.0030
	    sigma: 0.6414

Once again, we obtain the result as a wrapped normal distribution.

Unit Tests
----------

The unit tests for libDirectional are located in the `tests` subfolder. You can automatically run all unit tests by executing `runLibDirectionalUnitTests`. Tests that are very computationally expensive are skipped by default. You can enable computationally expensive tests by running `runLibDirectionalUnitTests(true)`, which is more thorough but takes much more time.

Externals
---------

libDrectional relies on the following external libraries, which are also included in the externals folder.

  * [Eigen](http://eigen.tuxfamily.org/) (MPL 2)
  * [fmath](https://github.com/herumi/fmath) (BSD)
  * [mhg](http://www-math.mit.edu/~plamen/software/mhgref.html) (GLPv2 or later)
  * [Faddeva](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) (MIT)
  * [Nonlinear Estimation Toolbox](https://bitbucket.org/nonlinearestimation/toolbox) (GPLv3), only the necessary subset is included in externals folder
  * [TesselateS3](https://github.com/gerhardkurz/TessellateS3) (GPLv3)
  * Modified versions of [Slepian Alpha](http://csdms.colorado.edu/wiki/Model:SLEPIAN_Alpha) and [Slepian Bravo](http://csdms.colorado.edu/wiki/Model:SLEPIAN_Bravo) (GPLv2 or later)
  * [Wigner3j](https://de.mathworks.com/matlabcentral/fileexchange/20619-wigner3j-symbol) (BSD)
  * A modified version of the [Recursive Zonal Equal Area Sphere Partitioning Toolbox](http://eqsp.sourceforge.net/) (MIT)
  * [UKF-M](https://de.mathworks.com/matlabcentral/fileexchange/71994-arrow-3d) (BSD)
  * The farrow function of [Arrow 3d](https://de.mathworks.com/matlabcentral/fileexchange/71994-arrow-3d) (BSD)

Furthermore, we use the script `circVMcdf` by Shai Revzen (GPLv3). We also use some code from [libBingham](https://github.com/sebastianriedel/bingham) (BSD), but this library is not in the externals folder as only small parts are used.

License
-------

libDirectional is licensed under the GPLv3 license.

Citation
--------

If you use libDirectional in your research, please cite the library using as follows.

    @Article{libdirectional,
    author = {Gerhard Kurz and Igor Gilitschenski and Florian Pfaff and
              Lukas Drude and Uwe D. Hanebeck and Reinhold Haeb-Umbach
              and Roland Y. Siegwart},
    title = {Directional Statistics and Filtering Using {libDirectional}},
    year = {2019},
    journal = {Journal of Statistical Software},
    volume = {89},
    number = {4},
    pages = {1--31},
    doi = {10.18637/jss.v089.i04},
    }


Contact
-------

Lead author: Gerhard Kurz

Mail: kurz.gerhard (at) gmail.com

Web: [https://www.gerhardkurz.de](https://www.gerhardkurz.de)

Contributors:

  * Igor Gilitschenski
  * Florian Pfaff
  * Lukas Drude
