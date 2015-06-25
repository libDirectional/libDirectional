% Example for comparing different parameter estimation methods for WN
% distributions.

% Initial WN distribution
wn = WNDistribution(4,0.4);

% Obtain stochastic samples
s = wn.sample(100);

% Estimate parameters using four different methods
wnMleJensen = WNDistribution.mleJensen(s)

wnMleNumerical = WNDistribution.mleNumerical(s)

wd = WDDistribution(s);
wnMoment = wd.toWN()

wnUnwrappingEM = wd.toWNunwrappingEM()

% Compare logLikelihood of the samples
wnMleJensen.logLikelihood(s)
wnMleNumerical.logLikelihood(s)
wnMoment.logLikelihood(s)
wnUnwrappingEM.logLikelihood(s)
