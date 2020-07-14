% This example evaluates stochastic sampling of circular densities.
%
% Three different sampling methods are compared, a native sample for the
% respective density, a generic Metropolis-Hastings sampler, and a sampler
% based on inverting the cumulative distribution function.

function stochasticSampling
    %dist = WNDistribution(3, 1.5);
    dist = VMDistribution(3, 1.5);
    %dist = WCDistribution(3, 1.5);
    
    nSamples = 1000;
    
    samples = dist.sample(nSamples);
    figure(1);
    compare(samples, dist, 'Native Sampling');

    figure(2);
    samples = dist.sampleMetropolisHastings(nSamples);
    compare(samples, dist, 'Metropolis-Hastings Sampling');

    figure(3);
    samples = dist.sampleCdf(nSamples);
    compare(samples, dist, 'Sampling Based on Inversion of CDF');
end

function compare(samples, distribution, titlestr)
    n = size(samples,2);
    steps = 50;
    s = zeros(1,steps);
    firstMomentError = zeros(1,steps);
    kuiper = zeros(1,steps);
    for i=1:steps
        s(i) = floor(i/steps*n);
        wd = WDDistribution(samples(:,1:s(i)));
        firstMomentError(i) = abs(wd.trigonometricMoment(1) - distribution.trigonometricMoment(1));
        kuiper(i) = wd.kuiperTest(distribution);
    end
    [hAx,~,~] = plotyy(s, firstMomentError, s, kuiper);
    xlabel('samples')
    ylabel(hAx(1),'first moment error') % left y-axis
    ylabel(hAx(2),'kuiper test') % right y-axis
    title(titlestr);

    wd = WDDistribution(samples);
    firstMomentErrorTotal = abs(wd.trigonometricMoment(1) - distribution.trigonometricMoment(1));
    kuiperTotal = wd.kuiperTest(distribution);
    
    fprintf('[%s] first moment error: %f\n', titlestr, firstMomentErrorTotal);
    fprintf('[%s] kuiper test: %f\n', titlestr, kuiperTotal);
    
    %llRatioWN = calculateLikelihoodRatio(samples, distribution, distribution.toWN())
    %llRatioVM = calculateLikelihoodRatio(samples, distribution, distribution.toVM())
    %llRatioWC = calculateLikelihoodRatio(samples, distribution, distribution.toWC())
    %llRatioUniform = calculateLikelihoodRatio(samples, distribution, CircularUniformDistribution)
end

function lr = calculateLikelihoodRatio(samples, distribution1, distribution2)
    lr = exp(sum(distribution1.logLikelihood(samples)) - sum(distribution2.logLikelihood(samples)));
end