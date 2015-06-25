classdef AbstractToroidalDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testAbstractToroidalDistribution(testCase)           
           mu = [1.3, 6.1]';
           C = [0.8 0.5; 0.5 0.9];
           twn = ToroidalWNDistribution(mu,C);
           
           %% test trigonometricMomentNumerical
           testCase.verifyEqual(twn.trigonometricMoment(0), twn.trigonometricMomentNumerical(0), 'RelTol', 1E-8);
           testCase.verifyEqual(twn.trigonometricMoment(1), twn.trigonometricMomentNumerical(1), 'RelTol', 1E-8);
           testCase.verifyEqual(twn.trigonometricMoment(2), twn.trigonometricMomentNumerical(2), 'RelTol', 1E-8);
           testCase.verifyEqual(twn.trigonometricMoment(3), twn.trigonometricMomentNumerical(3), 'RelTol', 1E-8);
           
           %% test circularMean
           testCase.verifyEqual(twn.circularMean, [twn.marginal(1).circularMean; twn.marginal(2).circularMean], 'RelTol', 1E-10);
           
           %% test correlation
           r1 = twn.circularCorrelationJupp();
           testCase.verifyGreaterThanOrEqual(r1, -1);
           testCase.verifyGreaterThanOrEqual(1, r1);
           
           r2 = twn.circularCorrelationJohnson;
           testCase.verifyGreaterThanOrEqual(r2, -1);
           testCase.verifyGreaterThanOrEqual(1, r2);
           
           testCase.verifyEqual(twn.circularCorrelationJammalamadaka(), twn.circularCorrelationJammalamadakaNumerical(),'RelTol', 1E-9);
           
           %% test mean4D
           me = twn.mean4D();
           mo = twn.trigonometricMoment(1);
           testCase.verifyEqual(me(1) + 1i*me(2),  mo(1));
           testCase.verifyEqual(me(3) + 1i*me(4),  mo(2));
           
           %% test covariance4DNumerical
           %compare to analytical solution
           CN = twn.covariance4DNumerical();
           testCase.verifyEqual(CN, CN');
           testCase.verifyGreaterThan(eig(CN), 0);
           testCase.verifyEqual(twn.covariance4D(), CN, 'RelTol', 1E-7);
           
           %compare to monte carlo solution
           rng default
           s = twn.sample(100000);
           S = [cos(s(1,:)); sin(s(1,:)); cos(s(2,:)); sin(s(2,:)); ];
           Cmontecarlo = cov(S');
           testCase.verifyEqual(Cmontecarlo, CN, 'RelTol', 1E-1); %quite inaccurate because it is a Monte Carlo solution
           
           %% test kld numerical
           testCase.verifyEqual(twn.kldNumerical(twn), 0);
           testCase.verifyGreaterThan(twn.kldNumerical(ToroidalWNDistribution(mu, eye(2,2))), 0);
           
           %% test integral numerical
           testCase.verifyEqual(twn.integral, 1, 'RelTol', 1E-9);
           testCase.verifyEqual(twn.integralNumerical, 1, 'RelTol', 1E-9);
           
           %% test sampling
           n = 10;
           s = twn.sampleMetropolisHastings(n);
           testCase.verifyEqual(size(s), [2, n]);
           testCase.verifyGreaterThanOrEqual(s, zeros(size(s)));
           testCase.verifyLessThan(s, 2*pi*ones(size(s)));
        end
    end
end