classdef WNDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWNDistribution(testCase)
            mu = 3;
            sigma = 1.5;
            wn = WNDistribution(mu,sigma);
            
            %% test pdf
            function p = truePdf(x)
                sum = 0;
                for k = -20:20
                    sum = sum + exp(-(x-mu + 2*pi*k).^2/2/sigma^2);
                end
                p = 1/sqrt(2*pi)/sigma * sum;
            end
            testCase.verifyEqual(wn.pdf(mu), truePdf(mu), 'RelTol', 1E-10);
            testCase.verifyEqual(wn.pdf(mu-1), truePdf(mu-1), 'RelTol', 1E-10);
            testCase.verifyEqual(wn.pdf(mu+2), truePdf(mu+2), 'RelTol', 1E-10);
            x = 0:6;
            testCase.verifyEqual(wn.pdf(x), arrayfun(@truePdf,x), 'RelTol', 1E-10);
            
            %% test pdf with large sigma
            wnLargeSigma = WNDistribution(0,100);
            fx = ones(size(x))/2/pi;
            testCase.verifyEqual(wnLargeSigma.pdf(x), fx, 'RelTol', 1E-10);
            
            %% test cdf
            fx = arrayfun(@(y) wn.integral(0, y), x);
            testCase.verifyEqual(wn.cdf(x), fx, 'RelTol', 1E-10);
            startingPoint = 3;
            fx = arrayfun(@(y) wn.integral(startingPoint, startingPoint+mod(y-startingPoint, 2*pi)), x);
            testCase.verifyEqual(wn.cdf(x,startingPoint), fx, 'RelTol', 1E-10);

            %% test integral
            testCase.verifyEqual(wn.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wn.integral(0,pi)+wn.integral(pi,2*pi), 1, 'RelTol', 1E-10);
            
            %% test angular moments
            testCase.verifyEqual(wn.trigonometricMoment(1), wn.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(wn.trigonometricMoment(2), wn.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(wn.trigonometricMoment(3), wn.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            testCase.verifyEqual(wn.circularMean, mu);
            
            %% test fromMoment
            wnFromMoment = WNDistribution.fromMoment(wn.trigonometricMoment(1));
            testCase.verifyClass(wnFromMoment, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnFromMoment.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wn.sigma, wnFromMoment.sigma, 'RelTol', 1E-10);
            
            %% test conversions based on angular moment matching
            wd2 = wn.toDirac2();
            testCase.verifyClass(wd2, 'WDDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), wd2.trigonometricMoment(1), 'RelTol', 1E-10);
            wd3 = wn.toDirac3();
            testCase.verifyClass(wd3, 'WDDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), wd3.trigonometricMoment(1), 'RelTol', 1E-10);
            wd5 = wn.toDirac5();
            testCase.verifyClass(wd5, 'WDDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), wd5.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(wn.trigonometricMoment(2), wd5.trigonometricMoment(2), 'RelTol', 1E-10);
            wd5 = wn.toDirac5(0.7);
            testCase.verifyClass(wd5, 'WDDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), wd5.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(wn.trigonometricMoment(2), wd5.trigonometricMoment(2), 'RelTol', 1E-10);
            vm = wn.toVM();
            testCase.verifyClass(vm, 'VMDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-10);
            
            %% test Gaussian conversion
            g = wn.toGaussian();
            testCase.verifyClass(g, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, wn.mu);
            testCase.verifyEqual(g.C, wn.sigma^2); % adjust this check, if a more sophisticated conversion is used
            
            %% test multiplication
            wnMul = wn.multiplyVM(wn);
            testCase.verifyClass(wnMul, 'WNDistribution');
            testCase.verifyEqual(wnMul.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wnMul.mu, wn.mu, 'RelTol', 1E-10);
            testCase.verifyLessThan(wnMul.sigma, wn.sigma);
            wnMul2 = wn.multiplyMomentBased(wn);
            testCase.verifyClass(wnMul2, 'WNDistribution');
            testCase.verifyEqual(wnMul2.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wnMul2.mu, wn.mu, 'RelTol', 1E-10);
            testCase.verifyLessThan(wnMul2.sigma, wn.sigma);
            wnMul3 = wn.multiplyTruncated(wn);
            testCase.verifyClass(wnMul3, 'WNDistribution');
            testCase.verifyEqual(wnMul3.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wnMul3.mu, wn.mu, 'RelTol', 1E-10);
            testCase.verifyLessThan(wnMul3.sigma, wn.sigma);
            
            %% test convolution
            wnConv = wn.convolve(wn);
            testCase.verifyClass(wnConv, 'WNDistribution');
            testCase.verifyEqual(wnConv.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wnConv.mu, wn.mu + wn.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wnConv.trigonometricMoment(1), wn.trigonometricMoment(1)*wn.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyGreaterThan(wnConv.sigma, wn.sigma);
            
            %% test sampling
            n = 10;
            s = wn.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
            
            %% test mle, mlenumerical
            s = wn.toDirac3().d;
            wnMle = wn.mleJensen(s);
            testCase.verifyClass(wnMle, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnMle.mu, 'RelTol', 1E-5);
            testCase.verifyEqual(wn.sigma, wnMle.sigma, 'RelTol', 0.2);
            
            wnMleNumerical = wn.mleNumerical(s);
            testCase.verifyClass(wnMleNumerical, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnMleNumerical.mu, 'RelTol', 1E-5);
            testCase.verifyEqual(wn.sigma, wnMleNumerical.sigma, 'RelTol', 0.2);
        end
        
        function testQuantization(testCase)
            wn = WNDistribution(2,0.7);
            n = 11;
            [s,w] = wn.sampleOptimalQuantization(n);
            testCase.verifyEqual(sum(w), 1, 'RelTol', 1E-10);
            testCase.verifySize(s, [1 n]);
            testCase.verifySize(w, [1 n]);
            wd = WDDistribution(s,w);
            testCase.verifyEqual(wd.trigonometricMoment(1), wn.trigonometricMoment(1), 'RelTol', 1E-2);
        end
    end
end