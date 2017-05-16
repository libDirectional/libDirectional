classdef WCDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWCDistribution(testCase)
            mu = 3;
            gamma = 1.5;
            wc = WCDistribution(mu,gamma);
            %% test pdf
            function p = truePdf(x)
                sum = 0;
                for k = -20000:20000
                    sum = sum + gamma/(pi * (gamma^2 + (x-mu +2*pi*k)^2));
                end
                p = sum;
            end
            testCase.verifyEqual(wc.pdf(mu), truePdf(mu), 'RelTol', 1E-5);
            testCase.verifyEqual(wc.pdf(mu-1), truePdf(mu-1), 'RelTol', 1E-5);
            testCase.verifyEqual(wc.pdf(mu+2), truePdf(mu+2), 'RelTol', 1E-3);
            
            %% test shift
            delta = 2;
            wcShifted = wc.shift(delta);
            testCase.verifyEqual(wc.pdf((1:10)-delta), wcShifted.pdf(1:10), 'RelTol', 1E-5);
            
            %% test cdf
            testCase.verifyEqual(wc.cdf(mu),wc.cdfNumerical(mu), 'RelTol', 1E-8);
            testCase.verifyEqual(wc.cdf(mu+1.3),wc.cdfNumerical(mu+1.3), 'RelTol', 1E-8);
            testCase.verifyEqual(wc.cdf(-10:0.1:10),wc.cdfNumerical(-10:0.1:10), 'RelTol', 1E-8);
            testCase.verifyEqual(wc.cdf(-10:10, 2),wc.cdfNumerical(-10:10, 2), 'RelTol', 1E-8);

            %% test integral
            testCase.verifyEqual(wc.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wc.integral(0,pi)+wc.integral(pi,2*pi), 1, 'RelTol', 1E-10);    
            
            %% test angular moments
            testCase.verifyEqual(wc.trigonometricMoment(1), wc.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(wc.trigonometricMoment(2), wc.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(wc.trigonometricMoment(3), wc.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            
            testCase.verifyEqual(wc.circularMean, mu);
            %% test fromMoment
            wcFromMoment = WCDistribution.fromMoment(wc.trigonometricMoment(1));
            testCase.verifyClass(wcFromMoment, 'WCDistribution');
            testCase.verifyEqual(wc.mu, wcFromMoment.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wc.gamma, wcFromMoment.gamma, 'RelTol', 1E-10);
            
            %% test conversions based on angular moment matching
            wd2 = wc.toDirac2();
            testCase.verifyClass(wd2, 'WDDistribution');
            testCase.verifyEqual(wc.trigonometricMoment(1), wd2.trigonometricMoment(1), 'RelTol', 1E-10);
            wd3 = wc.toDirac3();
            testCase.verifyClass(wd3, 'WDDistribution');
            testCase.verifyEqual(wc.trigonometricMoment(1), wd3.trigonometricMoment(1), 'RelTol', 1E-10);
            wd5 = wc.toDirac5();
            testCase.verifyClass(wd5, 'WDDistribution');
            testCase.verifyEqual(wc.trigonometricMoment(1), wd5.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(wc.trigonometricMoment(2), wd5.trigonometricMoment(2), 'RelTol', 1E-10);
            
            %% test convolution
            wcConv = wc.convolve(wc);
            testCase.verifyClass(wcConv, 'WCDistribution');
            testCase.verifyEqual(wcConv.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wcConv.mu, wc.mu + wc.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wcConv.trigonometricMoment(1), wc.trigonometricMoment(1)*wc.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyGreaterThan(wcConv.gamma, wc.gamma);
            
            %% test entropy
            testCase.verifyEqual(wc.entropy(), wc.entropyNumerical(), 'RelTol', 1E-10);
            
            %% test sampling
            n = 10;
            s = wc.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
            
            %% test periodicity
            testCase.verifyEqual(wc.pdf(linspace(-2*pi,0,100)),wc.pdf(linspace(0,2*pi,100)),'RelTol',1E-10);
        end
    end
end