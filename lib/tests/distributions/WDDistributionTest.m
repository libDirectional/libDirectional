classdef WDDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWDDistribution(testCase)
            d = [1 2 3];
            w = [0.2 0.6 0.2];
            wd = WDDistribution(d, w);
            
            %% test pdf
            testCase.verifyWarning(@wd.pdf,'PDF:UNDEFINED');
            
            %% test cdf
            x = 0:6;
            fx = arrayfun(@(y) wd.integral(0, y), x);
            testCase.verifyEqual(wd.cdf(x), fx, 'RelTol', 1E-10);
            startingPoint = 4;
            function y = cdf(x)
                if x>=startingPoint
                    y = wd.integral(startingPoint, x);
                else
                    y = wd.integral(startingPoint,2*pi) + wd.integral(0, x);
                end
            end
            fx = arrayfun(@cdf, x);
            testCase.verifyEqual(wd.cdf(x,startingPoint), fx, 'RelTol', 1E-10);
            
            %% test integral
            testCase.verifyEqual(wd.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0.5,1.5), w(1), 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0.5,2.5), w(1)+w(2), 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0.5,2*pi+1.5), w(1)+1, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0.5,4*pi+1.5), w(1)+2, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0,4*pi), 2, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0.5-2*pi,2.5), w(1)+w(2)+1, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(0,10*pi), 5, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(-pi,pi), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(-8*pi,6*pi), 7, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(3.1,2*pi+0.1), 0, 'RelTol', 1E-10);
            testCase.verifyEqual(wd.integral(2.5,2*pi+0.1), w(3), 'RelTol', 1E-10);
            
            %% test angular moments
            testCase.verifyEqual(wd.trigonometricMoment(1), sum(w.*exp(1i*d)), 'RelTol', 1E-10);
            testCase.verifyEqual(wd.circularMean, 2);
            
            %% test conversions based on angular moment matching
            wnMoment = wd.toWN();
            testCase.verifyClass(wnMoment, 'WNDistribution');
            testCase.verifyEqual(wd.trigonometricMoment(1), wnMoment.trigonometricMoment(1), 'RelTol', 1E-10);
            vm = wd.toVM();
            testCase.verifyClass(vm, 'VMDistribution');
            testCase.verifyEqual(wd.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-10);
            
            %% apply function
            same = wd.applyFunction(@(x) x);
            testCase.verifyEqual(same.trigonometricMoment(1), wd.trigonometricMoment(1));
            shiftOffset = 1.4;
            shifted = wd.applyFunction(@(x) x + shiftOffset);
            testCase.verifyClass(shifted, 'WDDistribution');
            shiftedMu = shifted.toWN().mu;
            originalMu = wnMoment.mu;
            testCase.verifyEqual(shiftOffset, abs(shiftedMu - originalMu), 'RelTol', 1E-10);
            
            %% test sampling
            n = 10;
            s = wd.sample(10);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
            
            %% Kuiper's test
            t = wd.kuiperTest(wnMoment);
            testCase.verifyGreaterThan(t,0);
            
            %% toPWC
            n = 10;
            pwc = wd.toPWC(n);
            testCase.verifyClass(pwc, 'PWCDistribution');
            testCase.verifySize(pwc.w, [1,n]);
            
            %% unwrapping EM
            wnUnwrap = wd.toWNunwrappingEM();
            testCase.verifyClass(wnUnwrap, 'WNDistribution');
            testCase.verifyEqual(wnUnwrap.mu, wnMoment.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wnUnwrap.sigma, wnMoment.sigma, 'RelTol', 0.3);
            
            %% reweigh
            f = @(x) x==2; %only dirac at 2 gets weight
            wdRew = wd.reweigh(f);
            testCase.verifyClass(wdRew, 'WDDistribution');
            testCase.verifyEqual(wdRew.d, wd.d);
            testCase.verifyEqual(wdRew.w, double(f(wd.d)));
            
            f = @(x) 2; %does not change anything because of renormalization
            wdRew = wd.reweigh(f);
            testCase.verifyClass(wdRew, 'WDDistribution');
            testCase.verifyEqual(wdRew.d, wd.d);
            testCase.verifyEqual(wdRew.w, wd.w);
            
            f = @(x) x;
            wdRew = wd.reweigh(f);
            testCase.verifyClass(wdRew, 'WDDistribution');
            testCase.verifyEqual(wdRew.d, wd.d);
            wNew = wd.d.*wd.w;
            testCase.verifyEqual(wdRew.w, wNew/sum(wNew));
            
            %% convolution
            wd2 = WDDistribution([0, 4, 5, 6], [1, 3, 4, 2]);
            wdConv = wd.convolve(wd2);
            testCase.verifyClass(wdConv, 'WDDistribution');
            testCase.verifyEqual(length(wdConv.d), length(wd.d)*length(wd2.d));
            
            %% entropy
            testCase.verifyWarning(@wd.entropy, 'ENTROPY:DISCRETE');
            %disable warning for the next function call
            warning('off','ENTROPY:DISCRETE');
            e = wd.entropy();
            %reenable warning
            warning('on','ENTROPY:DISCRETE');
            testCase.verifyGreaterThanOrEqual(e, 0);
        end
        
        function testShift(testCase)
            d = [0.5,3,4,6,6];
            w = [0.1 0.1 0.1 0.1 0.6];
            twd = WDDistribution(d,w);
            s = 6;
            twdShifted = twd.shift(s);
            testCase.verifyClass(twdShifted, 'WDDistribution');
            testCase.verifyEqual(twd.w, twdShifted.w);
            testCase.verifyEqual(twd.d, mod(twdShifted.d - repmat(s,1,size(d,2)),2*pi), 'RelTol', 1E-10);
        end        
    end
end