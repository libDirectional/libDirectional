classdef WDDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWDDistribution(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            d = [1 2 3];
            w = [0.2 0.6 0.2];
            wd = WDDistribution(d, w);
            
            %% test pdf
            testCase.verifyError(@wd.pdf,'PDF:UNDEFINED');
            
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
            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            wdRew = wd.reweigh(f);
            testCase.verifyClass(wdRew, 'WDDistribution');
            testCase.verifyEqual(wdRew.d, wd.d);
            testCase.verifyEqual(wdRew.w, wd.w);
            
            f = @(x) x(1,:);
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
            testCase.applyFixture(SuppressedWarningsFixture('ENTROPY:DISCRETE'));
            e = wd.entropy();
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
        
        function testToContinuousVoronoi(testCase)
            d=[0.1,5,6];
            w=[0.3,0.1,0.6];
            wd=WDDistribution(d,w);
            cd=wd.toContinuousVoronoi;
            testCase.verifyGreaterThan(cd.pdf(0.5*(d(2)+d(1))-0.1),cd.pdf(0.5*(d(2)+d(1))+0.1));
            testCase.verifyEqual(cd.pdf(2*pi-0.01),cd.pdf(0.01), 'RelTol', 1E-10);
            
            d=[0.5,5,6];
            w=[0.3,0.4,0.3];
            wd=WDDistribution(d,w);
            cd=wd.toContinuousVoronoi;
            testCase.verifyGreaterThan(cd.pdf(0.5*(d(2)+d(1))+0.1),cd.pdf(0.5*(d(2)+d(1))-0.1));
            testCase.verifyEqual(cd.pdf(2*pi-0.01),cd.pdf(0.01), 'RelTol', 1E-10);
            
            wd2=WDDistribution(d([2,3,1]),w([2,3,1]));
            cd2=wd2.toContinuousVoronoi;
            testCase.verifyEqual(cd2.totalVariationDistanceNumerical(cd),0, 'RelTol', 1E-10);
            
            wd3=WDDistribution([d d],[w/2 w/2]);
            cd3=wd3.toContinuousVoronoi;
            testCase.verifyEqual(cd3.totalVariationDistanceNumerical(cd),0, 'RelTol', 1E-10);
        end
        
        function testl2distanceCdf(testCase)
            d=[0,1,2,3,3.5,6];
            perm=[2,5,6,1,4,3];
            w=[0.3,0.1,0.15,0.2,0.1,0.05];            
            wd=WDDistribution(d(perm),w);
            startingPoint=1.5;
            
            vm=VMDistribution(0,10);
            fdId=FourierDistribution.fromDistribution(vm,101,'identity');
            fdSqrt=FourierDistribution.fromDistribution(vm,101,'sqrt');
            
            cvm=integral(@(x)(wd.cdf(x,startingPoint)-vm.cdf(x,startingPoint)).^2,0,2*pi);
            
            testCase.verifyEqual(wd.l2distanceCdfNumerical(vm,startingPoint),cvm,'RelTol',1E-6);
            testCase.verifyEqual(wd.l2distanceCdfNumerical(fdId,startingPoint),cvm,'RelTol',1E-6);
            testCase.verifyEqual(wd.l2distanceCdfNumerical(fdSqrt,startingPoint),cvm,'RelTol',1E-6);
            
            testCase.verifyEqual(vm.l2distanceCdfNumerical(wd,startingPoint),cvm,'RelTol',1E-6);
            testCase.verifyEqual(fdId.l2distanceCdfNumerical(wd,startingPoint),cvm,'RelTol',1E-6);
            testCase.verifyEqual(fdSqrt.l2distanceCdfNumerical(wd,startingPoint),cvm,'RelTol',1E-6);
            
            % Test wd-wd
            d2=rand(1,1000)*45/23*pi;
            perm2=randperm(numel(d2));
            w2=rand(size(d2));
            w2=w2/norm(w2,1);
            wd2=WDDistribution(d2(perm2),w2);
            cvm=integral(@(x)(wd.cdf(x,startingPoint)-wd2.cdf(x,startingPoint)).^2,0,2*pi);
            % AbsTol of 1E-4 is not that good, maybe a more thorough
            % comparison is necessary.
            testCase.verifyEqual(wd.l2distanceCdfNumerical(wd2,startingPoint),cvm,'AbsTol',1E-4);
        end
        function testFromDistribution(testCase)
            rng default
            wn = WNDistribution(1,1);
            wd = WDDistribution.fromDistribution(wn,100000);
            testCase.verifyEqual(wd.meanDirection, wn.meanDirection, 'AbsTol',0.001);
        end
    end
end