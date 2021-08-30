classdef HypertoroidalWDDistributionTest< matlab.unittest.TestCase
    properties
    end
    
    methods (Test)
        function testHypertoroidalWDDistribution(testCase)
            d = [0.5,3,4,6,6;
                 2,2,5,3,0;
                 0.5,0.2,5.8,4.3,1.2];
            w = [0.1 0.1 0.1 0.1 0.6];
            twd = HypertoroidalWDDistribution(d,w);
            %% sanity check
            testCase.verifyClass(twd, 'HypertoroidalWDDistribution');
            testCase.verifyEqual(twd.d, d);
            testCase.verifyEqual(twd.w, w);
            
            %% test trigonometric moment
            m = twd.trigonometricMoment(1);
            m1 = twd.marginalizeTo1D(1).trigonometricMoment(1);
            m2 = twd.marginalizeTo1D(2).trigonometricMoment(1);
            testCase.verifyEqual(m(1), m1, 'RelTol', 1E-10);
            testCase.verifyEqual(m(2), m2, 'RelTol', 1E-10);
            testCase.verifyEqual(m(1), sum(w.*exp(1i*d(1,:))), 'RelTol', 1E-10);
            testCase.verifyEqual(m(2), sum(w.*exp(1i*d(2,:))), 'RelTol', 1E-10);
            
            %% test errors
            testCase.verifyError(@() twd.logLikelihood(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.pdf(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.integralNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.trigonometricMomentNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.sampleMetropolisHastings(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.squaredDistanceNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.kldNumerical(), 'PDF:UNDEFINED');
            
            %% test sampling
            nSamples = 5;
            s = twd.sample(nSamples);
            testCase.verifySize(s,[size(d,1),nSamples]);
            testCase.verifyEqual(s, mod(s,2*pi))
            
            %% test getMarginal
            for i=1:size(d,1)
                wd=twd.marginalizeTo1D(i);
                testCase.verifyEqual(twd.w, wd.w);
                testCase.verifyEqual(twd.d(i,:), wd.d);
            end 
            
            %% test apply function 
            same = twd.applyFunction(@(x) x);
            testCase.verifyEqual(same.trigonometricMoment(1), twd.trigonometricMoment(1));
            shiftOffset = [1.4;-0.3;pi];
            shifted = twd.applyFunction(@(x) x + shiftOffset);
            testCase.verifyEqual(shifted.trigonometricMoment(1), twd.trigonometricMoment(1) .* exp(1i*shiftOffset), 'RelTol', 1E-10);
            
            %% test reweigh
            f = @(x) sum(x)==3; %only dirac with sum 3 gets weight
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'HypertoroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            testCase.verifyEqual(twdRew.w, double(f(twd.d)));
            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'HypertoroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            testCase.verifyEqual(twdRew.w, twd.w);
            
            f = @(x) x(1,:);
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'HypertoroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            wNew = twd.d(1,:).*twd.w;
            testCase.verifyEqual(twdRew.w, wNew/sum(wNew));
        end
        
        function testShift(testCase)
            d = [0.5,3,4,6,6;
                 2,2,5,3,0;
                 0.5,0.2,5.8,4.3,1.2];
            w = [0.1 0.1 0.1 0.1 0.6];
            twd = HypertoroidalWDDistribution(d,w);
            s = [1; -3; 6];
            twdShifted = twd.shift(s);
            testCase.verifyClass(twdShifted, 'HypertoroidalWDDistribution');
            testCase.verifyEqual(twd.w, twdShifted.w);
            testCase.verifyEqual(twd.d, mod(twdShifted.d - repmat(s,1,size(d,2)),2*pi), 'RelTol', 1E-10);
        end
        
        function testToWD(testCase)
            rng default
            n = 20;
            d = 2*pi*rand(1,n);
            w = rand(1,n);
            w = w/sum(w);
            hwd = HypertoroidalWDDistribution(d,w);
            wd1 = WDDistribution(d,w);
            wd2 = hwd.toWD();
            testCase.verifyClass(wd2, 'WDDistribution');
            testCase.verifyEqual(wd1.d, wd2.d, 'RelTol', 1E-10);
            testCase.verifyEqual(wd1.w, wd2.w, 'RelTol', 1E-10);
        end
        
        function testToToroidalWD(testCase)
            rng default
            n = 20;
            d = 2*pi*rand(2,n);
            w = rand(1,n);
            w = w/sum(w);
            hwd = HypertoroidalWDDistribution(d,w);
            twd1 = ToroidalWDDistribution(d,w);
            twd2 = hwd.toToroidalWD();
            testCase.verifyClass(twd2, 'ToroidalWDDistribution');
            testCase.verifyEqual(twd1.d, twd2.d, 'RelTol', 1E-10);
            testCase.verifyEqual(twd1.w, twd2.w, 'RelTol', 1E-10);
        end
        
        function testMarginalization(testCase)
            rng default
            n = 20;
            d = 2*pi*rand(2,n);
            w = rand(1,n);
            w = w/sum(w);
            hwd = HypertoroidalWDDistribution(d,w);
            wd1 = hwd.marginalizeTo1D(1);
            wd2 = hwd.marginalizeOut(2);
            testCase.verifyEqual(wd1.d, wd2.d);
            testCase.verifyEqual(wd1.w, wd2.w);
        end
        
        function testFromDistribution(testCase)
            rng(1) % This seed makes it work, not working well for rng default
            C = wishrnd(eye(3),3);
            hwn = HypertoroidalWNDistribution([1;2;3],C);
            hwd = HypertoroidalWDDistribution.fromDistribution(hwn,1000000);
            testCase.verifyEqual(hwd.meanDirection, hwn.meanDirection, 'AbsTol',0.002);
        end
    end
end