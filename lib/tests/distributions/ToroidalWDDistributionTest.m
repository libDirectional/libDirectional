classdef ToroidalWDDistributionTest< matlab.unittest.TestCase
    properties
    end
    
    methods (Test)
        function testToroidalWDDistribution(testCase)
            d = [0.5 3 4 6 6;
                 2.5 2 5 3 0];
            w = [0.1 0.1 0.1 0.1 0.6];
            twd = ToroidalWDDistribution(d,w);
            
            %% integral
            testCase.verifyError(@twd.integralNumerical, 'PDF:UNDEFINED');
            testCase.verifyEqual(twd.integral, 1, 'RelTol', 1E-10)
            testCase.verifyEqual(twd.integral([0;0], [pi;2*pi]) + twd.integral([pi;0], [2*pi; 2*pi]), 1, 'RelTol', 1E-10)
            testCase.verifyEqual(twd.integral([0;0], [2*pi;pi]) + twd.integral([0;pi], [2*pi; 2*pi]), 1, 'RelTol', 1E-10)
            testCase.verifyEqual(twd.integral([0;0], [2*pi;pi]) + twd.integral([0;pi], [2*pi; 0]), 0, 'RelTol', 1E-10)
            testCase.verifyEqual(twd.integral([0;0], [pi;pi]) + twd.integral([0; pi], [pi;2*pi]) + twd.integral([pi;0], [2*pi;pi]) + twd.integral([pi;pi],[2*pi;2*pi]), 1, 'RelTol', 1E-10)
            testCase.verifyEqual(twd.integral([0;0], [d(1,3);2*pi]) + twd.integral([d(1,3);0], [2*pi; 2*pi]), 1, 'RelTol', 1E-10) % check integration border on Dirac
            testCase.verifyEqual(twd.integral([0;0], [2*pi;d(2,3)]) + twd.integral([0;d(2,3)], [2*pi; 2*pi]), 1, 'RelTol', 1E-10) % check integration border on Dirac
            for i=1:size(d,2)
                testCase.verifyEqual(twd.integral([d(1,i);d(2,i)],[d(1,i)+0.1; d(2,i)+0.1]), w(i));
            end
            
            %% sanity check
            testCase.verifyClass(twd, 'ToroidalWDDistribution');
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
            
            %% test mean4D
            meanVal = twd.mean4D();
            testCase.verifyEqual(meanVal(1), real(m1), 'RelTol', 1E-10);
            testCase.verifyEqual(meanVal(2), imag(m1), 'RelTol', 1E-10);
            testCase.verifyEqual(meanVal(3), real(m2), 'RelTol', 1E-10);
            testCase.verifyEqual(meanVal(4), imag(m2), 'RelTol', 1E-10);
            
            %% test covariance4D
            rng default %fix rng to get deterministic test
            C = twd.covariance4D();
            testCase.verifyEqual(C,C', 'RelTol', 1E-10); 
            testCase.verifyGreaterThan(eig(C),[0; 0; 0; 0]); 
                                    
            %% test sampling
            nSamples = 5;
            s = twd.sample(nSamples);
            testCase.verifyEqual(size(s,1),2);
            testCase.verifyEqual(size(s,2),nSamples);
            testCase.verifyEqual(s, mod(s,2*pi));
            
            %% test correlation coefficients
            pm = twd.angularProductMoment(1);
            testCase.verifyGreaterThan(pm,-1);
            testCase.verifyLessThan(pm,1);
            rhoc = twd.circularCorrelationJammalamadaka();
            testCase.verifyGreaterThan(rhoc, -1);
            testCase.verifyLessThan(rhoc, 1);
            
            %% test conversions based on angular moment matching
            twn = twd.toToroidalWN();
            testCase.verifyClass(twn, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(twd.circularCorrelationJammalamadaka(), twn.circularCorrelationJammalamadaka(), 'RelTol', 1E-10);
            
            twn2 = twd.toToroidalWNjammalamadaka();
            testCase.verifyClass(twn2, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn2.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(real(twd.angularProductMoment(1)), twn2.angularProductMomentNumerical(1), 'RelTol', 1E-10);
            
            twn3 = twd.toToroidalWNjupp();
            testCase.verifyClass(twn3, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn3.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(twd.circularCorrelationJupp(), twn3.circularCorrelationJupp(), 'RelTol', 1E-10);
            
            twn4 = twd.toToroidalWNmixedMLE();
            testCase.verifyClass(twn4, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn4.trigonometricMoment(1), 'RelTol', 1E-10);
            %todo check maximum likelihood property?
            
            twn5 = twd.toToroidalWNunwrappingEM();
            testCase.verifyClass(twn5, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn4.trigonometricMoment(1), 'RelTol', 1E-10);
            %todo check correlation somehow
            
            twn6 = twd.toToroidalWNcovariance();
            testCase.verifyClass(twn6, 'ToroidalWNDistribution');
            testCase.verifyEqual(twd.trigonometricMoment(1), twn2.trigonometricMoment(1), 'RelTol', 1E-10);
            Ctwd = twd.covariance4D();
            C1 = twn.covariance4D();
            C2 = twn2.covariance4D();
            C3 = twn3.covariance4D();
            C4 = twn4.covariance4D();
            C5 = twn5.covariance4D();
            C6 = twn6.covariance4D();
            frob = @(x) norm(x(:)); % Frobenius norm
            testCase.verifyLessThanOrEqual( frob(Ctwd(1:2,3:4) - C6(1:2,3:4)), frob(Ctwd(1:2,3:4) - C1(1:2,3:4)));
            testCase.verifyLessThanOrEqual( frob(Ctwd(1:2,3:4) - C6(1:2,3:4)), frob(Ctwd(1:2,3:4) - C2(1:2,3:4)));
            testCase.verifyLessThanOrEqual( frob(Ctwd(1:2,3:4) - C6(1:2,3:4)), frob(Ctwd(1:2,3:4) - C3(1:2,3:4)));            
            testCase.verifyLessThanOrEqual( frob(Ctwd(1:2,3:4) - C6(1:2,3:4)), frob(Ctwd(1:2,3:4) - C4(1:2,3:4)));            
            testCase.verifyLessThanOrEqual( frob(Ctwd(1:2,3:4) - C6(1:2,3:4)), frob(Ctwd(1:2,3:4) - C5(1:2,3:4)));            
            
            %% test errors
            testCase.verifyError(@() twd.angularProductMomentNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() twd.covariance4DNumerical(), 'PDF:UNDEFINED');
                       
            %% test getMarginal
            wd1 = twd.marginalizeTo1D(1);
            wd2 = twd.marginalizeTo1D(2);
            testCase.verifyEqual(twd.w, wd1.w);
            testCase.verifyEqual(twd.w, wd2.w);
            testCase.verifyEqual(twd.d(1,:), wd1.d);
            testCase.verifyEqual(twd.d(2,:), wd2.d);         
            
            %% test apply function 
            same = twd.applyFunction(@(x) x);
            testCase.verifyEqual(same.trigonometricMoment(1), twd.trigonometricMoment(1));
            shiftOffset = [1.4; -0.3];
            shifted = twd.applyFunction(@(x) x + shiftOffset);
            testCase.verifyEqual(shifted.trigonometricMoment(1), twd.trigonometricMoment(1) .* exp(1i*shiftOffset), 'RelTol', 1E-10);
            
            %% test reweigh
            f = @(x) sum(x)==3; %only dirac with sum 3 gets weight
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'ToroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            testCase.verifyEqual(twdRew.w, double(f(twd.d)));
            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'ToroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            testCase.verifyEqual(twdRew.w, twd.w);
            
            f = @(x) x(1,:);
            twdRew = twd.reweigh(f);
            testCase.verifyClass(twdRew, 'ToroidalWDDistribution');
            testCase.verifyEqual(twdRew.d, twd.d);
            wNew = twd.d(1,:).*twd.w;
            testCase.verifyEqual(twdRew.w, wNew/sum(wNew));
        end
        
        function testShift(testCase)
            d = [0.5,3,4,6,6;
                 2,2,5,3,0];
            w = [0.1 0.1 0.1 0.1 0.6];
            twd = ToroidalWDDistribution(d,w);
            s = [-3; 6];
            twdShifted = twd.shift(s);
            testCase.verifyClass(twdShifted, 'ToroidalWDDistribution');
            testCase.verifyEqual(twd.w, twdShifted.w);
            testCase.verifyEqual(twd.d, mod(twdShifted.d - repmat(s,1,size(d,2)),2*pi), 'RelTol', 1E-10);
        end        
    end
end