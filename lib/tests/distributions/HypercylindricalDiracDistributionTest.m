classdef HypercylindricalDiracDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testlinD2boundD1(testCase)
            d = [ 1 2 3 4 5 6;
                  2 4 0 0.5 1 1;
                  0 10 20 30 40 50];
            w = [ 1 2 3 1 2 3];
            w = w/sum(w);
            pwd = HypercylindricalDiracDistribution(1, d, w);
            
            %% test mean and maginalization
            mean = pwd.hybridMoment();
            wd = pwd.marginalizeLinear();
            testCase.verifyEqual(mean(1), real(wd.trigonometricMoment(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(wd.trigonometricMoment(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(3), sum(w.*d(2,:),2), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(4), sum(w.*d(3,:),2), 'RelTol', 1E-10);

            %% test covariance
            %C = pwd.covariance4D;
            %testCase.verifyEqual(size(C), [4 4]);
            
            Clin = pwd.linearCovariance();
            testCase.verifySize(Clin, [pwd.linD pwd.linD]);
            
            %% test apply function 
            same = pwd.applyFunction(@(x) x);
            testCase.verifyEqual(pwd.d, same.d);
            testCase.verifyEqual(pwd.w, same.w);
            testCase.verifyEqual(pwd.linD, same.linD);
            testCase.verifyEqual(pwd.boundD, same.boundD);
            shiftOffset = [1.4; -0.3; 1];
            shifted = pwd.applyFunction(@(x) x + shiftOffset);
            testCase.verifyEqual(shifted.marginalizeLinear.trigonometricMoment(1), pwd.marginalizeLinear.trigonometricMoment(1) .* exp(1i*shiftOffset(1)), 'RelTol', 1E-10);
            
            %% reweigh
            f = @(x) sum(x)==3; %only dirac with sum 3 gets weight
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'HypercylindricalDiracDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            testCase.verifyEqual(pwdRew.w, double(f(pwd.d)));
            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'HypercylindricalDiracDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            testCase.verifyEqual(pwdRew.w, pwd.w);
            
            f = @(x) x(1,:);
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'HypercylindricalDiracDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            wNew = pwd.d(1,:).*pwd.w;
            testCase.verifyEqual(pwdRew.w, wNew/sum(wNew));
                        
            %% test sampling
            rng default
            n = 10;
            s = pwd.sample(n);
            testCase.verifyEqual(size(s,1), 3);
            testCase.verifyEqual(size(s,2), n);
            s = s(1,:);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
        
        function testFromDistribution(testCase)
            rng default % Could fail randomly otherwise
            C = wishrnd(eye(4),4);
            hwn = HypercylindricalWNDistribution([1;2;3;4],C,2);
            hddist = HypercylindricalDiracDistribution.fromDistribution(hwn,1000000);
            testCase.verifyEqual(hddist.hybridMoment, hwn.hybridMoment, 'AbsTol', 0.002);
            testCase.verifyEqual(hddist.hybridMean, hwn.hybridMean, 'AbsTol', 0.02);
        end
    end
end