classdef SE2PWDDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSE2PWDDistribution(testCase)
            d = [ 1 2 3 4 5 6;
                  2 4 0 0.5 1 1;
                  0 10 20 30 40 50];
            w = [ 1 2 3 1 2 3];
            w = w/sum(w);
            pwd = SE2PWDDistribution(d,w);
            
            % test mean and maginalization
            mean = pwd.mean4D();
            wd = pwd.marginalizeLinear();
            testCase.verifyEqual(mean(1), real(wd.trigonometricMoment(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(wd.trigonometricMoment(1)), 'RelTol', 1E-10);

            % test error
            testCase.verifyError(@() pwd.pdf(), 'PDF:UNDEFINED');
            
            % test covariance
            C = pwd.covariance4D;
            testCase.verifyEqual(size(C), [4 4]);
            
            % test apply function 
            same = pwd.applyFunction(@(x) x);
            testCase.verifyEqual(pwd.d, same.d);
            testCase.verifyEqual(pwd.w, same.w);
            shiftOffset = [1.4; -0.3; 1];
            shifted = pwd.applyFunction(@(x) x + shiftOffset);
            testCase.verifyEqual(shifted.marginalizeLinear.trigonometricMoment(1), pwd.marginalizeLinear.trigonometricMoment(1) .* exp(1i*shiftOffset(1)), 'RelTol', 1E-10);
            
            % reweigh
            f = @(x) sum(x)==3; %only dirac with sum 3 gets weight
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'SE2PWDDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            testCase.verifyEqual(pwdRew.w, double(f(pwd.d)));
            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'SE2PWDDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            testCase.verifyEqual(pwdRew.w, pwd.w);
            
            f = @(x) x(1,:);
            pwdRew = pwd.reweigh(f);
            testCase.verifyClass(pwdRew, 'SE2PWDDistribution');
            testCase.verifyEqual(pwdRew.d, pwd.d);
            wNew = pwd.d(1,:).*pwd.w;
            testCase.verifyEqual(pwdRew.w, wNew/sum(wNew));
                        
            % test sampling
            rng default
            n = 10;
            s = pwd.sample(n);
            testCase.verifyEqual(size(s,1), 3);
            testCase.verifyEqual(size(s,2), n);
            s = s(1,:);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
            
            % test initialization with equal weights
            pwd2 = SE2PWDDistribution(d);
            testCase.verifyEqual(pwd2.w, ones(1,size(d,2))/size(d,2), 'RelTol', 1E-10)
        end
    end
end