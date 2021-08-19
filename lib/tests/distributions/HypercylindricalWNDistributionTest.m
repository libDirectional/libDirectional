classdef HypercylindricalWNDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testL2C1(testCase)
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            testpoints = 2*pi*rand(3,100);
            
            mu = [2, 3, 4]';
            si1 = 1.9;
            si2 = 1.5;
            si3 = 0.7;
            rho12= 0.5;
            rho13= 0.3;
            rho23 = 0.4;
            C = [si1^2, si1*si2*rho12, si1*si3*rho13;
                si1*si2*rho12, si2^2, si2*si3*rho23;
                si1*si3*rho13, si2*si3*rho23, si3^2];
            pwn = HypercylindricalWNDistribution(mu, C, 1);
            
            %% test pdf
            function p = truePdf(x)
                p = 0;
                for k = -5:5
                    p = p + mvnpdf(x', mu' + [2*pi*k, 0, 0], C)';
                end
            end
            testCase.verifyEqual(pwn.pdf(testpoints), truePdf(testpoints), 'RelTol', 1E-5);
            
            %% test pdf with large uncertainty           
            %pwnLargeUncertainty = HypercylindricalWNDistribution(mu, 100*eye(3,3),2);
            %for t=1:6
            %    testCase.verifyEqual(pwnLargeUncertainty.pdf(mu), pwnLargeUncertainty.pdf(mu + [0,0,t]'), 'RelTol', 1E-10);
            %end

            %% test integral
            if enableExpensive            
                testCase.verifyEqual(pwn.integral, 1, 'RelTol', 1E-5);
            end
            
            %% test mean and maginalization
            mean = pwn.hybridMoment();
            hwn = pwn.marginalizeLinear();
            gauss = pwn.marginalizeCircular();
            testCase.verifyClass(hwn, 'HypertoroidalWNDistribution');
            testCase.verifyClass(gauss, 'GaussianDistribution');
            mom = hwn.trigonometricMoment(1);
            testCase.verifyEqual(mean(1), real(mom(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(mom(1)), 'RelTol', 1E-10);            
            testCase.verifyEqual(mean(2*pwn.boundD+1:end), gauss.mu, 'RelTol', 1E-10);
            
            %% test covariance
%             rng default %fix rng to get deterministic test
%             testCase.verifyEqual(pwn.covariance4D, pwn.covariance4DNumerical, 'RelTol', 0.3);
            
            %% test convolution
            pwnConv = pwn.convolve(pwn);
            testCase.verifyClass(pwnConv, 'HypercylindricalWNDistribution');
            correctMu = pwn.mu + pwn.mu;
            correctMu(1:pwn.boundD,:) = mod(correctMu(1:pwn.boundD,:), 2*pi);
            testCase.verifyEqual(pwnConv.mu, correctMu, 'RelTol', 1E-10);
            testCase.verifyEqual(pwnConv.C, pwn.C + pwn.C, 'RelTol', 1E-10);
            
            %% test sampling
            rng default
            n = 10;
            s = pwn.sample(n);
            testCase.verifyEqual(size(s,1), pwn.linD + pwn.boundD);
            testCase.verifyEqual(size(s,2), n);
            s = s(pwn.linD+1:end,:);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));                        
        end
        
        function testL1C2(testCase)
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            testpoints = 2*pi*rand(3,100);
            
            mu = [2, 3, 4]';
            si1 = 1.9;
            si2 = 1.5;
            si3 = 0.7;
            rho12= 0.5;
            rho13= 0.3;
            rho23 = 0.4;
            C = [si1^2, si1*si2*rho12, si1*si3*rho13;
                si1*si2*rho12, si2^2, si2*si3*rho23;
                si1*si3*rho13, si2*si3*rho23, si3^2];
            pwn = HypercylindricalWNDistribution(mu, C, 2);
            
            %% test pdf
            function p = truePdf(x)
                p = 0;
                for k = -5:5
                    for j = -5:5
                        p = p + mvnpdf(x', mu' + [2*pi*j, 2*pi*k, 0], C)';
                    end
                end
            end
            testCase.verifyEqual(pwn.pdf(testpoints), truePdf(testpoints), 'RelTol', 1E-5);
            
            %% test pdf with large uncertainty           
            %pwnLargeUncertainty = HypercylindricalWNDistribution(mu, 100*eye(3,3),2);
            %for t=1:6
            %    testCase.verifyEqual(pwnLargeUncertainty.pdf(mu), pwnLargeUncertainty.pdf(mu + [0,0,t]'), 'RelTol', 1E-10);
            %end

            %% test integral
            if enableExpensive            
                testCase.verifyEqual(pwn.integral, 1, 'RelTol', 1E-5);
            end
            
            %% test mean and maginalization
            mean = pwn.hybridMoment();
            hwn = pwn.marginalizeLinear();
            gauss = pwn.marginalizeCircular();
            testCase.verifyClass(hwn, 'HypertoroidalWNDistribution');
            testCase.verifyClass(gauss, 'GaussianDistribution');
            mom = hwn.trigonometricMoment(1);
            testCase.verifyEqual(mean(1), real(mom(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(mom(1)), 'RelTol', 1E-10);            
            testCase.verifyEqual(mean(3), real(mom(2)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(4), imag(mom(2)), 'RelTol', 1E-10);            
            testCase.verifyEqual(mean(2*pwn.boundD+1:end), gauss.mu, 'RelTol', 1E-10);
            %% test covariance
%             rng default %fix rng to get deterministic test
%             testCase.verifyEqual(pwn.covariance4D, pwn.covariance4DNumerical, 'RelTol', 0.3);
            
            %% test convolution
            pwnConv = pwn.convolve(pwn);
            testCase.verifyClass(pwnConv, 'HypercylindricalWNDistribution');
            correctMu = pwn.mu + pwn.mu;
            correctMu(1:pwn.boundD,:) = mod(correctMu(1:pwn.boundD,:), 2*pi);
            testCase.verifyEqual(pwnConv.mu, correctMu, 'RelTol', 1E-10);
            testCase.verifyEqual(pwnConv.C, pwn.C + pwn.C, 'RelTol', 1E-10);
            
            %% test sampling
            rng default
            n = 10;
            s = pwn.sample(n);
            testCase.verifyEqual(size(s,1), pwn.linD + pwn.boundD);
            testCase.verifyEqual(size(s,2), n);
            s = s(pwn.linD+1:end,:);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));                        
        end
        function testHybridMean(testCase)
            mu = [2, 3, 4]';
            si1 = 1.9;
            si2 = 1.5;
            si3 = 0.7;
            rho12= 0.5;
            rho13= 0.3;
            rho23 = 0.4;
            C = [si1^2, si1*si2*rho12, si1*si3*rho13;
                si1*si2*rho12, si2^2, si2*si3*rho23;
                si1*si3*rho13, si2*si3*rho23, si3^2];
            pwn = HypercylindricalWNDistribution(mu, C, 2);
            
            testCase.verifyEqual(pwn.hybridMean, mu);
        end
    end
end