classdef SE2PWNDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testSE2PWNDistribution(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end           
            mu = [2, 3, 4]';
            si1 = 0.9;
            si2 = 1.5;
            si3 = 1.7;
            rho12= 0.5;
            rho13= 0.3;
            rho23 = 0.4;
            C = [si1^2, si1*si2*rho12, si1*si3*rho13;
                si1*si2*rho12, si2^2, si2*si3*rho23;
                si1*si3*rho13, si2*si3*rho23, si3^2];
            pwn = SE2PWNDistribution(mu,C);
            
            % test pdf
            function p = truePdf(x)
                p = 0;
                for k = -20:20
                    p = p + mvnpdf(x', mu' + [2*pi*k, 0 0], C);
                end
            end
            testCase.verifyEqual(pwn.pdf(mu), truePdf(mu), 'RelTol', 1E-10);
            testCase.verifyEqual(pwn.pdf(mu-1), truePdf(mu-1), 'RelTol', 1E-10);
            testCase.verifyEqual(pwn.pdf(mu+2), truePdf(mu+2), 'RelTol', 1E-10);
            x = rand(3,20);
            testCase.verifyEqual(pwn.pdf(x), truePdf(x)', 'RelTol', 1E-10);
            
            %% test pdf with large uncertainty           
            pwnLargeUncertainty = SE2PWNDistribution(mu, 100*eye(3,3));
            for t=1:6
                testCase.verifyEqual(pwnLargeUncertainty.pdf(mu), pwnLargeUncertainty.pdf(mu + [t,0,0]'), 'RelTol', 1E-10);
            end

            % test integral
            if enableExpensive
                testCase.verifyEqual(pwn.integral, 1, 'RelTol', 1E-5);
            end
            
            % test mean and maginalization
            mean = pwn.mean4D();
            wn = pwn.marginalizeLinear();
            gauss = pwn.marginalizeCircular();
            testCase.verifyEqual(mean(1), real(wn.trigonometricMoment(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(wn.trigonometricMoment(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(3:4), gauss.mu, 'RelTol', 1E-10);
            
            % test covariance
            rng default %fix rng to get deterministic test
            testCase.verifyEqual(pwn.covariance4D, pwn.covariance4DNumerical, 'RelTol', 0.3);
            
            % test convolution
            pwnConv = pwn.convolve(pwn);
            testCase.verifyClass(pwnConv, 'SE2PWNDistribution');
            testCase.verifyEqual(pwnConv.mu, pwn.mu + pwn.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(pwnConv.C, pwn.C + pwn.C, 'RelTol', 1E-10);
            
            % test toGaussian
            g = pwn.toGaussian();
            testCase.verifyClass(g, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, pwn.mu);
            testCase.verifyEqual(g.C, pwn.C);            
            
            % test sampling
            rng default
            n = 10;
            s = pwn.sample(n);
            testCase.verifyEqual(size(s,1), 3);
            testCase.verifyEqual(size(s,2), n);
            s = s(1,:);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
                        
            %todo proper test
            pwd = pwn.toSE2PWD();
            pwn2 = SE2PWNDistribution.fromSamples(pwd.d); %todo do not ignore weights
            testCase.verifyClass(pwn2, 'SE2PWNDistribution');
        end
    end
end