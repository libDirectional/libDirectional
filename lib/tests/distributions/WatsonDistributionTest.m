classdef WatsonDistributionTest< matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testWatsonDistribution(testCase)
            mu = [1, 2, 3]';
            mu = mu/norm(mu);
            kappa = 2;
            w = WatsonDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(w, 'WatsonDistribution');
            testCase.verifyEqual(w.mu, mu);
            testCase.verifyEqual(w.kappa, kappa);
            testCase.verifyEqual(w.dim, length(mu));

            %% test mode
            % respect antipodal symmetry
            modeNumerical = w.modeNumerical;
            difference = min(abs(modeNumerical-w.mode), abs(modeNumerical+w.mode));
            testCase.verifyEqual(difference, zeros(w.dim,1), 'AbsTol',1E-5)

            %% test integral
            testCase.verifyEqual(w.integral(), 1, 'RelTol', 1E-5);
            
            %% test pdf
            kappa = 0.1;
            x = [0:10; 1:11; 5:15];
            x = bsxfun(@rdivide, x, sqrt(sum(x.^2,1))); %normalize
            wZero = WatsonDistribution([0,0,1]', kappa);
            testCase.verifyEqual(wZero.pdf(x), wZero.C*exp(kappa*(x(3,:)).^2), 'RelTol', 1E-10);
            testCase.verifyEqual(wZero.pdf(x), wZero.C*exp(kappa)*exp(-kappa* x(1,:).^2 + -kappa* x(2,:).^2  ), 'RelTol', 1E-10);
            B = BinghamDistribution(-kappa*[1, 1, 0]', eye(3,3));
            testCase.verifyEqual(wZero.pdf(x), B.pdf(x), 'RelTol', 1E-3); 
            
            %% test toBingham
            B2 = w.toBingham();
            B2.F = BinghamDistribution.computeF(B2.Z, 'mhg');
            testCase.verifyEqual(w.pdf(x), B2.pdf(x), 'RelTol', 1E-3);
            
            %% test stochastic sampling
            rng default
            n = 1000;
            samples = w.sample(n);
            testCase.verifySize(samples, [length(mu), n]);
            testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);
        end
        
        function testWatsonDistributionNegativeKappa(testCase)
            mu = [4, 1, 5]';
            mu = mu/norm(mu);
            kappa = -1.6;
            w = WatsonDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(w, 'WatsonDistribution');
            testCase.verifyEqual(w.mu, mu);
            testCase.verifyEqual(w.kappa, kappa);
            testCase.verifyEqual(w.dim, length(mu));

            %% test integral
            testCase.verifyEqual(w.integral(), 1, 'RelTol', 1E-5);
            
            %% test pdf
            kappa = 0.1;
            x = [0:10; 1:11; 5:15];
            x = bsxfun(@rdivide, x, sqrt(sum(x.^2,1))); %normalize
            wZero = WatsonDistribution([0,0,1]', kappa);
            testCase.verifyEqual(wZero.pdf(x), wZero.C*exp(kappa*(x(3,:)).^2), 'RelTol', 1E-10);
            testCase.verifyEqual(wZero.pdf(x), wZero.C*exp(kappa)*exp(-kappa* x(1,:).^2 + -kappa* x(2,:).^2  ), 'RelTol', 1E-10);
            B = BinghamDistribution(-kappa*[1, 1, 0]', eye(3,3));
            testCase.verifyEqual(wZero.pdf(x), B.pdf(x), 'RelTol', 1E-3); 
            
            %% test toBingham
            %B2 = w.toBingham();
            %B2.F = BinghamDistribution.computeF(B2.Z, 'mhg');
            %testCase.verifyEqual(w.pdf(x), B2.pdf(x), 'RelTol', 1E-3);
            
            %% test stochastic sampling
            %rng default
            %n = 1000;
            %samples = w.sample(n);
            %testCase.verifySize(samples, [length(mu), n]);
            %testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);
        end
    end
end