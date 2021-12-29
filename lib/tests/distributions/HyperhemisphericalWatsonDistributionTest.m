classdef HyperhemisphericalWatsonDistributionTest< matlab.unittest.TestCase
   
    methods (Test)
        function testHyperhemisphericalWatsonDistribution(testCase)
            mu = [1, 2, 3]';
            mu = mu/norm(mu);
            kappa = 2;
            w = HyperhemisphericalWatsonDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(w, 'HyperhemisphericalWatsonDistribution');
            testCase.verifyEqual(w.mode(), mu);
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
            wZeroHemi = HyperhemisphericalWatsonDistribution([0,0,1]', kappa);
            testCase.verifyEqual(wZeroHemi.pdf(x), 2*wZero.pdf(x), 'RelTol', 1E-10);
            B = BinghamDistribution(-kappa*[1, 1, 0]', eye(3,3));
            testCase.verifyEqual(wZeroHemi.pdf(x), 2*B.pdf(x), 'RelTol', 1E-3); 
            
            %% test stochastic sampling
            rng default
            n = 1000;
            samples = w.sample(n);
            testCase.verifySize(samples, [length(mu), n]);
            testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);
        end
        
        function testHyperhemisphericalWatsonDistributionNegativeKappa(testCase)
            mu = [4, 1, 5]';
            mu = mu/norm(mu);
            kappa = -1.6;
            w = HyperhemisphericalWatsonDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(w, 'HyperhemisphericalWatsonDistribution');
            testCase.verifyEqual(w.mode(), mu);
            testCase.verifyEqual(w.dim, length(mu));

            %% test integral
            testCase.verifyEqual(w.integral(), 1, 'RelTol', 1E-5);
            
            %% test pdf
            kappa = 0.1;
            x = [0:10; 1:11; 5:15];
            x = bsxfun(@rdivide, x, sqrt(sum(x.^2,1))); %normalize
            wZero = WatsonDistribution([0,0,1]', kappa);
            wZeroHemi = HyperhemisphericalWatsonDistribution([0,0,1]', kappa);
            testCase.verifyEqual(wZeroHemi.pdf(x), 2*wZero.pdf(x), 'RelTol', 1E-10);
            B = BinghamDistribution(-kappa*[1, 1, 0]', eye(3,3));
            testCase.verifyEqual(wZeroHemi.pdf(x), 2*B.pdf(x), 'RelTol', 1E-3); 
        end
        function testShifting(testCase)
            dist = HyperhemisphericalWatsonDistribution([0; 0; 0; 1], 0.5);
            distShifted = dist.shift([1;2;3;4]/norm([1;2;3;4]));
            testCase.verifyEqual(distShifted.mu, [1;2;3;4]/norm([1;2;3;4]), 'RelTol', 1E-10);
        end
    end
end