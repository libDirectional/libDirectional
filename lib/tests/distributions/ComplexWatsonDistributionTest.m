classdef ComplexWatsonDistributionTest < matlab.unittest.TestCase
        properties
		N
		D
        
        cW_highKappa
		Z_highKappa
        
        cW_lowKappa
		Z_lowKappa
    end
 
    methods(TestClassSetup)
        function createSampling(tc)
			%% Define mean shape (different ideas)
			mu = [1+.5i -1+.5i -1 -1-.5i].';
			mu = mu/norm(mu);

			tc.D = 4;
            
			%% Define concentration parameter
			% ToDo: test with high kappas
			lowKappa = 20;
			highKappa = 200;
			tc.N = 5000;
            
            tc.cW_highKappa = ComplexWatsonDistribution(mu, highKappa);
            tc.cW_lowKappa = ComplexWatsonDistribution(mu, lowKappa);
            
            
			%% Sampling
			tc.Z_highKappa = ...
                tc.cW_highKappa.sample(tc.N);
            tc.Z_lowKappa = ...
                tc.cW_lowKappa.sample(tc.N);
            
        end
    end
    methods (Test)
        function sanityCheck(tc)
            %% sanity check
            mu = [1+.5i -1+.5i -1-.5i 1-.5i].';
            mu = mu/norm(mu);
            kappa = 20;
            w = ComplexWatsonDistribution(mu, kappa);
            
            %% verification
            tc.verifyClass(w, 'ComplexWatsonDistribution');
            tc.verifyEqual(w.mu, mu);
            tc.verifyEqual(w.kappa, kappa);
            tc.verifyEqual(w.dim, length(mu));
        end
        function testIntegral(tc)
            %% test integral
            tc.verifyEqual(tc.cW_lowKappa.integral(), 1, 'RelTol', 0.1);
            
        end
        function checkSizeOfSampleOutput(tc)
            %% verification
            verifySize(tc, tc.Z_highKappa, [tc.D, tc.N]);
        end
        
        function testFittingHighConcentration(tc)
            
            cW_fit = ComplexWatsonDistribution.fit(tc.Z_highKappa);
            
            %% verification
            
            tc.verifyEqual(cW_fit.kappa, ...
                tc.cW_highKappa.kappa, 'RelTol', 0.05);
            
			normalizedEstimatedMu = ...
                ComplexWatsonDistributionTest.normalizedPhase(cW_fit.mu);
			normalizedTrueMu = ...
                ComplexWatsonDistributionTest.normalizedPhase( ...
                tc.cW_highKappa.mu);
			
			tc.verifyEqual(normalizedEstimatedMu, ...
                normalizedTrueMu, 'AbsTol', 1e-1);
            
        end
        function testFittingLowConcentration(tc)
            %%
            % fails if accurate estimation of kappa is not implemented yet
            
            cW_fit = ComplexWatsonDistribution.fit(tc.Z_lowKappa);
            
            %% verification
            
            tc.verifyEqual(cW_fit.kappa, ...
                tc.cW_lowKappa.kappa, 'RelTol', 0.05);
            
			normalizedEstimatedMu = ...
                ComplexWatsonDistributionTest.normalizedPhase(cW_fit.mu);
			normalizedTrueMu = ...
                ComplexWatsonDistributionTest.normalizedPhase( ...
                tc.cW_lowKappa.mu);
			
			tc.verifyEqual(normalizedEstimatedMu, ...
                normalizedTrueMu, 'AbsTol', 1e-1);
            
        end
        function testFittingWithWeights(tc)
            weights = 0.5 + rand(1, tc.N);
            
            cW_fit = ComplexWatsonDistribution.fit(tc.Z_highKappa, weights);
            
            %% verification
            
            tc.verifyEqual(cW_fit.kappa, ...
                tc.cW_highKappa.kappa, 'RelTol', 0.05);
            
			normalizedEstimatedMu = ...
                ComplexWatsonDistributionTest.normalizedPhase(cW_fit.mu);
			normalizedTrueMu = ...
                ComplexWatsonDistributionTest.normalizedPhase( ...
                tc.cW_highKappa.mu);
			
			tc.verifyEqual(normalizedEstimatedMu, ...
                normalizedTrueMu, 'AbsTol', 1e-1);
            
        end
    end
    methods (Static)
        function Z_out = normalizedPhase(Z_in)
            
            phase = angle( Z_in(1, :) );
            Z_out = bsxfun(@times, Z_in, exp(-1i*phase));
            
        end
    end
end

