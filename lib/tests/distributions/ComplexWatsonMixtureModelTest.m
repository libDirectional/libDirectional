classdef ComplexWatsonMixtureModelTest < matlab.unittest.TestCase
    properties
		N
		D
        K
        
        cWMM_highConcentrations
		Z_highConcentrations
		label_highConcentrations
        
        cWMM_lowConcentrations
		Z_lowConcentrations
		label_lowConcentrations
    end
 
    methods(TestClassSetup)
        function createSampling(tc)
			%% Define mean shape (different ideas)
			mu1 = [1+.5i -1+.5i -1 -1-.5i].';
			mu1 = mu1/norm(mu1);

			tc.D = 4;
            tc.K = 2;
			mu2 = exp((0:tc.D-1)'*2/tc.D*pi*1i);
			mu2 = mu2/norm(mu2);

			modeSet = [mu1, mu2];

			%% Define concentration parameter
			lowConcentrations = [20; 10];
			highConcentrations = [200; 100];
			tc.N = 5000;
            
			mixtureWeights = [0.7; 0.3];
            
            tc.cWMM_highConcentrations = ...
                ComplexWatsonMixtureModel( ...
                modeSet, highConcentrations, mixtureWeights);
            
            tc.cWMM_lowConcentrations = ...
                ComplexWatsonMixtureModel( ...
                modeSet, lowConcentrations, mixtureWeights);
            
			%% Sampling
			[tc.Z_highConcentrations, tc.label_highConcentrations] = ...
                tc.cWMM_highConcentrations.sample(tc.N);
            [tc.Z_lowConcentrations, tc.label_lowConcentrations] = ...
                tc.cWMM_lowConcentrations.sample(tc.N);
            
        end
    end
    methods (Test)
        function testIntegralLowConcentrations(tc)
            %% test integral
            tc.verifyEqual( ...
                tc.cWMM_lowConcentrations.integral(), 1, 'RelTol', 0.1);
            
        end
        function checkSizeOfSampleOutput(tc)
            %% verification
            verifySize(tc, tc.Z_highConcentrations, [tc.D, tc.N]);
        end
        function testFittingHighConcentrations(tc)
            %%
            
            cWMM_fit = ...
                ComplexWatsonMixtureModel.fit(tc.K, tc.Z_highConcentrations);
            
            %% verification
            
            [estimatedConcentrations, idx] = ...
                sort(cWMM_fit.concentrations, 'descend');
            estimatedMixtureWeights = cWMM_fit.mixtureWeights(idx);
            estimatedModeSet = cWMM_fit.modeSet(:, idx);
            
            tc.verifyEqual(estimatedConcentrations, ...
                tc.cWMM_highConcentrations.concentrations, 'RelTol', 0.05);
            tc.verifyEqual(estimatedMixtureWeights, ...
                tc.cWMM_highConcentrations.mixtureWeights, 'AbsTol', 0.05);
            
			normalizedEstimatedModeSet = ...
                ComplexWatsonMixtureModel.normalizedPhase(estimatedModeSet);
			normalizedTrueModeSet = ...
                ComplexWatsonMixtureModel.normalizedPhase( ...
                tc.cWMM_highConcentrations.modeSet);
			
			tc.verifyEqual(normalizedEstimatedModeSet, ...
                normalizedTrueModeSet, 'AbsTol', 1e-1);
            
        end
        function testFittingLowConcentrations(tc)
            %%
            
            cWMM_fit = ...
                ComplexWatsonMixtureModel.fit(tc.K, tc.Z_lowConcentrations);
            
            %% verification
            
            [estimatedConcentrations, idx] = ...
                sort(cWMM_fit.concentrations, 'descend');
            estimatedMixtureWeights = cWMM_fit.mixtureWeights(idx);
            estimatedModeSet = cWMM_fit.modeSet(:, idx);
            
            tc.verifyEqual(estimatedConcentrations, ...
                tc.cWMM_lowConcentrations.concentrations, 'RelTol', 0.05);
            tc.verifyEqual(estimatedMixtureWeights, ...
                tc.cWMM_lowConcentrations.mixtureWeights, 'AbsTol', 0.05);
            
			normalizedEstimatedModeSet = ...
                ComplexWatsonMixtureModel.normalizedPhase(estimatedModeSet);
			normalizedTrueModeSet = ...
                ComplexWatsonMixtureModel.normalizedPhase( ...
                tc.cWMM_lowConcentrations.modeSet);
			
			tc.verifyEqual(normalizedEstimatedModeSet, ...
                normalizedTrueModeSet, 'AbsTol', 1e-1);
            
        end
    end
end