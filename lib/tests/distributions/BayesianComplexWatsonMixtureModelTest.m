classdef BayesianComplexWatsonMixtureModelTest < matlab.unittest.TestCase
    properties
        W
        piVector
        kappaVector
        label
        Z
        N
        D
    end
    
    methods(TestClassSetup)
        function createSampling(tc)
            %% Define mean shape (different ideas)
            mu1 = [1+.5i -1+.5i -1 -1-.5i 1-.5i 1].'; %'
            mu1 = mu1/norm(mu1);
            
            tc.D = 6;
            mu2 = exp((0:tc.D-1).'*2/tc.D*pi*1i);
            mu2 = mu2/norm(mu2);
            
            tc.W = [mu1, mu2];
            
            %% Define concentration parameter
            tc.kappaVector = [50; 50];
            tc.N = 1000;
            
            tc.piVector = [0.3; 0.7];
            
            % Ensure sum of probabilites is correct
            tc.piVector = tc.piVector/sum(tc.piVector);
            
            %% Sampling
            cWMM = ComplexWatsonMixtureModel( ...
                tc.W, tc.kappaVector, tc.piVector);
            
            [tc.Z, tc.label] = cWMM.sample(tc.N);
        end
    end
    
    methods(Test)
        function em_algorithm_test_with_ground_true(tc)
            %% Initialized with ground true
            parameters.initial.B = cat(3, ...
                eye(tc.D) + tc.W(:, 1) * tc.W(:, 1)', ...
                eye(tc.D) + tc.W(:, 2) * tc.W(:, 2)' ...
                );
            
            % force B to be hermitian (Generalization of 1/2*(B+B')).
            parameters.initial.B = 1/2*(parameters.initial.B + ...
                permute(conj(parameters.initial.B ), ...
                [2 1 3:ndims(parameters.initial.B )]));
            
            parameters.initial.kappa = tc.kappaVector;
            parameters.initial.alpha = [300 700].';
            
            %% Set an uninformative prior
            parameters.prior.B = zeros(6, 6, 2);
            parameters.prior.alpha = [300 700].';
            
            %% Set further EM parameters
            parameters.I = 100;
            
            %% Execute EM_Algorithm, test if result is similar to ground true
            posterior = BayesianComplexWatsonMixtureModel.estimatePosterior( ...
                tc.Z, parameters);
            
            %% Check estimate for mixture weights
            piVectorEstimated = posterior.alpha ./ sum(posterior.alpha);
            tc.verifyEqual(piVectorEstimated, tc.piVector, 'AbsTol', 0.1);
            
            %% Check concentration parameters
            tc.verifyEqual(posterior.kappa, tc.kappaVector, 'AbsTol', 5);
            
            WSetEstimated = zeros(6, 2);
            for k = 1:2
                [eigenVectors, eigenValues] = eig(posterior.B(:, :, k));
                [~, index] = max(diag(eigenValues));
                WSetEstimated(:, k) = eigenVectors(:, index);
                tc.verifyEqual(abs(WSetEstimated(:, k)' * tc.W(:, k))^2, 1, 'AbsTol', 0.1);
            end
        end
        
        function em_algorithm_test_convergence(tc)
            
            B_init = eye(tc.D) + (tc.W(:, 1) * tc.W(:, 1)' + ...
                + tc.W(:, 2) * tc.W(:, 2)')/2;
            
            parameters.initial.B = repmat(B_init, 1, 1, 2);
            
            % force B to be hermitian (Generalization of 1/2*(B+B')).
            parameters.initial.B = 1/2*(parameters.initial.B + ...
                permute(conj(parameters.initial.B ), ...
                [2 1 3:ndims(parameters.initial.B )]));
            
            parameters.initial.kappa = [20; 20];
            parameters.initial.alpha = [490 510].';
            
            %% Set an uninformative prior
            parameters.prior.B = zeros(6, 6, 2);
            parameters.prior.alpha = [400 600].';
            
            %% Set further EM parameters
            parameters.I = 200;
            
            %% Execute EM_Algorithm, test if result is similar to ground true
            posterior = BayesianComplexWatsonMixtureModel.estimatePosterior( ...
                tc.Z, parameters);
            
            %% Check estimate for mixture weights
            piVectorEstimated = posterior.alpha ./ sum(posterior.alpha);
            tc.verifyEqual(piVectorEstimated, tc.piVector, 'AbsTol', 0.1);
            
            %% Check concentration parameters
            tc.verifyEqual(posterior.kappa, tc.kappaVector, 'AbsTol', 5);
            
            WSetEstimated = zeros(6, 2);
            for k = 1:2
                [eigenVectors, eigenValues] = eig(posterior.B(:, :, k));
                [~, index] = max(diag(eigenValues));
                WSetEstimated(:, k) = eigenVectors(:, index);
                tc.verifyEqual(abs(WSetEstimated(:, k)' * tc.W(:, k))^2, 1, 'AbsTol', 0.1);
            end
        end
        
        function testInstantiation(testCase)
            B = 2;
            concentrations = 3;
            alpha = 4;
            bcWMM = BayesianComplexWatsonMixtureModel(B, concentrations, alpha);
            testCase.verifyClass(bcWMM, 'BayesianComplexWatsonMixtureModel');
        end
    end
end
