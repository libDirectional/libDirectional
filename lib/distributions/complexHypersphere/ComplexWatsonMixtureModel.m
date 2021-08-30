classdef ComplexWatsonMixtureModel < AbstractComplexHypersphericalDistribution
    % ComplexWatsonMixtureModel represents a complex Watson Mixture model
    
    properties
        
        % D times K matrix of complex column mode vectors for each class
        % dependent Watson distribution
        modeSet
        
        % Column vector of concentration parameters
        concentrations
        
        % Number of mixture components
        K
        
        % Column vector of weights for each mixture component
        mixtureWeights
    end
    
    methods
        function cWMM = ComplexWatsonMixtureModel( ...
                modeSet_, concentrations_, mixtureWeights_)
            
            epsilon = 1E-6;
            assert(size(modeSet_,2) == size(mixtureWeights_, 1), ...
                'dimensions of modeSet_ and mixtureWeights_ do not fit');
            assert(size(modeSet_,2) == size(concentrations_, 1), ...
                'dimensions of modeSet_ and concentrations_ do not fit');
            assert(all(sum(abs(modeSet_).^2, 1) - 1<epsilon), ...
                'modeSet_ must be normalized');
            
            cWMM.modeSet = modeSet_;
            cWMM.concentrations = concentrations_;
            cWMM.mixtureWeights = mixtureWeights_;
            
            cWMM.dim = size(modeSet_, 1);
            cWMM.K = size(mixtureWeights_, 1);
            
        end
        
        function [marginalPdf, jointPdf, conditionalPdf] = pdf(this, za)
            %PDF Evaluates pdf at each column of za
            % Parameters:
            %   za (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   *Pdf (1 x n row vector)
            %       values of pdf at each column of za
            
            logC = ...
                ComplexWatsonDistribution.logNorm(this.dim, this.concentrations);
            
            exponents = ...
                bsxfun(@times, this.concentrations, abs(this.modeSet' * za).^2);
            
            conditionalPdf = exp(bsxfun(@plus, logC, exponents));
            
            jointPdf = bsxfun(@times, this.mixtureWeights, ...
                conditionalPdf);
            
            marginalPdf = sum(jointPdf, 1);
            
        end
        
        function [Z, C] = sample(this, N)
            % Sample from a complex Watson mixture model
            %
            % :param N: Number of requested samples
            % :return Z: Sampled observation vectors,
            %   complex matrix with size #dimensions times #samples
            % :return C: Row vector of N class labels
            
            % Sample all category labels
            C = discretesample(this.mixtureWeights, N);
            % C = sort(C);
            
            % Sample Watson distribution for each class
            Z = zeros(size(this.modeSet, 1), N);
            for k = 1:this.K
                index = (C == k);
                N_k = sum(index);
                
                % Only set the random vectors which are labeled with class k
                cW = ComplexWatsonDistribution( ...
                    this.modeSet(:, k), this.concentrations(k));
                
                Z(:, index) = cW.sample(N_k);
            end
        end
    end
    
    methods (Static)
        function [cWMM, gamma, gamma_kmean] = fit(K, Z, maxIterations)
            if nargin < 3
                maxIterations = 100;
            end
            
            % Init
            % normalize Z
            Z = bsxfun(@rdivide, Z, sqrt(sum(abs(Z).^2, 1)));
            
            N = size(Z, 2);
            D = size(Z, 1);
            
            % kmeans
            Z_kmeans = ComplexWatsonMixtureModel.normalizedPhase(Z);
            label_kmeans = (kmeans(vertcat(real(Z_kmeans), imag(Z_kmeans))', ...
				                   K, 'start', 'sample'))';
            gamma_kmean = zeros(K, N);
            for ii = 1 : K
                gamma_kmean(ii, :) = label_kmeans == ii;
            end
            N_k = sum(gamma_kmean, 2);
            InitMixtureWeights = N_k / sum(N_k);
            
            InitModeSet = zeros(D, K);
            InitConcentrations = zeros(K, 1);

%           I = randperm(N, K);
%           InitModeSet = Z(:, I, :);
%           InitMixtureWeights = zeros(K, 1) + 1 / K;
%           InitConcentrations = zeros(K, 1) + 20;

            for k = 1:K
                [InitModeSet(:, k), InitConcentrations(k)] = ...
                    ComplexWatsonDistribution.estimateParameters( ...
                    Z, gamma_kmean(k, :));
            end
            
            cWMM = ComplexWatsonMixtureModel( ...
                InitModeSet, InitConcentrations, InitMixtureWeights);
            
            for iterations = 1 : maxIterations
                %% E-Step
                [marginalPdf, jointPdf] = cWMM.pdf(Z);
                gamma = bsxfun(@rdivide, jointPdf, marginalPdf);
                    
                %% M-Step
                % Calculate the expected Number of samples from source k
                N_k = sum(gamma, 2);
                
                % Calculate the mixtureWeights
                cWMM.mixtureWeights = N_k / sum(N_k);
                
                for k = 1:cWMM.K
                    [cWMM.modeSet(:, k), cWMM.concentrations(k)] = ...
                        ComplexWatsonDistribution.estimateParameters( ...
                        Z, gamma(k, :));
                end
                
            end

            if nargout > 1
                gamma = gamma';
            end
            if nargout > 2
                gamma_kmean = gamma_kmean';
            end
        end
        
        function Z_out = normalizedPhase(Z_in)
            phase = angle( Z_in(1, :) );
            Z_out = bsxfun(@times, Z_in, exp(-1i*phase));            
        end
    end
end
