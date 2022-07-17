classdef BayesianComplexWatsonMixtureModel < AbstractComplexHypersphericalDistribution
    % BayesianComplexWatsonMixtureModel represents an extension to the complex
    % Watson mixture model by adding a prior distribution for the mode vectors
    % and a prior distribution for the mixture weights.
    %
    properties
        
        % set of complex bingham parameter matrices
        B
        
        % Column vector of concentration parameters
        concentrations
        
        % Number of mixture components
        K
        
        % Dirichlet parameter vector
        alpha
    end
    
    methods
        function bcWMM = BayesianComplexWatsonMixtureModel( ...
                B_, concentrations_, alpha_)
            
            assert(isequal(B_, permute(conj(B_), [2 1 3:ndims(B_)])), ...
                'B must be hermitian')
            
            assert(size(B_, 3) == size(alpha_, 1), '')
            assert(size(concentrations_, 1) == size(alpha_, 1), '')
            
            bcWMM.K = size(alpha_, 1);
            bcWMM.dim = size(B_, 1);
            
            bcWMM.B = B_;
            bcWMM.concentrations = concentrations_;
            bcWMM.alpha = alpha_;
        end
        
        function pdf(~, ~)
            error('not implemented yet');
        end
    end
    
    methods (Static)
        function [bcWMM, posterior] = fit(Z, parameters)
            
            posterior = BayesianComplexWatsonMixtureModel.estimatePosterior(Z, parameters);
            bcWMM = BayesianComplexWatsonMixtureModel(...
                posterior.B, posterior.kappa, posterior.alpha);
        end
        
        function [bcWMM, posterior] = fitDefault(Z, K)
            % :param Z: Complex observations with dimensions N times D.
            % :param K: Number of mixture components
            % assume, that the observation vector has fewer than 100 components. 
            assert(size(Z, 1) < 100) 
            
            D = size(Z, 1);
            parameters = BayesianComplexWatsonMixtureModel.parametersDefault(D, K);

            [bcWMM, posterior] = BayesianComplexWatsonMixtureModel.fit(Z, parameters);
            
        end

    function [parameters] = parametersDefault(D, K)
            % :param D:
            % :param K: Number of mixture components
            % assume, that the observation vector has fewer than 100 components.

            parameters = struct();
            parameters.initial = struct();
            parameters.initial.B = zeros(D, D, K);
            % Alternative: parameters.initial.B = rand(D, D, K);
            %              parameters.initial.B = conj(parameters.initial.B) .* permute(parameters.initial.B, [2 1 3]);
            parameters.initial.kappa = 20;
            parameters.initial.alpha = 1/K + linspace(-0.14/K, 0.14/K, K)';
            % Alternative: parameters.initial.alpha = ones(K, 1) / K;
            parameters.prior = struct();
            parameters.prior.B = zeros(D, D, K);
            parameters.prior.alpha = ones(K, 1) / K;
            parameters.prior.saliencies = 1;
            parameters.I = 40;
            parameters.uniformComponent = false;
        end

        function posterior = estimatePosterior(Z, parameters)
            % This is a quick and dirty implementation of a complex Watson mixture model
            % with a complex Bingham prior for the mode vectors and a Dirichlet prior
            % for the mixture weights.
            %
            % :param Z: Complex observations with dimensions N times D.
            % :param uniformComponent: If true, estimate the last component as
            %    uniform (means kappa(conzentration is zero). Example usecase: One component is noise.
            
            if ~ isfield(parameters, 'uniformComponent')
                parameters.uniformComponent = false;
            end
                        
            assert(isfield(parameters, 'initial'));
            assert(isfield(parameters.initial, 'B'));
            assert(isfield(parameters.initial, 'alpha'));
            assert(isfield(parameters.initial, 'kappa'));
            
            assert(isfield(parameters, 'prior'));
            assert(isfield(parameters.prior, 'B'));
            assert(isfield(parameters.prior, 'alpha'));
            
            assert(isequal(parameters.initial.B, ...
                permute(conj(parameters.initial.B), ...
                [2 1 3:ndims(parameters.initial.B)])), ...
                'initial B must be hermitian')
            
            assert(isfield(parameters, 'I'));
            
            [D, N] = size(Z);
            K_ = size(parameters.initial.alpha, 1);
            
            posterior.B = parameters.initial.B;
            posterior.alpha = parameters.initial.alpha;
            posterior.kappa = parameters.initial.kappa;
            posterior.gamma = zeros(N, K_);
            
            % Shiftdim to add another dimension in front of all, repmat to expand
            % this new dimension to D to yield a DxDxFxT array.
            dyadicProducts = repmat(shiftdim(Z, -1), [D, 1, 1]);
            
            % Transmute equals exchange the first and second dimension
            dyadicProducts = conj(dyadicProducts) .* permute(dyadicProducts, [2 1 3 ]);
            
            ZZ = reshape(bsxfun(@times, permute(Z, [1, 3, 2]), permute(conj(Z), [3 1 2])), D*D, N);
            
            if isfield (parameters.prior, 'saliencies')
                ln_saliencies = log(max(parameters.prior.saliencies, 1e-7));
            else
                ln_saliencies = log(1e-7);
            end
            
            if size(ln_saliencies, 1) == N && size(ln_saliencies, 2) == 1
                ln_saliencies = repmat(ln_saliencies, 1, K_);
            end
            
            for i = 1:parameters.I
                %% E-step
                % Contribution due to the observations
                posterior.gamma = ln_saliencies;
                
                posterior.gamma = posterior.gamma + bsxfun(@times, ...
                    posterior.kappa.', ...
                    BayesianComplexWatsonMixtureModel.quadraticExpectation( ...
                    dyadicProducts, posterior.B));
                
                % Contribution due to the complex Watson normalization constant
                posterior.gamma = bsxfun( ...
                    @plus, ...
                    posterior.gamma, ...
                    ComplexWatsonDistribution.logNorm(D, posterior.kappa).' ...
                    );
                
                % Expectation of logarithm of mixture weights
                posterior.gamma = bsxfun( ...
                    @plus, ...
                    posterior.gamma, ...
                    psi(posterior.alpha).' ...
                    );
                
                % Set maximum value to have highest precision for the most
                % relevant posterior probabilities.
                posterior.gamma = bsxfun( ...
                    @minus, ...
                    posterior.gamma, ...
                    max(posterior.gamma, [], 2) ...
                    );
                
                posterior.gamma = exp(posterior.gamma);
                
                % Since the posterior for the class labels is a categorial
                % distribution, each gamma has to sum up to one
                posterior.gamma = bsxfun( ...
                    @rdivide, ...
                    posterior.gamma, ...
                    sum(posterior.gamma, 2) ...
                    );
                
                assert(~anynan(posterior.gamma(:)));
                
                %% M-step
                numberOfObservationsPerClass = sum(posterior.gamma, 1);
                
                posterior.alpha = parameters.prior.alpha + ...
                    numberOfObservationsPerClass.';
                
                concentrationMax = 500;
                % covarianceMatrix = zeros(D, D, K_);
                % for k = 1:K_
                    % weights = permute(posterior.gamma(:, k), [3 1 2]);
                    % weights = repmat(weights / sum(weights), D, 1);
                    % covarianceMatrix(:, :, k) = (weights .* Z) *Z';
                % end               
                
                covarianceMatrix = reshape(...
                    bsxfun(@rdivide, ZZ * posterior.gamma, sum(posterior.gamma)), ...
                    [D, D, K_]);
                
                posterior.B = bsxfun(@times, ...
                    permute( posterior.kappa .* ...
                    numberOfObservationsPerClass.', [2 3 1]), ...
                    covarianceMatrix) + parameters.prior.B;
                
                % force B to be hermitian (Generalization of 1/2*(B+B')).
                posterior.B = 1/2*(posterior.B + ...
                    permute(conj(posterior.B), [2 1 3:ndims(posterior.B)]));
                
                quadraticExpectation = zeros(K_, 1);
                for k = 1:K_
                    quadraticExpectation(k) = ...
                        BayesianComplexWatsonMixtureModel.quadraticExpectation( ...
                        covarianceMatrix(:, :, k), posterior.B(:, :, k));
                end
                
                posterior.kappa = ...
                    hypergeometricRatioInverse(quadraticExpectation, D, concentrationMax);
                
                if parameters.uniformComponent
                    posterior.kappa(end) = 0;
                end
            end
        end
        function E = quadraticExpectation(dyadicProducts, B)
            % This function calculates the quadratic estimate E{X^H B X} and
            % E{tr(Phi_XX B)}.
            %
            % :param dyadicProducts: X * X^H or Phi_XX
            % :param B: Complex Bingham parameter matrix of the underlying distribution.
            
            assert(isequal(B, permute(conj(B), [2 1 3:ndims(B)])), ...
                'initial B must be hermitian')
            
            N = size(dyadicProducts, 3);
            D = size(B, 1);
            K_ = size(B, 3);
            E = zeros(N, K_);
            
            dyadicProductsReshape = reshape(dyadicProducts, [], N);
            
            evalString = sprintf('@(x) cBinghamGradNormDividedByNorm%d(transpose(x))', D);
            firstOrderMoments = str2func(evalString);
            
            for k = 1:K_
                [U, Lambda] = eig(B(:, :, k), 'Vector');
                
                assert(isreal(Lambda), ...
                    'all Eigenvalues of B have to be real, since B is hermitian')
    
                %% Calculate first order moments roughly
                % Covariance matrix without rotation
                if sum(Lambda > 1)
                    Lambda = Lambda.'+(1:D)*1e-2;
                    C = diag(firstOrderMoments(Lambda-max(Lambda)));
                else
                    C = 1/D*eye(D);
                end
                
                % Covariance matrix with rotation
                covarianceMatrix = U * C * U';
                
                % E(n, k) = real(trace(dyadicProducts(:, :, n) * covarianceMatrix));
                % http://stackoverflow.com/questions/8031628/octave-matlab-efficient-calc-of-frobenius-inner-product
%               E(:, k) = real(sum(sum(bsxfun(@times, dyadicProducts, covarianceMatrix.'))));
                E(:, k) = real(dyadicProductsReshape.' * conj(covarianceMatrix(:)));

            end
        end
        
    end
end
