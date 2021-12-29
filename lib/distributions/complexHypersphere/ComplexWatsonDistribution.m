classdef ComplexWatsonDistribution < AbstractComplexHypersphericalDistribution
    % ComplexWatsonDistribution represents a complex Watson distribution
    %
    % Mardia, K. V. & Dryden, I. L. 
    % The Complex Watson Distribution and Shape Analysis 
    % Journal of the Royal Statistical Society: Series B 
    % (Statistical Methodology), 
    % Blackwell Publishers Ltd., 1999, 61, 913-926
    
    properties
        kappa   % concentration (scalar)
        mu      % mode vector
        C       % normalization constant
        logC    % logarithmic normalization constant
    end
    
    methods
        function W = ComplexWatsonDistribution(mu_, kappa_)
            %% Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            epsilon = 1E-6;
            assert(size(mu_,2) == 1, 'mu must be a row vector');
            assert(abs(norm(mu_) - 1)<epsilon, 'mu must be normalized');
            %assert(kappa_>0, 'kappa must be postive');
            
            W.mu = mu_;
            W.kappa = kappa_;
            
            W.dim = size(mu_,1);
            
            W.logC = ComplexWatsonDistribution.logNorm(W.dim, W.kappa);
            W.C = exp(W.logC);
        end

        function m = mean(this)
            m = this.mu;
        end
        
        function p = pdf(this, za)
            % Evaluates pdf at each column of za
            % Parameters:
            %   za (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of za
            
            p = exp(this.logC + this.kappa * abs(this.mu' * za).^2);
        end
        
        function Z = sample(this, N)
            % Use special case of complex Bingham distribution to sample from
            % a complex Watson distribution according to
            % Mardia, K. V. & Jupp, P. E. Directional statistics Wiley, 2009, page 336
            %
            % :param N: Number of requested samples
            % :type N: Integer scalar
            % :return Z: Sampled observation vectors,
            %   complex matrix with size #dimensions times #samples
            
            % Calculate the Bingham parameter matrix
            B = - this.kappa * (eye(length(this.mu)) - this.mu*this.mu');
            
            if ~ismatrix(B) || size(B, 1) ~= size(B, 2)
                error('sampleComplexBingham:B', 'Matrix B is not square.');
            end
            
            P = size(B, 1); % Z will have dimension P*N
            if P < 2
                error('sampleComplexBingham:B', 'Matrix B is too small.');
            end
            
            % Calculate eigenvalues labda_i of -B (negative)
            % force B to be hermitian
            B = 1/2*(B+B');
            [V, D] = eig(-B);
            
            Lambda = diag(D);
            
            % Sort eigenvalues in descending order and arrange corresponding
            % eigenvectors
            [Lambda, IX] = sort(Lambda, 'descend');
            V = V(:, IX);
            
            % Perform the eigenvalue shift. This does not affect the eigenvectors
            % since A and A+kappa*I represent the same distribution.
            Lambda = Lambda - Lambda(end);
            
            temp1 = -(1./Lambda(1:end-1));
            temp2 = (1-exp(-Lambda(1:end-1)));
            
            Z = zeros(P, N);
            for n = 1:N % parfor
                S = zeros(P, 1);
                while true
                    % Sample uniform random variable U_j ~ U(0, 1), j = 1, ... k-1
                    U = rand(P-1, 1);
                    
                    % Transform so that result are independent Texp(Lambda_i) random
                    % variables
                    for i = 1:P-1
                        % Exponential distribution is nearly uniform for small Lambda
                        % This speeds up the sampling a lot
                        % if abs(diff(exppdf([0, 1], 1/Lambda(i)))) < 1e-3
                        if Lambda(i) < 0.03 % equivalent but faster condition
                            S(i) = U(i);
                        else
                            S(1:end-1) = temp1 ...
                                .* log(1 - U.*temp2);
                            % S(1:end-1) = -(1./Lambda(1:end-1)) ...
                            %     .* log(1 - U.*(1-exp(-Lambda(1:end-1))));
                        end
                    end
                    
                    % Reject S if the sum is larger or equal to 1.
                    if sum(S) < 1
                        break;
                    end
                end
                
                % The last entry was still missing.
                S(end) = 1 - sum(S);
                
                % Generate idependent angles
                theta = 2*pi*rand(P, 1);
                
                W = zeros(P, 1);
                for i = 1:P
                    W(i) = sqrt(S(i))*exp(1i*theta(i));
                end
                
                % Matrix V consists of right eigenvectors of B
                Z(:, n) = V*W;
            end
        end
    end
    
    methods (Static)
        function cW = fit(Z, weights)
            % according to "The complex Watson distribution and shape
            % analysis" by K.V. Mardia and I.L. Dryden
            % Z: D x N
            % weight: 1 x N
            if nargin<2
                [mu_hat, kappa_hat] = ...
                    ComplexWatsonDistribution.estimateParameters(Z);
                cW = ComplexWatsonDistribution(mu_hat, kappa_hat);
            else
                [mu_hat, kappa_hat] = ...
                    ComplexWatsonDistribution.estimateParameters(Z, weights);
                cW = ComplexWatsonDistribution(mu_hat, kappa_hat);
            end
        end
        
        function [mu_hat, kappa_hat] = estimateParameters(Z, weights)
            % according to "The complex Watson distribution and shape
            % analysis" by K.V. Mardia and I.L. Dryden
            % Z: D x N
            % weight: 1 x N
            
            N = size(Z, 2);           
            
            if nargin<2
                S = Z*Z';
            else
                assert(size(weights,1)==1, 'weights is no row vector');
                assert(size(weights,2)==N, ...
                    'dimensions of Z and weights mismatch');
                S = bsxfun(@times,weights ,Z) *Z'* N/sum(weights);
            end
            
            % force S to be hermitian
            S = 1/2*(S+S');
            
            [V, lambda] = eig(S, 'Vector');
            
            assert(all(lambda > 0));
            % Sort eigenvalues in descending order and arrange corresponding
            % eigenvectors
            [lambda, IX] = sort(lambda, 'descend');
            V = V(:, IX);
            
            mu_hat = V(:, 1);
            
            D = size(S, 1);
            
            % approximation: only valid under high concentration
            kappa_hat = N*(D-1)/(N - lambda(1));
            
            approximation_threshold = 200;
            
            if kappa_hat < approximation_threshold
                concentrationMax = 1000;
                normed_lamda = lambda(1)/N;
                kappa_hat = hypergeometricRatioInverse(normed_lamda, D, concentrationMax);
            end
        end
        
        function log_c = logNorm(D, kappa)
            % Calculates the normalization constant for the compex Watson distribution.
            % The only available variant is the reformulation of Mardia1999Watson in
            % logarithmic domain.
            %
            % The logNorm function is more robust for very high concentrations.
            %
            % :param D: Dimension of feature space.
            % :type D: Scalar.
            % :param kappa: Concentration parameter.
            % :type kappa: Any size, numerical.
            % :return log_c: Normalization factor, size is equal to sie of kappa.
            
            inputDimensions = size(kappa);
            kappa = kappa(:).';
            log_c = zeros(size(kappa));
            
            % Mardia1999Watson Page 917 (no equation number)
            log_c_high = log(2) + D*log(pi) + (1-D)*log(kappa) + kappa;
            
            % Mardia1999Watson Equation (3)
            log_c_medium = log_c_high + log(1-sum(bsxfun(@rdivide, ...
                bsxfun(@times, ...
                bsxfun(@power, kappa, (0:D-2).'), ...
                exp(-kappa)), ...
                factorial(0:D-2).'), 1));
            
            % Mardia1999Watson Equation (4), Taylor series up to power of 10
            log_c_low = log(2) + D*log(pi) - log(factorial(D-1)) ...
                + log(1 + sum(cumprod( ...
                bsxfun(@rdivide, repmat(kappa, 10, 1), (D:D+10-1).'), 1)));
            
            % A good boundary between low and medium is kappa = 1/D. This can be
            % motivated by plotting each term for very small values and for
            % values from [0.1, 1].
            % A good boundary between medium and high is when
            % kappa * exp(kappa) < epsilon. Choosing kappa = 100 is sufficient.
            log_c(kappa < 1/D) = log_c_low(kappa < 1/D);
            log_c((1/D <= kappa) & (kappa < 100)) ...
                = log_c_medium((1/D <= kappa) & (kappa < 100));
            log_c(100 <= kappa) = log_c_high(100 <= kappa);
            
            log_c = -reshape(log_c, inputDimensions);
        end
    end
end