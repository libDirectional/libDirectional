classdef ComplexBinghamDistribution < AbstractComplexHypersphericalDistribution
    % ComplexBinghamDistribution represents a complex Bingham distribution
    %
    % Kent, J. T. 
    % The Complex Bingham Distribution and Shape Analysis 
    % Journal of the Royal Statistical Society. Series B (Methodological), 
    % JSTOR, 1994, 285-299
    
    properties
        B       % parameter matrix
        C       % normalization constant
        logC    % logarithmic normalization constant
    end
    
    methods
        function cB = ComplexBinghamDistribution(B_)
            %% Constructor
            %
            % Parameters:
            %   B_ (d x d): parameter matrix
            
            assert(all(all(B_ == B_')), 'B must be hermitian');
            
            cB.B = B_;
            cB.dim = size(B_, 1);
            
            cB.logC = ComplexBinghamDistribution.logNorm(cB.B);
            cB.C = exp(cB.logC);
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
            
            p = exp(this.logC + abs(dot(za, this.B * za)));
        end
        
        function Z = sample(this, N)
            % Samples from a complex Bingham distribution according to the paper
            % Kent, J. T.; Constable, P. D. L. & Er, F.,
            % Simulation for the complex Bingham distribution Statistics and Computing,
            % Springer, 2004, 14, 53-57
            %
            % :param N: Number of requested samples
            % :type N: Integer scalar
            % :return Z: Sampled observation vectors,
            %   complex matrix with size #dimensions times #samples
            
            if ~ismatrix(this.B) || size(this.B, 1) ~= size(this.B, 2)
                error('sampleComplexBingham:B', 'Matrix B is not square.');
            end
            
            P = size(this.B, 1); % Z will have dimension P*N
            if P < 2
                error('sampleComplexBingham:B', 'Matrix B is too small.');
            end
            
            % Calculate eigenvalues lambda_i of -B (negative)
            [V, D] = eig(-this.B);
            
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
        function cB = fit(Z)
            N = size(Z, 2);
            S = Z*Z'/N;
            % S_test = random.covariance(Z, [], 2, 1);
            B_ = ComplexBinghamDistribution.estimateParameterMatrix(S);
            % B_test = random.cBingham.maximumLikelihoodFit(S)
            cB = ComplexBinghamDistribution(B_);
        end
        
        function [B, eigenvectorsB, eigenvaluesB] = estimateParameterMatrix(S)
            %maximumLikelihoodFit Calculate ML-estimate of complex Bingham parameter matrix
            %   [Kent1994Complex] describes an ML-estimation for the complex Bingham
            %   distribution.
            %
            %   S ... Sensors * Sensors * Speakers Sum of squares and products matrix.
            %
            %   This function is only able to determine the ...
            
            options = optimoptions(@lsqnonlin, ...
                'TolFun', 1e-20, ...
                'MaxFunEvals', 1e5, ...
                'MaxIter', 1e4, ...
                'Display', 'off' ...
                );
            
            D = size(S, 1);
            K = size(S, 3);
            
            eigenvectorsB = zeros(D, D, K);
            eigenvaluesB = zeros(D, K);
            B = zeros(D, D, K);
            
            % TODO: Better initialization possible with eq. (3.3)
            X0 = 100*(-D+1:-1)/D;
            for k = 1:K
                [eigenvectorsS, temp] = eig(S(:, :, k));
                eigenvaluesS = diag(temp);
                
                % TODO: Remove this unorthodox correction:
                % eigenvaluesS = eigenvaluesS + 1e-4*(0:D-1).';
                
                % TODO: Allow arbitrary dimensions
                
                functionHandle = ...
                    str2func(['cBinghamGradLogNorm', num2str(D)]);
                
                f = @(X) functionHandle([X; 0]).' - eigenvaluesS;
                
                myX = lsqnonlin(f, X0.', ...
                    -1000*ones(1, D-1).', ...
                    -1e-2*ones(1, D-1).', ...
                    options);
                eigenvaluesB(:, k) = [myX; 0];
                eigenvectorsB(:, :, k) = eigenvectorsS;
                
                % TODO: Check this: Only correct if order of eigenvalues and order of
                % eigenvectors is corresponding.
                B(:, :, k) ...
                    = eigenvectorsB(:, :, k) ...
                    * diag(eigenvaluesB(:, k)) ...
                    * eigenvectorsB(:, :, k)';
                B(:, :, k) = 1/2*(B(:, :, k) + B(:, :, k)');
            end
        end
        
        function C = logNorm(B, variant, N)
            % Calculate normalization constant of a complex Bingham distribution.
            % There are different variants available.
            % It is still to be evaluated, which numerical variant is most robust.
            %
            % :param B: Hermitian Bingham parameter matrix or matrices.
            % :type B: D times D times K, where D is the dimension of the feature space.
            % :type variant: String
            %   (Default: 'SymbolicToolbox', Alternatives: 'MonteCarlo')
            % :return C: Normalization constant, dimension K times 1 column vector.
            %
            % If you choose the Monte Carlo variant, you need to provide the number of
            % samples as the last parameter N.
            %
            % The variant 'SymbolicToolbox' is from Kent1994Complex Equation (2.3).
            % The eigenvalue shift is taken from Kent1994Complex Page 286 as well.
            
            D = size(B, 1);
            K = size(B, 3);
            C = zeros(K, 1);
            
            if nargin == 1
                variant = 'SymbolicToolbox';
            end
            
            assert(isequal(B, permute(conj(B), [2 1 3:ndims(B)])), ...
                'B must be hermitian')
            
            switch variant
                case 'SymbolicToolbox'
                    for k = 1:K
                        [~, eigenvalues] = eig(B(:, :, k));
                        eigenvalues = diag(eigenvalues);
                        eigenvalue_shift = max(eigenvalues);
                        eigenvalues = sort(eigenvalues) - eigenvalue_shift;
                        eigenvalues = makeSureEigenvaluesAreNotTooClose(eigenvalues);
                        
                        if all(abs(eigenvalues) < 1e-3)
                            C(k) = log(2*pi^D / factorial(D-1));
                        else
                            functionHandle = ...
                                str2func(['cBinghamNorm', num2str(D)]);
                            C(k) = functionHandle(eigenvalues);
                        end
                        
                        % Correction term for eigenvalue shift according to
                        % Kent1994Complex page 286.
                        C(k) = log(C(k)) + eigenvalue_shift;
                    end
                case 'MonteCarlo'
                    volume = 2*pi^D / factorial(D-1);
                    for k = 1:K
                        %use monte carlo integration
                        Z = complex(randn(D, N), randn(D, N));
                        Z = bsxfun(@rdivide, Z, sqrt(sum(Z .* conj(Z), 1)));
                        p = exp(abs(dot(Z, B(:, :, k) * Z)));
                        
                        C(k) = log(sum(p)/N * volume);
                        %average value of pdf times surface area of unit sphere
                    end
            end
            
            if ~isreal(C(k)) && (C(k) > 0)
                error(['Unexpected value occured ', ...
                    'for complex Bingham normalization constant.']);
            end
            
            C = -C; % apply an inverted definition
            
            function lambda = makeSureEigenvaluesAreNotTooClose(lambda)
                % Workaround
                [lambda, indexMap] = sort(lambda, 'descend');
                b = diff(lambda);
                b = min(b, -1e-2 * ones(size(b)));
                
                % This index map is responsible to maintain the ordering of the input.
                % Thus, the sort process is compensated.
                lambda(indexMap) = cumsum([lambda(1); b]);
            end
        end
        
        function D = CauchySchwarzDivergence(cB1, cB2)
            % Calculate Cauchy Schwarz Divergence of two complex Bingham
            % distributions
            if isa(cB1, 'ComplexBinghamDistribution') ...
                    && isa(cB2, 'ComplexBinghamDistribution')
                B1 = cB1.B;
                B2 = cB2.B;
            elseif isa(cB1, 'double') ...
                    && isa(cB2, 'double')
                B1 = cB1;
                B2 = cB2;
            else
                error('wrong input')
            end
            
            assert(isequal(B1, permute(conj(B1), [2 1 3:ndims(B1)])), ...
                'B1 must be hermitian')
            
            assert(isequal(B2, permute(conj(B2), [2 1 3:ndims(B2)])), ...
                'B2 must be hermitian')
            
            log_c1 = ComplexBinghamDistribution.logNorm(2*B1);
            log_c2 = ComplexBinghamDistribution.logNorm(2*B2);
            log_c3 = ComplexBinghamDistribution.logNorm(B1 + B2);
            
            D = log_c3 - 1/2*(log_c1+log_c2);            
        end
    end
end