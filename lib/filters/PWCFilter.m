classdef PWCFilter < AbstractCircularFilter
    % A filter based on the piecewise constant distributions, basically a
    % Wonham filter.
    %
    % Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
    % Discrete Recursive Bayesian Filtering on Intervals and the Unit Circle
    % Proceedings of the 2016 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2016), 
    % Baden-Baden, Germany, September 2016.    
    
    properties
        pwc
    end
    
    methods
        function this = PWCFilter(n)
            % Constructor
            % Parameters:
            %   n (scalar)
            %       number of discretization steps
            assert(isscalar(n));
            this.setState(PWCDistribution(ones(1, n)/n));
        end
        
        function setState(this, dist)
            assert(isa(dist, 'AbstractCircularDistribution'));
            if isa(dist, 'PWCDistribution')
                this.pwc = dist;
            else
                this.pwc = PWCDistribution(PWCDistribution.calculateParametersNumerically(@dist.pdf, length(this.pwc.w)));
            end
        end
        
        function predict(this, A)
            % Prediction step based on transition matrix
            % Parameters:
            %   A (L x L matrix)
            %       system matrix
            w = this.pwc.w;
            w = A*w';
            w = w'/sum(w);
            this.pwc = PWCDistribution(w);
        end
        
        function update(this, H, z)
            % Measurement update based on measurement matrix
            % Parameters:
            %   H (Lw x L matrix)
            %       measurement matrix 
            %   z (scalar)
            %       measurement 
            %get right column from H
            assert(size(H,2) == length(this.pwc.w));
            assert(isscalar(z));
            Lw = size(H,1);
            row = floor(1 + z/2/pi*Lw);
            w = this.pwc.w;
            w = H(row,:).*w;
            this.pwc = PWCDistribution(w);
        end
        
        function updateLikelihood(this, likelihood, z)
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            
            assert(isa(likelihood,'function_handle'));
            L = length(this.pwc.w);
            tmp = zeros(1,L);
            for i=1:L
                tmp(i) = integral(@(x) likelihood(z,x), PWCDistribution.leftBorder(i,L), PWCDistribution.rightBorder(i,L));
            end
            this.pwc = PWCDistribution(tmp.*this.pwc.w);
        end
        
        function pwc = getEstimate(this)
            pwc = this.pwc;
        end
    end
    
    methods (Static)
        function A = calculateSystemMatrixNumerically(L, a, noiseDistribution)
            % Obtains system matrix by 2d numerical integration from system
            % function
            % Parameters:
            %   L (scalar)
            %       number of discretization intervals
            %   a (function handle)
            %       system fuction x_k+1 = a(x_k,w_k), needs to be
            %       vectorized
            %   noiseDistribution (AbstractCircularDistribution)
            %       noise (assumed to be defined on [0,2pi)
            % Returns;
            %   A (L x L matrix)
            %       system matrix
            assert(isscalar(L));
            assert(isa (noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(a,'function_handle'));

            A = zeros(L,L);
            for i=1:L
                l1 = PWCDistribution.leftBorder(i,L);
                r1 = PWCDistribution.rightBorder(i,L);
                for j=1:L
                    l2 = PWCDistribution.leftBorder(j,L);
                    r2 = PWCDistribution.rightBorder(j,L);
                    indicator = @(x) double(l2 < x & x < r2);
                    % todo this usually causes warnings
                    warning('off', 'MATLAB:integral2:maxFunEvalsFail');
                    A(j,i) = integral2(@(x,w) reshape(noiseDistribution.pdf(w(:)').*indicator(a(x(:)',w(:)')), size(x,1),size(x,2)), l1, r1, 0, 2*pi)*L/2/pi;
                    warning('on', 'MATLAB:integral2:maxFunEvalsFail');
                end
            end
        end
        
        function A = calculateSystemMatrixNumericallySampledNoise(L, a, noiseSamples, noiseWeights)
            % Obtains system matrix by 1d numerical integration from system
            % function
            % Parameters:
            %   L (scalar)
            %       number of discretization intervals
            %   a (function handle)
            %       system fuction x_k+1 = a(x_k,w_k), needs to be
            %       vectorized
            %   noiseSamples (1 x L matrix)
            %       L samples of the noise 
            %   noiseWeights (1 x L vector)
            %       weight of each sample
            % Returns;
            %   A (L x L matrix)
            %       system matrix
            assert(size(noiseWeights,1) == 1, 'weights most be row vector')
            assert(size(noiseSamples,2) == size(noiseWeights,2), 'samples and weights must match in size');
            assert(isa(a,'function_handle'));
            
            A = zeros(L,L);
            for i=1:L
                l1 = PWCDistribution.leftBorder(i,L);
                r1 = PWCDistribution.rightBorder(i,L);
                for j=1:L
                    l2 = PWCDistribution.leftBorder(j,L);
                    r2 = PWCDistribution.rightBorder(j,L);
                    indicator = @(x) double(l2 < x & x < r2);
                    A(j,i) = integral(@(x) sum(repmat(noiseWeights', 1, size(x,2)).*indicator(a(repmat(x,size(noiseSamples,2),1),repmat(noiseSamples', 1, size(x,2))))), l1, r1)*L/2/pi;
                end
            end
        end
        
        function H = calculateMeasurementMatrixNumerically(L, Lmeas, h, noiseDistribution)
            % Obtains system matrix by 2d numerical integration from
            % measurement function
            % Parameters:
            %   L (scalar)
            %       number of discretization intervals for state
            %   Lmeas (scalar)
            %       number of discretization intervals for measurement
            %   h (function handle)
            %       system fuction z_k = h(x_k,v_k), needs to be
            %       vectorized
            %   noiseDistribution (AbstractCircularDistribution)
            %       noise (assumed to be defined on [0,2pi)
            % Returns;
            %   H (Lmeas x L matrix)
            %       measurement matrix            
            assert(isscalar(L));
            assert(isscalar(Lmeas));
            assert(isa (noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(h,'function_handle'));
            
            H = zeros(Lmeas,L);
            for i=1:Lmeas
                l1 = PWCDistribution.leftBorder(i,Lmeas);
                r1 = PWCDistribution.rightBorder(i,Lmeas);
                for j=1:L
                    l2 = PWCDistribution.leftBorder(j,L);
                    r2 = PWCDistribution.rightBorder(j,L);
                    indicator = @(x) double(l1 < x & x < r1);
                    % todo this usually causes warnings
                    warning('off', 'MATLAB:integral2:maxFunEvalsFail');
                    H(i,j) = integral2(@(x,v) reshape(noiseDistribution.pdf(v(:)').*indicator(h(x(:)',v(:)')), size(x,1),size(x,2)), l2, r2, 0, 2*pi)*L/2/pi;
                    warning('on', 'MATLAB:integral2:maxFunEvalsFail');
                end
            end
        end
        
        function H = calculateMeasurementMatrixNumericallySampledNoise(L, Lmeas, h, noiseSamples, noiseWeights)
            % Obtains system matrix by 1d numerical integration from
            % measurement function
            % Parameters:
            %   L (scalar)
            %       number of discretization intervals for state
            %   Lmeas (scalar)
            %       number of discretization intervals for measurement
            %   h (function handle)
            %       system fuction z_k = h(x_k,v_k), needs to be
            %       vectorized
            %   noiseSamples (1 x L matrix)
            %       L samples of the noise 
            %   noiseWeights (1 x L vector)
            %       weight of each sample
            % Returns;
            %   H (Lmeas x L matrix)
            %       measurement matrix            
            assert(isscalar(L));
            assert(isscalar(Lmeas));
            assert(isa(h,'function_handle'));
            assert(size(noiseWeights,1) == 1, 'weights most be row vector')
            assert(size(noiseSamples,2) == size(noiseWeights,2), 'samples and weights must match in size');
            
            H = zeros(Lmeas,L);
            for i=1:Lmeas
                l1 = PWCDistribution.leftBorder(i,Lmeas);
                r1 = PWCDistribution.rightBorder(i,Lmeas);
                for j=1:L
                    l2 = PWCDistribution.leftBorder(j,L);
                    r2 = PWCDistribution.rightBorder(j,L);
                    indicator = @(x) double(l1 < x & x < r1);
                    H(i,j) = integral(@(x) sum(repmat(noiseWeights', 1, size(x,2)).*indicator(h(repmat(x,size(noiseSamples,2),1),repmat(noiseSamples', 1, size(x,2))))), l2, r2)*L/2/pi;
                end
            end
        end
        
        function H = calculateMeasurementMatrixNumericallyFromLikelihood(L, Lmeas, likelihood)
            % Obtains measurement matrix by 2d numerical integration
            % Parameters:
            %   L (scalar)
            %       number of discretization intervals for state
            %   Lmeas (scalar)
            %       number of discretization intervals for measurement
            %   likelihood (function handle)
            %       likelihood fuction likelihood(z,x) = f(z|x)
            % Returns;
            %   H (Lmeas x L matrix)
            %       measurement matrix     
            H = zeros(Lmeas,L);
            for i=1:Lmeas
                l1 = PWCDistribution.leftBorder(i,Lmeas);
                r1 = PWCDistribution.rightBorder(i,Lmeas);
                for j=1:L
                    l2 = PWCDistribution.leftBorder(j,L);
                    r2 = PWCDistribution.rightBorder(j,L);
                    H(i,j) = integral2(@(z,x) reshape(likelihood(z(:)',x(:)'), size(x,1), size(x,2)), l1, r1, l2, r2)*L/2/pi;
                end
            end
        end
    end
    
end

