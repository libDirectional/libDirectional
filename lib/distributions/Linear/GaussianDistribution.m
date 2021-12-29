classdef GaussianDistribution < AbstractLinearDistribution
    % Gaussian distribution (possibly multidimensional).
    
    properties
        mu (:,1) double % mean
        C (:,:) double % covariance, C=sigma^2 in 1D
    end
    
    methods
        function this = GaussianDistribution(mu_, C_, checkValidity)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (mean vector)
            %   C_ (d x d)
            %       covariance matrix
            arguments
                mu_ (:,1) double {mustBeNonempty}
                C_ (:,:) double {mustBeNonempty}
                checkValidity (1,1) logical = true % Can disable because it is not a cheap operation
            end
            assert(size(mu_,1)==size(C_,1), 'size of C invalid')
            assert(size(mu_,1)==size(C_,2), 'size of C invalid')
            this.dim = size(mu_,1);
            this.mu = mu_;
            
            % check that C is positive definite
            dim = length(mu_);
            if checkValidity&&dim==1
                assert(C_>0, 'C must be positive definite');
            elseif checkValidity&&dim==2
                assert(C_(1,1)>0 && det(C_)>0, 'C must be positive definite');
            elseif checkValidity
                chol(C_); % will fail if C_ is not pos. def.
            end
            
            this.C = C_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (d x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            arguments
                this (1,1) GaussianDistribution
                xa double {mustBeNonempty}
            end
            assert(size(xa,1)==size(this.mu,1), 'dimension incorrect')
            
            p = mvnpdf(xa',this.mu',this.C)'; % just for me...
        end
        
        function g = shift(this, offsets)
            arguments
                this (1,1) GaussianDistribution
                offsets (:,1) double
            end
            assert(numel(offsets)==this.dim);
            g = this;
            g.mu = this.mu + offsets;
        end
        
        function mu = mean(this)
            arguments
                this (1,1) GaussianDistribution
            end
            mu = this.mu;
        end
        
        function m = mode(this)
            % Determines the mode of the distribution, i.e., the point
            % where the pdf is largest.
            %
            % Returns:
            %   m (vector)
            %       the mode
            arguments
                this (1,1) GaussianDistribution
            end
            m = this.mu;
        end

        function dist = setMode(this, newMode)
            arguments
                this (1,1) GaussianDistribution
                newMode (:,1) double
            end
            dist = this;
            dist.mu = newMode;
        end

        function C = covariance(this)
            arguments
                this (1,1) GaussianDistribution
            end
            C = this.C;
        end
        
        function dist = multiply(this, other)
            arguments
                this (1,1) GaussianDistribution
                other (1,1) GaussianDistribution
            end
            dist = this; % Skip constructor (attention, due to numerical issues, this may lead to non spd matrices
            K = this.C / (this.C + other.C);
            dist.mu = this.mu + K * (other.mu - this.mu);
            dist.C = this.C - K*this.C;
        end
        
        function dist = marginalizeOut(this,dimensions)
            arguments
                this (1,1) GaussianDistribution
                dimensions (1,:) double {mustBeInteger,mustBePositive}
            end
            assert(all(dimensions<=this.dim));
            dimsSorted = unique(dimensions); % Mainly for next assertion
            assert(numel(dimsSorted)==numel(dimensions));
            remaingDims = 1:this.dim;
            remaingDims(dimsSorted) = [];
            dist = this;
            dist.mu = dist.mu(remaingDims);
            dist.C = dist.C(remaingDims, remaingDims);
            dist.dim = numel(remaingDims);
        end

        function wn = toWN(this)
            % Convert to WN (only for 1D Gaussians)
            %
            % Returns:
            %   wn (WNDistribution)
            %       WN with same parameters
            %
            % this is a simple conversion that just keeps the parameters
            % for large sigma, better conversions are possible
            arguments
                this (1,1) GaussianDistribution
            end
            assert(this.dim == 1, 'conversion to WN is only possible for 1D Gaussians');
            wn = WNDistribution(this.mu, sqrt(this.C));
        end
        
        function samples = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (d x n)
            %       one sample per column
            arguments
                this (1,1) GaussianDistribution
                n (1,1) {mustBeInteger,mustBePositive}
            end
            samples = mvnrnd(this.mu, this.C, n)';
        end
        
        function [s,w] = sampleDeterministic(this)
            % Obtain UKF samples from the distribution.
            %
            % Returns:
            %   s (d x n)
            %       one sample per column
            %   w (1 x n)
            %       weights (always uniform for now)
            arguments
                this (1,1) GaussianDistribution
            end
            g = GaussianSamplingUKF;
            [s,w,numGaussianSamples] = g.getSamples(Gaussian(this.mu, this.C));
            if isscalar(w)
                w = repmat(w, 1, numGaussianSamples);
            end
        end
        
        function [s,w] = sampleDeterministicS2KF(this, n)
            % Obtain S2KF samples from the distribution (using LCD).
            %
            % Returns:
            %   s (d x n)
            %       one sample per column
            %   w (1 x n)
            %       weights (always uniform for now)
            arguments
                this (1,1) GaussianDistribution
                n (1,1) {mustBeInteger,mustBePositive}
            end
            g = GaussianSamplingLCD;
            g.setNumSamples(n);
            [s,w,numGaussianSamples] = g.getSamples(Gaussian(this.mu, this.C));
            if isscalar(w)
                w = repmat(w, 1, numGaussianSamples);
            end
        end        
    end
    
    methods (Static)
        function gaussian = fromDistribution(dist)
            arguments
                dist (1,1) AbstractLinearDistribution
            end
            if isa(dist,'GaussianMixtureDistribution')
                gaussian = dist.toGaussian;
            else
                gaussian = GaussianDistribution(dist.mean, dist.covariance);
            end
        end
    end
end

