classdef GaussianDistribution
    % Gaussian distribution (possibly multidimensional).
    
    properties
        mu  % mean
        C   % covariance, C=sigma^2 in 1D
    end
    
    methods
        function this = GaussianDistribution(mu_, C_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (mean vector)
            %   C_ (d x d)
            %       covariance matrix
            assert(size(mu_,2)==1, 'mu must have 1 column')
            assert(size(mu_,1)==size(C_,1), 'size of C invalid')
            assert(size(mu_,1)==size(C_,2), 'size of C invalid')
            this.mu = mu_;
            
            % check that C is positive definite
            dim = length(mu_);
            if dim==1
                assert(C_>0, 'C must be positive definite');
            elseif dim==2
                assert(C_(1,1)>0 && det(C_)>0, 'C must be positive definite');
            else
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
            assert(size(xa,1)==size(this.mu,1), 'dimension incorrect')
            
            p = mvnpdffast(xa',this.mu',this.C)';
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
            assert(length(this.C) == 1, 'conversion to WN is only possible for 1D Gaussians');
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
            samples = mvnrnd(this.mu, this.C, n)';
        end
        
    end
    
end

