classdef WatsonDistribution < AbstractHypersphericalDistribution
    % Represents a Watson distribution.
    %
    % Watson, G. S. 
    % Equatorial Distributions on a Sphere
    % Biometrika, 1965, 52, 193-201
    
    
    properties
        kappa   % concentration (scalar)
        mu      % mean as vector
        C       % normalization constant
    end
    
    methods
        function W = WatsonDistribution(mu_, kappa_)
            %% Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            epsilon = 1E-6;
            assert(size(mu_,2) == 1, 'mu must be a row vector');
            assert(abs(norm(mu_) - 1)<epsilon, 'mu must be a normalized');
            assert(isscalar(kappa_));
            
            W.mu = mu_;
            W.kappa = kappa_;
            
            W.dim = size(mu_,1);
            W.C = gamma(W.dim/2)/2/pi^(W.dim/2)/hypergeom(1/2, W.dim/2, W.kappa);
        end
        
        function p = pdf(this, xa)
            % Evaluates pdf at each column of xa
            % Parameters:
            %   xa (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of xa
            p = zeros(1, size(xa,2));
            for i=1:size(xa,2)
                x = xa(:,i);
                p(i) = this.C * exp( this.kappa * (this.mu' * x).^2);
            end
        end      
        
        function B = toBingham(this)
            % Converts this distribution to a Bingham distribution. The
            % conversion is exact because a Watson distribution is an
            % isotropic Bingham distribution.
            %
            % Returns:
            %   B (BinghamDstribution)
            %       the Bingham distribution with identical pdf
            if this.kappa<0
                error ('conversion to Bingham is not implemented yet for kappa<0');
            end
            
            M = repmat(this.mu, 1, this.dim);
            % make vectors in M linear independent
            E = eye(this.dim, this.dim); 
            E(1,1)=0;
            M = M + E;
            % make vectors in M orthogonal
            [Q,~] = qr(M);
            M = [Q(:,2:end)  Q(:,1)]; % put mu at in the last column
            Z = [repmat(-this.kappa,this.dim-1,1); 0];
            B = BinghamDistribution(Z,M);
        end
        
        function s = sample(this,n)
            % Generate samples from Watson distribution.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            s = this.toBingham.sample(n); % use Bingham sampling
        end
        
        function m = mode(this)
            % Calculate the mode of a Watson distribution
            % Returns:
            %   m (column vector)
            %       mode of the distribution (note that -m is the mode as well)
            if this.kappa>=0
                m = this.mu; %todo: this is only correct for kappa>=0
            else
                error ('mode for kappa<0 is not implemented yet');
            end
        end
    end
end