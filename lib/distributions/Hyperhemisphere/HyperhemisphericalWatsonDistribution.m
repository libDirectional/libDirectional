classdef HyperhemisphericalWatsonDistribution < AbstractHyperhemisphericalDistribution
    % Watson distribution for the hyperhemisphere
    
    properties
        distFullSphere WatsonDistribution
    end
    
    methods
        function this = HyperhemisphericalWatsonDistribution(mu_, kappa_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            arguments
                mu_ (:,1) double
                kappa_ (1,1) double
            end
            assert(mu_(end)>=0)
            this.distFullSphere = WatsonDistribution(mu_, kappa_);
            this.dim = this.distFullSphere.dim;
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
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
                xa (:,:) double
            end
            p = 2*this.distFullSphere.pdf(xa);
        end   

        function w = setMode(this, mu)
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
                mu (:,1) double
            end
            w = this;
            w.distFullSphere.mu = mu;
        end
        function s = sample(this,n)
            % Generate samples from Watson distribution.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   s (d x n matrix)
            %       generated samples (one sample per column)
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
                n (1,1) double {mustBeInteger, mustBePositive}
            end
            sFull = this.distFullSphere.sample(n);
            s = sFull.*(-1).^(sFull(end,:)<0); % Mirror to upper hemisphere
        end
        
        function m = mu(this)
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
            end
            m = this.distFullSphere.mu;
        end

        function k = kappa(this)
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
            end
            k = this.distFullSphere.kappa;
        end

        function m = mode(this)
            % Calculate the mode of a Watson distribution
            % Returns:
            %   m (d x 1 column vector)
            %       mode of the distribution
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
            end
            m = this.distFullSphere.mu;
        end

        function distShifted = shift(this, offsets)
            % There is no true shifting for the hyperhemisphere. This is a function for compatibility and only works when mu is [0,0,...,1].
            arguments
                this (1,1) HyperhemisphericalWatsonDistribution
                offsets (:,1) double {mustBeNonempty}
            end
            assert(isequal(this.distFullSphere.mu,[zeros(this.dim-1,1);1]), 'There is no true shifting for the hypersphere. This is a function for compatibility and only works when mu is [0,0,...,1].');
            distShifted = this;
            distShifted.distFullSphere.mu = offsets;
        end
    end
end