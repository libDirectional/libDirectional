classdef HyperhemisphericalBinghamDistribution < AbstractHyperhemisphericalDistribution    % Represents a Watson distribution.
    
    properties
        distFullSphere BinghamDistribution
    end
    
    methods
        function this = HyperhemisphericalBinghamDistribution(Z_, M_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            arguments
                Z_ (:,1) double
                M_ (:,:) double
            end
            this.distFullSphere = BinghamDistribution(Z_, M_);
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
                this (1,1) HyperhemisphericalBinghamDistribution
                xa (:,:) double
            end
            p = 2*this.distFullSphere.pdf(xa);
        end   

        function ZBingham = Z(this)
            ZBingham = this.distFullSphere.Z;
        end

        function MBingham = M(this)
            MBingham = this.distFullSphere.M;
        end

        function FBingham = F(this)
            FBingham = this.distFullSphere.F;
        end

        function dFBingham = dF(this)
            dFBingham = this.distFullSphere.dF;
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
                this (1,1) HyperhemisphericalBinghamDistribution
                n (1,1) double {mustBeInteger, mustBePositive}
            end
            sFull = this.distFullSphere.sample(n);
            s = sFull.*(-1).^(sFull(end,:)<0); % Mirror to upper hemisphere
        end
        
        function B = multiply(this, B2)
            arguments
                this (1,1) HyperhemisphericalBinghamDistribution
                B2 (1,1) HyperhemisphericalBinghamDistribution
            end
            B = this;
            B.distFullSphere = this.distFullSphere.multiply(B2.distFullSphere);
        end

        function B = compose(this, B2)
            arguments
                this (1,1) HyperhemisphericalBinghamDistribution
                B2 (1,1) HyperhemisphericalBinghamDistribution
            end
            B = this;
            B.distFullSphere = this.distFullSphere.compose(B2.distFullSphere);
        end

        function m = mode(this)
            arguments
                this (1,1) HyperhemisphericalBinghamDistribution
            end
            m = this.distFullSphere.mode;
        end

        function m = meanAxis(this)
            arguments
                this (1,1) HyperhemisphericalBinghamDistribution
            end
            m = this.distFullSphere.meanAxis;
        end
    end
end