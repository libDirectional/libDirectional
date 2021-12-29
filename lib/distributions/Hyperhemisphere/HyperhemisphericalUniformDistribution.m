classdef HyperhemisphericalUniformDistribution < AbstractHyperhemisphericalDistribution & AbstractUniformDistribution
    methods
        function this = HyperhemisphericalUniformDistribution(dim_)
            arguments
                dim_ (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(dim_, 2)}
            end
            this.dim = dim_;
        end
        
        function s = sample(this, n)
            % Stocahastic sampling
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column     
            %
            arguments
                this (1,1) HyperhemisphericalUniformDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            s = HypersphericalUniformDistribution(this.dim).sample(n);
            % Mirror ones with negative on the last dimension up for hemisphere. This
            % may give a disadvantage to ones with exactly zero at the first dimension but
            % since this is hit with quasi probability zero, we neglect
            % this.
            s = (1-2*(s(end,:)<0)).*s;
        end
    end 
end