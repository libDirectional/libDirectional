classdef HypersphericalUniformDistribution < AbstractHypersphericalDistribution & AbstractUniformDistribution
    methods
        function this = HypersphericalUniformDistribution(dim_)
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
            % General algorithm based on "A note on a method for generating 
            % points uniformly on n-dimensional spheres" by Mervin E. Muller April 1959.
            % Algorithm for 2-sphere based on "Spherical sampling by archimedes' theorem"
            % by Min-Zhi Shao and Norman Badler, 1996
            arguments
                this (1,1) HypersphericalUniformDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            
            if this.dim == 3
                s = NaN(this.dim,n);
                phi = 2*pi*rand(1,n);
                s(3,:) = rand(n,1)*2-1;
                r = sqrt(1-s(3,:).^2);
                s(1,:) = r.*cos(phi);
                s(2,:) = r.*sin(phi);
            else
                samplesUnnorm = randn(this.dim,n);
                s = samplesUnnorm./vecnorm(samplesUnnorm);
            end
        end
    end 
end