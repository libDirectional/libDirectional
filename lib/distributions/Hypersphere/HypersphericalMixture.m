classdef HypersphericalMixture < AbstractHypersphericalDistribution & AbstractMixture
    % Mixture of multiple hyperspherical distributions. The distributions may belong
    % to different classes.
    
    methods
        function this = HypersphericalMixture(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractHypersphericalDistribution'),dists)),...
                'dists must be a cell array of hyperspherical distributions');
            if all(cellfun(@(dist)isa(dist,'AbstractSphericalHarmonicDistribution'),dists))
                warning('Creating a mixture of Spherical Harmonics may not be necessary.');
            end
            this = this@AbstractMixture(dists, w);
        end
    end
    methods (Sealed)
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       n samples on the dim-dimensional hypersphere
            
            % Sample component first, then sample from the chosen component
            d = discretesample(this.w,n);
            s = zeros(this.dim,n);
            occurrences=histc(d,1:length(this.dists));
            count=1;
            for i=1:length(this.dists)
                s(:,count:count+occurrences(i)-1) = this.dists{i}.sample(occurrences(i));
                count=count+occurrences(i);
            end
            [~,order]=sort(d);
            s(:,order)=s;
        end
    end
end

