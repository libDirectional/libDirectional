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
end

