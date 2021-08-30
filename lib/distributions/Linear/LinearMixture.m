classdef LinearMixture < AbstractMixture & AbstractLinearDistribution
    methods
        function this = LinearMixture(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractLinearDistribution'),dists)),...
                'dists must be a cell array of linear distributions');
            if all(cellfun(@(dist)isa(dist,'GaussianDistribution'),dists))
                warning('LinearMixture:AllGaussians','For mixtures of Gaussians, consider using GaussianMixture.');
            end
            this = this@AbstractMixture(dists, w);
        end 
    end
end