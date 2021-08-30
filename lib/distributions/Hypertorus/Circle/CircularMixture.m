classdef CircularMixture < AbstractCircularDistribution & HypertoroidalMixture
    % Mixture of multiple circular distributions. The distributions may belong
    % to different classes.
    
    methods
        function this = CircularMixture(dists, w)
            % Constructor
            this@HypertoroidalMixture(dists, w);
            assert(isa(dists,'cell') && all(cellfun(@(cd)isa(cd,'AbstractCircularDistribution'),dists)),...
                'dists must be a cell array of circular distributions');
            assert(all(size(dists) == size(w)),'size of dists and w must be equal');
            
            if all(cellfun(@(cd)isa(cd,'FourierDistribution'),dists))
                warning('Mixtures of Fourier distributions can be built by combining the Fourier coefficients so using a mixture may not be necessary');
            end
            if all(cellfun(@(cd)isa(cd,'WDDistribution'),dists))
                warning('Mixtures of WDDistributions can usually be combined into one WDDistribution.');
            end
            this.dists = dists;
            this.w = w/sum(w);
        end
    end
end

