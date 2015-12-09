classdef ToroidalMixture < HypertoroidalMixture & AbstractToroidalDistribution
    % Mixture of multiple toroidal distributions. The distributions may belong
    % to different classes.
    methods
        function this = ToroidalMixture(hds, w)
            assert(isa(hds,'cell') && all(cellfun(@(cd)isa(cd,'AbstractToroidalDistribution'),hds)),...
                    'hds must be a cell array of toroidal distributions');
            this@HypertoroidalMixture(hds, w);
        end
    end
end