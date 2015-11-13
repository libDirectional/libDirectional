classdef GeneralToroidalMixture < GeneralHypertoroidalMixture & AbstractToroidalDistribution
    methods
        function this=GeneralToroidalMixture(hds, w)
            assert(isa(hds,'cell') && all(cellfun(@(cd)isa(cd,'AbstractToroidalDistribution'),hds)),...
                    'hds must be a cell array of toroidal distributions');
            this@GeneralHypertoroidalMixture(hds, w);
        end
    end
end