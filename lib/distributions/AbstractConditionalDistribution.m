classdef (Abstract) AbstractConditionalDistribution < AbstractDistribution
    methods
        function m = mean(~) %#ok<STOUT> 
            error('Mean currently not supported.')
        end
    end
end