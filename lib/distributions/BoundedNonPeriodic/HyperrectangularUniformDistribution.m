classdef HyperrectangularUniformDistribution < AbstractHyperrectangularDistribution & AbstractUniformDistribution
    methods
        function this = HyperrectangularUniformDistribution(bounds)
            arguments
                bounds (2,:) double {mustBeNonempty}
            end
            this.bounds = bounds;
        end
    end
end