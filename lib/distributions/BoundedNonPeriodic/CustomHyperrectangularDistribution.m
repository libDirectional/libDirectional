classdef CustomHyperrectangularDistribution < AbstractHyperrectangularDistribution & CustomDistribution
    methods
        function this = CustomHyperrectangularDistribution(f, bounds)
            dim = size(bounds,2);
            this@CustomDistribution(f,dim);
            this.bounds = bounds;
        end
    end
end