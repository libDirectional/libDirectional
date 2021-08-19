classdef SE2WNDistribution < HypercylindricalWNDistribution & AbstractSE2Distribution
    methods
        function this = SE2WNDistribution(mu_, C_)
            arguments
                mu_ (:,1) double
                C_ (:,:) double
            end
            this@HypercylindricalWNDistribution(mu_, C_, 1);
        end
    end
end